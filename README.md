## Retrieve data from SRA
#### For each dataset (project), a table named "SraRunInfo.csv" were downloaded from SRA.
#### 
```bash
RAWPATH="path_to_raw_fg.gz" ## directory for raw data, according to ENV

cd ${RAWPATH}

## download
awk -F "," '/http/{print $10}' SraRunInfo.csv | \
    xargs -i -P 10 bash -c "wget -c {}"

## decrypt
ls | grep "[EDSZ]RR" | \
    xargs -i bash -c "echo {}; \
    fasterq-dump --split-3 -e 10 -f -O ./ ./{}"

## compress
ls *fastq | xargs -i -P 10 bash -c "echo {};gzip {}"
```

## QC
#### quality control, trimmomatic
```bash
SOFTPATH="path_to_trimmomatic-0.38" ## according to ENV
CLEANPATH="path_to_clean_fg.gz" ## according to ENV
SRA="sampleID_or_runID"

## Illumina GA, PAIRED END
java -jar ${SOFTPATH}/trimmomatic-0.38.jar PE -threads 4 \
    -basein ${RAWPATH}/${SRA}_1.fastq.gz \
    -baseout ${CLEANPATH}/${SRA}.trimmed.fastq.gz \
    ILLUMINACLIP:${SOFTPATH}/adapters/TruSeq2-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

## Illumina GA, SINGLE END
java -jar ${SOFTPATH}/trimmomatic-0.38.jar SE -threads 4 \
    ${RAWPATH}/${SRA}.fastq.gz \
    ${CLEANPATH}/${SRA}.trimmed.fastq.gz \
    ILLUMINACLIP:${SOFTPATH}/adapters/TruSeq2-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

## Illumina HiSeq, PAIRED END
java -jar ${SOFTPATH}/trimmomatic-0.38.jar PE -threads 4 \
    -basein ${RAWPATH}/${SRA}_1.fastq.gz \
    -baseout ${CLEANPATH}/${SRA}.trimmed.fastq.gz \
    ILLUMINACLIP:${SOFTPATH}/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

## Illumina HiSeq, SINGLE END
java -jar ${SOFTPATH}/trimmomatic-0.38.jar SE -threads 4 \
    ${RAWPATH}/${SRA}.fastq.gz \
    ${CLEANPATH}/${SRA}.trimmed.fastq.gz \
    ILLUMINACLIP:${SOFTPATH}/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

## other platform, PAIRED END
fastp \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --length_required 25 \
    --detect_adapter_for_pe \
    --cut_by_quality5 3 \
    --cut_by_quality3 3 \
    --cut_window_size 4 \
    --cut_mean_quality 15 \
    -i ${RAWPATH}/${SRA}_1.fastq.gz \
    -I ${RAWPATH}/${SRA}_2.fastq.gz \
    -o ${CLEANPATH}/${SRA}.trimmed_1P.fastq.gz \
    -O ${CLEANPATH}/${SRA}.trimmed_2P.fastq.gz

## other platform, SINGLE END
fastp \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --length_required 25 \
    --cut_by_quality5 3 \
    --cut_by_quality3 3 \
    --cut_window_size 4 \
    --cut_mean_quality 15 \
    -i ${RAWPATH}/${SRA}.fastq.gz \
    -o ${CLEANPATH}/${SRA}.trimmed.fastq.gz
```

## Align
```bash
STARINDEX="path_to_STAR_index" ## according to ENV
STAROUT="path_to_star_outdir" ## according to ENV

## SINLE END
STAR --runThreadN 8 \
    --genomeDir ${STARINDEX} \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMcompression 9 \
    --limitBAMsortRAM 100000000000 \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --readFilesIn ${CLEANPATH}/${SRA}.trimmed.fastq.gz \
    --outFileNamePrefix ${STAROUT}/${SRA}.

## PAIRED END
STAR --runThreadN 8 \
    --genomeDir ${STARINDEX} \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMcompression 9 \
    --limitBAMsortRAM 100000000000 \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --readFilesIn \
        ${CLEANPATH}/${SRA}.trimmed_1P.fastq.gz \
        ${CLEANPATH}/${SRA}.trimmed_2P.fastq.gz \
    --outFileNamePrefix ${STAROUT}/${SRA}.
```

## Infer experiment
```bash
BED12="path_to_bed12" ## bed12 file of corresponding annotation

infer_experiment.py \
    -s 2000000 \
    -i ${STAROUT}/${SRA}.Aligned.sortedByCoord.out.bam \
    -r ${BED12}
```

## Quantify
```bash
RSEMINDEX="path_to_RSEM_index" ## index of RSEM
STRAND="none_reverse_or_forward" ## strand-specificity of library, accroding to result of infer experiment
RSEMOUT="path_to_RSEM_outdir"

## PAIRED END
rsem-calculate-expression \
    --paired-end \
    --strandedness ${STRAND} \
    --alignments -p 8 \
    --no-bam-output \
    --quiet \
    ${STAROUT}/${SRA}.Aligned.toTranscriptome.out.bam \
    ${RSEMINDEX} \
    ${RSEMOUT}/${SRA}.

## SINGLE END
rsem-calculate-expression \
    --paired-end \
    --strandedness ${STRAND} \
    --alignments -p 8 \
    --no-bam-output \
    --quiet \
    ${STAROUT}/${SRA}.Aligned.toTranscriptome.out.bam \
    ${RSEMINDEX} \
    ${RSEMOUT}/${SRA}.
```

## Filter SJ
#### exclude low reliable SJs
```R
## R
source("./Rscripts/filter_SJs.R")
raw_SJs_path <- "path_to_HomSap" ## STAR output directory for a species, STAROUT
filtered_SJs <- "HomSap_filtered_SJs.tsv" ## SJs for calculate PSI

res <- 
    multi_merge(
        path = ,
        pattern = raw_SJs_path,
        minSJ = 10,
        minSJs = 100,
        minSamps = 2,
        uniqueMapOnly = TRUE,
        SJtype = "allAnnotatedAndCanonicallyNovel",
        cores = 20
    )
  
fwrite(res, filtered_SJs, sep = "\t", row.names = F, col.names = F)
```
#### generate new SJ.out.tab of filtered SJs 
```bash
## for HomSap
FILTEREDSJ="path_to_filtered_SJ_tabs"
grep -F -f HomSap_filtered_SJs.tsv ${STAROUT}/${SRA}.SJ.out.tab > ${FILTEREDSJ}/${SRA}.SJ.out.tab
```

## Calculate PSI of exonic part
#### prepare exonic part gff
```bash
python2 path/to/dexseq_prepare_annotation.py path/to/HomSap.gtf HomSap.reduce.gtf

awk '{OFS="\t"}{if ($3 == "exonic_part") print $1,$2,$3,$4,$5,$6,$7,$8,$14":"$12}' HomSap.reduce.gtf | \
    sed 's=[";]==g' | \
    sort -k1,1 -k2,2n > HomSap_Exonic_part.gff
```

#### call PSI
```bash
LEN="mapped_reads_length_of_STAR" ## integer

bash path/to/ExonicPartPSI_2.sh \
    path/to/bedtools2.23/bedtools \
    path/to/HomSap_Exonic_part.gff \
    ${SRA}.Aligned.toTranscriptome.out.bam \
    ${LEN} \
    ${FILTEREDSJ}/${SRA}.SJ.out.tab \
    ${SRA}

## for Monodelphis domestica (large chromsome size, over 500Mb), needs genome file
awk '{print $1"\t"$2}' MonDom.fa.fai | sort -k1,1 -k2,2n > MonDom.genome

bash path/to/ExonicPartPSI_2.sh \
    path/to/bedtools2.23/bedtools \
    path/to/MonDom_Exonic_part.gff \
    MonDom.genome \
    ${SRA}.Aligned.toTranscriptome.out.bam \
    ${LEN} \
    ${FILTEREDSJ}/${SRA}.SJ.out.tab \
    ${SRA}
```
