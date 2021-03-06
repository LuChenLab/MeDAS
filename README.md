# MeDAS: a Metazoan Developmental Alternative Splicing database 
This repository contains the pipeline for generating exonic PSI used in MeDAS.

## MeDAS URL
[MeDAS at https://das.chenlulab.com](https://das.chenlulab.com)

## Pipeline overview
<center>
<img src="./doc/MeDAS_pipeline.png" width="75%" height="75%" />
</center>

## Retrieve data from SRA
For each dataset (project), a table named "SraRunInfo.csv" was downloaded from SRA.
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
Trimmomatic
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
```

## Align
STAR
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
RSeQC
```bash
BED12="path_to_bed12" ## bed12 file of corresponding annotation

infer_experiment.py \
    -s 2000000 \
    -i ${STAROUT}/${SRA}.Aligned.sortedByCoord.out.bam \
    -r ${BED12}
```

## Quantify
RSEM
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
```R
## R
## exclude low reliable SJs
source("./Rscripts/filter_SJs.R")
raw_SJs_path <- "path_to_HomSap" ## STAR output directory for a species, STAROUT
filtered_SJs <- "HomSap_filtered_SJs.tsv" ## SJs for calculate PSI

res <- 
    multi_merge(
        path = raw_SJs_path,
        pattern = ".SJ.out.tab",
        minSJ = 10,
        minSJs = 100,
        minSamps = 2,
        uniqueMapOnly = TRUE,
        SJtype = "allAnnotatedAndCanonicallyNovel",
        cores = 20
    )
  
fwrite(res, filtered_SJs, sep = "\t", row.names = F, col.names = F)
```
```bash
## bash
## generate new SJ.out.tab of filtered SJs 
## for HomSap
FILTEREDSJ="path_to_filtered_SJ_tabs"
grep -F -f HomSap_filtered_SJs.tsv ${STAROUT}/${SRA}.SJ.out.tab > ${FILTEREDSJ}/${SRA}.SJ.out.tab
```

## Calculate PSI of exonic parts
DEXseq, bedtools
```bash
## prepare exonic part gff
python2 path/to/dexseq_prepare_annotation.py path/to/HomSap.gtf HomSap.reduce.gtf

awk '{OFS="\t"}{if ($3 == "exonic_part") print $1,$2,$3,$4,$5,$6,$7,$8,$14":"$12}' HomSap.reduce.gtf | \
    sed 's=[";]==g' | \
    sort -k1,1 -k2,2n > HomSap_Exonic_part.gff
```
```bash
## call PSI
LEN="mapped_reads_length_from_STAR" ## integer

bash path/to/ExonicPartPSI_2.sh \
    path/to/bedtools2.23/bedtools \
    path/to/HomSap_Exonic_part.gff \
    ${SRA}.Aligned.toTranscriptome.out.bam \
    ${LEN} \
    ${FILTEREDSJ}/${SRA}.SJ.out.tab \
    ${SRA}

## for Monodelphis domestica (large chromosome size, over 500Mb), needs genome file
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
```bash
## merge PSI outputs
Rscript ExonicPart_PSI/mergePSI.R -h
```

## AS type  
SUPPA2
```bash
GTF="gtf_file"
ABBR="species_abbreviation"
suppa.py generateEvents \
  -i ${GTF} -o \
  ${ABBR} \
  -f ioe \
  -e SE SS MX RI FL
```

## Call time course expression  
edgeR, maSigPro
```r
## R
## For example
source("./Rscripts/cpm_maSigPro.R")

gene_res_musmus <- 
    run_maSigPro(
        Species = "MusMus",
        Meta = "all_Meta.tsv",
        Outpath = "./",
        Count = "MusMus_Gene_expected_count.tsv",
        TPM = "MusMus_Gene_TPM.tsv"
    )
```

## Call correlation and KW test of PSI
```r
## R
source("./Rscripts/call_cor_tau_kw.R")

psi_res_musmus <-
    sum_foo_psi(
        mat,    ## matrix or df, with rownames (ExonicPartName) and colnames (sample ID)
        colDat, ## metainfo, with rownames (sample ID), for all samples of specific tissue
        col,    ## column name of stage-culumn, ordered factor (rank)
        keep = 0.3 ## maximum fraction of NAs to keep when call correlation
    )
```
