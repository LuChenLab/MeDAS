## Exonic Part PSI
#### Reference: [doi: 10.1002/0471142905.hg1116s87](./doc/ExonicPart_PSI.pdf)
- Alternative Splicing Signatures in RNA-seq Data: Percent Spliced in (PSI), Curr. Protoc. Hum. Genet.

##### PSI.sh derived from 
- https://github.com/lehner-lab/Scaling_Law/blob/master/001_Exon_inclusion_levels_in_different_animals/PSI.sh
- http://www.currentprotocols.com/protocol/hg1116
- adding <path_to_bedtools2.23> for specify bedtools2.23

##### ExonicPartPSI_*.sh is organized from [doi: 10.1002/0471142905.hg1116s87](./doc/ExonicPart_PSI.pdf)
- add "awk '{if ($2 < 0) $2 = 0}{print $0}'" in Filtering junction
- add -sorted -g genomefile for large size chromosome

#### Dependency
- python2
- DEXseq, python2
- bedtools2.23

## Example, RNAseq sample from human
#### Make exonic part gff file
```bash
python2 path/to/dexseq_prepare_annotation.py path/to/Homo_sapiens.GRCh38.93.gtf Homo_sapiens.GRCh38.93_reduce.gtf

awk '{OFS="\t"}{if ($3 == "exonic_part") print $1,$2,$3,$4,$5,$6,$7,$8,$14":"$12}' Homo_sapiens.GRCh38.93_reduce.gtf | sed 's=[";]==g' > Homo_sapiens.GRCh38.93_Exonic_part.gff

rm Homo_sapiens.GRCh38.93_reduce.gtf
```

#### Make SJ bed file from all samples
```bash
awk 'BEGIN{OFS="\t"}{print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7, ($4 == 1)? "+":"-",$2-20-1, $3+20, "255,0,0", 2, "20,20", "0,300" }' example_1.SJ.out.tab > junctions.bed
```

#### Calculate PSI
##### Method 1, for junction.bed
```bash
## only count SJ of all annotated exon
bash path/to/PSI_1.bash StartPSI path/to/Homo_sapiens.GRCh38.93_Exonic_part.gff <reads_length> example_1.bam junctions.bed example_1 path/to/bedtools2.23

## count SJ of all in SJ.out.tab
bash path/to/PSI_2.bash StartPSI path/to/Homo_sapiens.GRCh38.93_Exonic_part.gff <reads_length> example_1.bam junctions.bed example_1 path/to/bedtools2.23
```

##### Method 2, for STAR output
```bash
## only count SJ of all annotated exon
 bash path/to/ExonicPartPSI_1.sh path/to/bedtools2.23/bedtools path/to/Homo_sapiens.GRCh38.93_Exonic_part.gff example_1.bam <reads_length> example_1.SJ.out.tab example_1

## count SJ of all in SJ.out.tab
 bash path/to/ExonicPartPSI_2.sh path/to/bedtools2.23/bedtools path/to/Homo_sapiens.GRCh38.93_Exonic_part.gff example_1.bam <reads_length> example_1.SJ.out.tab example_1

## only count SJ of all annotated exon (large chromosome)
 bash path/to/ExonicPartPSI_1g.sh path/to/bedtools2.23/bedtools path/to/Homo_sapiens.GRCh38.93_Exonic_part.gff genomefile example_1.bam <reads_length> example_1.SJ.out.tab example_1

## count SJ of all in SJ.out.tab (large chromosome)
 bash path/to/ExonicPartPSI_2g.sh path/to/bedtools2.23/bedtools path/to/Homo_sapiens.GRCh38.93_Exonic_part.gff genomefile example_1.bam <reads_length> example_1.SJ.out.tab example_1
```
