#!/bin/bash
# creatation: 2023-08-04

# Stop on error
set -e



############################################
#ATACseq_pipleline
#Usage: 
#Example: 
############################################


#####	Parameters for ATAC-seq Data	###
bowtie2_index_path=$1
threads=$2
basename=$3
file1=$4
file2=$5
file1base=$6
file2base=$7
assembly=$8
spike_in=$9

# 1. Adapter Trimming
/share/home/zhongyiting/pipeline/ATAC-seq/code/adapterTrimmingModified $file1 $file2

# 2. Mapping
bowtie2 -p $threads -X 2000 --very-sensitive -x $bowtie2_index_path -1 ${file1base}.trim.fastq -2 ${file2base}.trim.fastq -S ${basename}.sam >> ${basename}_ATAC-seq_mapping_summary.txt 2>&1

# 3. ChrM 
awk '$3!="chrM"' ${basename}.sam |samtools view -S -b -f 0x2 -q 10 |samtools sort -o ${basename}.pe.q10.sort.bam
samtools view -Sb ${basename}.sam  > ${basename}.bam
awk '$3=="chrM" || NF<10 ' ${basename}.sam |samtools view -S -b > ${basename}.chrM.bam

# 4. Picard (duplicate removal) 
picard MarkDuplicates INPUT=${basename}.pe.q10.sort.bam \
OUTPUT=${basename}.pe.q10.rmdup.bam \
METRICS_FILE=${basename}.Picard_Metrics_unfiltered_bam.txt REMOVE_DUPLICATES=true >> ${basename}_ATAC-seq_mapping_summary.txt 2>&1

# 5. bedgraph, normalized bedgraph, normalized bw 
perl /share/home/zhongyiting/pipeline/ATAC-seq/code/bam2bed_shift.pl  ${basename}.pe.q10.rmdup.bam
bedtools slop -i ${basename}.pe.q10.rmdup.bed -g /share/Genomes/${assembly}/Sequence/${assembly}.chrom.sizes -b 0 > ${basename}.pe.q10.rmdup.cut.bed
genomeCoverageBed -bg -split -i ${basename}.pe.q10.rmdup.cut.bed -g /share/Genomes/${assembly}/Sequence/${assembly}.chrom.sizes > ${basename}.bg

##### normalize #####
bowtie2 -p $threads -X 2000 --very-sensitive -x /share/Genomes/${spike_in}/bowtie2/${spike_in} -1 ${file1base}.trim.fastq -2 ${file2base}.trim.fastq -S ${basename}.spike.sam >> ${basename}_ATAC-seq_mapping_summary.txt 2>&1
samtools view -Sb ${basename}.spike.sam  > ${basename}.spike.bam
picard MarkDuplicates INPUT=${basename}.spike.bam \
OUTPUT=${basename}.spike.rmdup.bam \
METRICS_FILE=${basename}.spike.Picard_Metrics_unfiltered_bam.txt REMOVE_DUPLICATES=true >> ${basename}_ATAC-seq_mapping_summary.txt 2>&1
perl /share/home/zhongyiting/pipeline/ATAC-seq/code/bam2bed_shift.pl  ${basename}.spike.rmdup.bam 
nor=$(wc -l ${basename}.spike.rmdup.bed | cut -d ' ' -f 1)
awk '{printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$4*100000000/'${nor}')}' ${basename}.bg > ${basename}.norm.bg
bedSort ${basename}.norm.bg ${basename}.norm.sort.bg
bedGraphToBigWig ${basename}.norm.sort.bg /share/Genomes/${assembly}/Sequence/${assembly}.chrom.sizes ${basename}.norm.bw

#perl /share/home/zhongyiting/pipeline/ATAC-seq/code/norm_bedGraph.pl ${basename}.bedGraph ${basename}.norm.bedGraph 

echo "count of total reads after QC filter" >> ${basename}_ATAC-seq_mapping_summary.txt
bedtools bamtobed -i ${basename}.pe.q10.sort.bam | wc -l >> ${basename}_ATAC-seq_mapping_summary.txt
echo "count of chrM reads" >> ${basename}_ATAC-seq_mapping_summary.txt
bedtools bamtobed -i ${basename}.chrM.bam | wc -l >> ${basename}_ATAC-seq_mapping_summary.txt
echo "final mapped reads" >> ${basename}_ATAC-seq_mapping_summary.txt
wc -l ${basename}.pe.q10.rmdup.bed >> ${basename}_ATAC-seq_mapping_summary.txt
echo "count of total reads before QC filter" >> ${basename}_ATAC-seq_mapping_summary.txt
bedtools bamtobed -i ${basename}.bam | wc -l >> ${basename}_ATAC-seq_mapping_summary.txt
echo "Total number of reads in the black list" >> ${basename}_ATAC-seq_mapping_summary.txt
bedtools intersect -a ${basename}.pe.q10.rmdup.bed -b /share/Genomes/${assembly}/Annotation/${assembly}.blacklist.bed -u | wc -l >> ${basename}_ATAC-seq_mapping_summary.txt



###### prepare the file for UCSC  ##########
lines=$(wc -l ${basename}.norm.sort.bg | cut -d ' ' -f 1)
if [ $lines -le 50000000 ]; then
sed -i '1i\track type=bedGraph name='${basename}' description='${basename}' color=0,0,0' ${basename}.norm.sort.bg
gzip ${basename}.norm.sort.bg

elif [ $lines -gt 50000000 ]; then
l=$(expr $line + 1)
a=$(expr $l / 2)
head ${basename}.norm.sort.bg -n $a > ${basename}.1.norm.sort.bg
tail ${basename}.norm.sort.bg -n $a > ${basename}.2.norm.sort.bg
sed -i '1i\track type=bedGraph name='${basename}.1' description='${basename}.1' color=0,0,0' ${basename}.1.norm.sort.bg
sed -i '1i\track type=bedGraph name='${basename}.2' description='${basename}.2' color=0,0,0' ${basename}.2.norm.sort.bg
gzip ${basename}.1.norm.sort.bg
gzip ${basename}.2.norm.sort.bg
fi



##### QC #####
samtools index ${basename}.pe.q10.rmdup.bam
#python /share/home/zhongyiting/pipeline/ATAC-seq/code/pyMakeVplot.py -a ${basename}.pe.q10.rmdup.bam -b /share/home/zhongyiting/data/reference/human/hg38/hg38.refGene.TSS.bed -p ends -e 2000 -u -v  -o ${basename}.TSSenrich

#perl /share/home/zhongyiting/pipeline/ATAC-seq/code/fragment_length_dist.pl ${basename}.pe.q10.rmdup.bam ${basename}.fragL.txt
#sort -n ${basename}.fragL.txt | uniq -c > ${basename}.frag.sort.txt
#Rscript /share/home/zhongyiting/pipeline/ATAC-seq/code/fragment_length_dist.R ${basename}.fragL.txt ${basename}.frag.sort.txt ${basename}.fragment_length_distribution.pdf ${basename}.fragment_length_distribution.txt
#rm ${basename}.fragL.txt
#rm ${basename}.frag.sort.txt

##### Peak Calling   #####
mkdir -p peak
macs2  callpeak -t ${basename}.pe.q10.rmdup.bed -f BED  -g hs -q 0.01 -n peak/${basename} --nomodel  --shift 0  
bedtools intersect -a peak/${basename}_peaks.narrowPeak -b /share/Genomes/${assembly}/Annotation/${assembly}.blacklist.bed -v > ${basename}_peaks.filterBL.bed

