#!/bin/bash

#chmod +x your_script.sh
bam_file="$1"
#head -n 10 ../ENCFF609LFX.bam
# 2. sort and index the BAM file
samtools sort -o sorted_input.bam ${bam_file}
samtools index sorted_input.bam 
# 3. make genome_file.txt
samtools view -H ${bam_file} | grep "@SQ" | awk -F'\t' '{split($2,chr,":"); split($3,len,":"); print chr[2]"\t"len[2]}' > genome_file.txt
# 4. creating genomic windows with 'bedtools'
bedtools makewindows -g genome_file.txt -w 1000 > genome_windows.bed
# 5. counting reads in each window
bedtools coverage -a genome_windows.bed -b sorted_input.bam > coverage_in_windows.bed
# 6. peak calling using macs2
macs2 callpeak -t sorted_input.bam -f BAM -g hs -n output --outdir peak_calling_output

#!/usr/bin/env python
# coding: utf-8
# 1. create a small BAM file
# Get the total number of reads
#total_reads=$(samtools view -c ENCFF609LFX.bam)
# Calculate 20% of the total reads
#small_reads=$((total_reads / 5))
# Extract the first 20% of reads
#samtools view -h ENCFF609LFX.bam | head -n $((small_reads + 1000)) | samtools view -b -o small_ENCFF609LFX.bam

# 7. visualizing the peaks
# echo "[genes]
# file = /home/jis036/final-project/Mus_musculus.GRCm38.104.gtf
# height = 4
# color = black
# labels = true
# fontsize = 10

# [bam file]
# file = /home/jis036/final-project/sorted_small_ENCFF609LFX.bam
# height = 2
# color = blue
# min_value = 0
# max_value = 50

# [bed file]
# file = /home/jis036/final-project/peak_calling_output/small_ENCFF609LFX_peaks.narrowPeak
# height = 2
# color = red" > tracks.ini


# In[ ]:


# 8. generate the plot using IGV
# Select sorted_small_ENCFF609LFX.bam and sorted_small_ENCFF609LFX.bam.bai
# Select peak_calling_output/small_ENCFF609LFX_peaks.narrowPeak.

