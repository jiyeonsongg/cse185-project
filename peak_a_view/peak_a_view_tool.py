import pandas as pd
import pybedtools
import subprocess

def bam_to_bed(bam_file, bed_file):
    subprocess.run(f"bedtools bamtobed -i {bam_file} > {bed_file}", shell=True)

def create_genome_file(bam_file, genome_file):
    subprocess.run(f"samtools view -H {bam_file} | grep '@SQ' | awk -F'\\t' '{{split($2,chr,\":\"); split($3,len,\":\"); print chr[2]\"\\t\"len[2]}}' > {genome_file}", shell=True)

def create_genome_windows(genome_file, windows_file, window_size=1000):
    subprocess.run(f"bedtools makewindows -g {genome_file} -w {window_size} > {windows_file}", shell=True)

def peak_a_view(input_bam, output_file=None, window_size=1000):
    # Step 1: Convert BAM to BED
    bed_file = "output.bed"
    bam_to_bed(input_bam, bed_file)
    
    # Load BED file
    test_bed = pd.read_csv(bed_file, sep='\t', names=['Chr', 'start', 'end', 'read_name', 'quality_score', 'strand'])
    # this BED file will contain the reads from the BAM file
    test_bed = test_bed.iloc[:,:3] # only need chromosome number, start, end
    test_bed.to_csv('reads.bed', sep = '\t' ,index=False, header=None) # we will count overlaps of reads with the windows
        
    # Step 2: Create genome file and windows BED file
    genome_file = "genome_file.txt"
    create_genome_file(input_bam, genome_file)
    windows_file = "genome_windows.bed"
    create_genome_windows(genome_file, windows_file, window_size)
    
    # Step 3: Count overlaps in each window
    # use pybedtools to count overlaps instead (more efficient)
    # load reads and windows from their respective BED files
    reads = pybedtools.BedTool('reads.bed')
    windows = pybedtools.BedTool('genome_windows.bed')
    
    # Find overlaps using intersect
    overlaps = windows.intersect(reads, c=True)
    
    # Convert to DataFrame
    overlap_df = overlaps.to_dataframe(names=['chr', 'start', 'end', 'count'])
    
    # Step 4: Identify peaks
    # Calculate peaks using empirical rule (top 1%)
    counts = overlap_df['count']
    threshold = counts.quantile(0.99)  # Top 1% of the overlap counts
    significant_peaks = overlap_df[overlap_df['count'] >= threshold]
    
    # Step 5: Export peaks to BED file
    output_file = output_file if output_file else 'peaks.bed'
    with open(output_file, 'w') as file:
        for i,peak in significant_peaks.iterrows():
            file.write(f"{peak['chr']}\t{peak['start']}\t{peak['end']}\n")
    
    print(f"Peak calling completed and results saved to {output_file}.")

# Example usage
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Peak calling tool for ChIP-seq data.")
    parser.add_argument("bam_file", type=str, help="Input BAM file.")
    parser.add_argument("-o", "--output", type=str, help="Output BED file.", default="peaks.bed")
    parser.add_argument("-w", "--window_size", type=int, help="Window size for peak calling.", default=1000)
    args = parser.parse_args()

    peak_a_view(args.bam_file, args.output, args.window_size)
