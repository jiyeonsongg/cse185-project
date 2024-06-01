import pandas as pd
import subprocess

def bam_to_bed(bam_file, bed_file):
    subprocess.run(f"bedtools bamtobed -i {bam_file} > {bed_file}", shell=True)

def create_genome_file(bam_file, genome_file):
    subprocess.run(f"samtools view -H {bam_file} | grep '@SQ' | awk -F'\\t' '{{split($2,chr,\":\"); split($3,len,\":\"); print chr[2]\"\\t\"len[2]}}' > {genome_file}", shell=True)

def create_genome_windows(genome_file, windows_file, window_size=1000):
    subprocess.run(f"bedtools makewindows -g {genome_file} -w {window_size} > {windows_file}", shell=True)

def count_overlaps(reads, windows):
    overlap_counts = {index: 0 for index in range(len(windows))}
    for i, window in windows.iterrows():
        for _, read in reads.iterrows():
            if read['Chr'] == window[0] and not (read['end'] < window[1] or read['start'] > window[2]):
                overlap_counts[i] += 1
    return overlap_counts

def find_peaks(overlap_counts, threshold):
    peaks = [index for index, count in overlap_counts.items() if count >= threshold]
    return peaks

def main(input_bam):
    # Step 1: Convert BAM to BED
    bed_file = "output.bed"
    bam_to_bed(input_bam, bed_file)
    
    # Load BED file
    test_bed = pd.read_csv(bed_file, sep='\t', names=['Chr', 'start', 'end', 'read_name', 'quality_score', 'strand'])
    
    # Step 2: Combine chromosomes 1 to 12
    combined_df = pd.DataFrame(columns=['chrom', 'start', 'end'])
    for i in range(1, 13):
        chrom = f'chr{i}'
        chrom_df = test_bed[test_bed['Chr'] == chrom]
        if not chrom_df.empty:
            combined_row = {
                'chrom': chrom,
                'start': chrom_df['start'].min(),
                'end': chrom_df['end'].max()
            }
            combined_df = combined_df.append(combined_row, ignore_index=True)
    
    # Save the combined result to a new BED file
    combined_df.to_csv('combined_chroms.bed', sep='\t', header=False, index=False)
    
    # Step 3: Create genome file and windows
    genome_file = "genome_file.txt"
    create_genome_file(input_bam, genome_file)
    windows_file = "genome_windows.bed"
    create_genome_windows(genome_file, windows_file)
    
    # Load windows
    windows = pd.read_csv(windows_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
    
    # Step 4: Count overlaps in each window
    overlap_counts = count_overlaps(combined_df, windows)
    
    # Step 5: Identify peaks
    mean_count = sum(overlap_counts.values()) / len(overlap_counts)
    std_dev = (sum((x - mean_count) ** 2 for x in overlap_counts.values()) / len(overlap_counts)) ** 0.5
    threshold = mean_count + 2 * std_dev  # Adjust threshold if needed
    peaks = find_peaks(overlap_counts, threshold)
    
    # Step 6: Export peaks to BED file
    with open('peaks.bed', 'w') as file:
        for peak in peaks:
            p = windows.loc[peak]
            file.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\n")
    
    print("Peak calling completed and results saved to peaks.bed.")

# Example usage
if __name__ == "__main__":
    main("ENCFF609LFX.bam")
