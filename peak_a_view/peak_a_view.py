import argparse
import subprocess
# samtools, bedtools, macs2
def main():
    parser = argparse.ArgumentParser(
        prog="peakfind",
        description="Command-line script to perform peak finding of BAM files"
    )
    
    parser.add_argument("bam", help="Indexed BAM file", type = str)
    
    parser.add_argument("-o","--out", help="Write output to file. Default: stdout", metavar="FILE", type=str, required=False)
    
    args = parser.parse_args()
    bam_file = args.bam
    
    subprocess.run(['./peak_a_view/peakfind.sh', bam_file])
    

    
if __name__=="__main__":
    main()

