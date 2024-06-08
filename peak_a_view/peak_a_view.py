import argparse
import subprocess
import os

# ANSI escape sequences for text formatting
RESET = "\033[0m"
BOLD = "\033[1m"
BLUE = "\033[34m"
PINK = "\033[35m"
CYAN = "\033[36m"

def main():
    parser = argparse.ArgumentParser(
        prog=f"{BLUE}peak_a_view{RESET}",
        description=(
            f"{PINK}============================================================={RESET}\n"
            f"{BOLD}                Peak Calling Command Line Tool                {RESET}\n"
            f"{PINK}============================================================={RESET}\n"
            f"This tool performs peak calling on BAM files.\n"
            "It processes the BAM file to identify regions of interest (peaks) \n"
            "using custom algorithms and generates output suitable for \n"
            "visualization in IGV.\n"
            f"{PINK}-------------------------------------------------------------{RESET}\n"
            "Example usage:\n"
            f"{BLUE}peak_a_view{RESET} {CYAN}input.bam{RESET} {PINK}-o output_peaks.bed -w 500{RESET}\n"
            f"{PINK}============================================================={RESET}"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("bam", help="Indexed BAM file", type=str)
    parser.add_argument("-o", "--out", help="Write output to file. Default: output_peaks", metavar="FILE", type=str, required=False)
    parser.add_argument("-w", "--window_size", help="Window size for peak calling. Default: 1000", metavar="WINDOW_SIZE", type=int, default=1000)
    
    args = parser.parse_args()
    bam_file = args.bam
    output_file = args.out if args.out else "peaks.bed"
    window_size = args.window_size

    # Get the absolute path of the peak_a_view_tool.py script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    tool_path = os.path.join(script_dir, 'peak_a_view_tool.py')

    # Print paths for debugging
    print(f"Script directory: {script_dir}")
    print(f"Tool path: {tool_path}")

    # Call the peak_a_view function from peak_a_view_tool.py
    subprocess.run(['python', tool_path, bam_file, '-o', output_file, '-w', str(window_size)], check=True)

if __name__ == "__main__":
    main()
