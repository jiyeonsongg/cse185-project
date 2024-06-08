# PEAK-A-VIEW ( •̀ ω •́ )✧

## Project Description
CSE 185 Final Project: Peak Calling Implementation for Mouse Liver ChIP-seq

## Title/Description: PEAK-A-VIEW
```peak-a-view``` is designed to process ChIP-Seq data to identify peaks corresponding to key regions representing protein-DNA interactions.

## Table of Contents
- [Contributions](#contributions)
- [Installation](#installation)
  - [Setting Up the Environment](#setting-up-the-environment)
- [Usage](#usage)
  - [Basic Usage Instructions](#basic-usage-instructions)
  - [Complete Usage Instructions](#complete-usage-instructions)
- [Project Structure](#project-structure)

## Contributions
- Jiyeon Song: [jis036@ucsd.edu](mailto:jis036@ucsd.edu)
- Ivana Roque: [iroque@ucsd.edu](mailto:iroque@ucsd.edu)
- Jacob Ketchum: [dketchum@ucsd.edu](mailto:dketchum@ucsd.edu)

## Installation

### Setting Up the Environment
1. Clone the repository:
    ```sh
    git clone https://github.com/jiyeonsongg/peak-a-view.git
    cd peak-a-view
    ```
2. Set up the requirements:
    ```sh
    pip install -r requirements.txt
    ```

3. Install `peak_a_view.py`:
    ```sh
    python setup.py install
    ```

4. Verify the installation:
    ```sh
    which peak_a_view  # or `where peak_a_view` on Windows
    ```
5. Get started!
   ```sh
   peak_a_view --help
   ```
   EXAMPLE:
   ```sh
   peak_a_view small_ENCFF609LFX.bam -o output_peaks.bed -w 500
   # peak_a_view INPUT.bam -o OUTPUT.bed -w WINDOW_SIZE_NUMBER
   ```
   
## Usage

### Basic Usage Instructions
A peak calling tool using a Python script and visualization methods. Compare with the peak calling algorithm findPeaks.

### Complete Usage Instructions
Once we export the output_peaks.bed after running the command line, open IGV(Integrative Genomics Viewer) and upload the file to see peak calling visualization. &nbsp; 
Recommendation: Download IGV Desktop version - https://igv.org/

## Project Structure
```sh
peak-a-view
|
|--analysis
|--benchmark
|--peak_a_view
|  |  __init__.py
|  |  peak_a_view.py
|  |  peak_a_view_tool.py
|  |  test_peak_a_view.py
|--.gitignore
|--README.md
|--requirements.txt
|--setup.py
```

## Testing
To run unit tests: 
- Ensure you are in peak-a-view/peak_a_view directory
- Run: python test_peak_a_view.py
