# CSE185-project

## Project Description
CSE 185 Final Project: Peak Calling Implementation for Mouse Liver ChIP-seq

## Title/Description: PEAK-A-VIEW
A peak calling tool to compete with findPeaks.

## Table of Contents
- [Contributions](#contributions)
- [Installation](#installation)
  - [Setting Up the Environment](#setting-up-the-environment)
- [Usage](#usage)
  - [Basic Usage Instructions](#basic-usage-instructions)
  - [Complete Usage Instructions](#complete-usage-instructions)
- [Project Structure](#project-structure)
- [Credits](#credits)

## Contributions
- Jiyeon Song: [jis036@ucsd.edu](mailto:jis036@ucsd.edu)
- Ivana Roque: [iroque@ucsd.edu](mailto:iroque@ucsd.edu)
- Jacob Ketchum: [dketchum@ucsd.edu](mailto:dketchum@ucsd.edu)

## Installation

### Prerequisites
- Python 3.8
- Required libraries: check the `requirements.txt`

### Setting Up the Environment
1. Clone the repository:
    ```sh
    git clone https://github.com/jiyeonsongg/peak-a-view.git
    cd peak-a-view
    ```

2. Create a virtual environment named `peak-a-view`:
    ```sh
    python -m venv peak-a-view
    ```

3. Activate the virtual environment:
    - On Windows:
      ```sh
      .\peak-a-view\Scripts\activate
      ```
    - On macOS/Linux:
      ```sh
      source peak-a-view/bin/activate
      ```

4. Install the required libraries:
    ```sh
    pip install -r requirements.txt
    ```

5. Install `peak_a_view.py`:
    ```sh
    python setup.py install
    ```

   Note: If you do not have root access, you can install locally:
    ```sh
    pip install --user -r requirements.txt
    python setup.py install --user
    ```

6. Verify the installation:
    ```sh
    peak_a_view --help
    ```

## Usage

### Basic Usage Instructions
A peak calling tool using a Python script and visualization methods. Compare with the peak calling algorithm findPeaks.

### Complete Usage Instructions
TBD

## Project Structure
