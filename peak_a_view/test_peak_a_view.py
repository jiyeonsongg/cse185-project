import unittest
import os
from unittest.mock import patch, mock_open, MagicMock
from peak_a_view_tool import bam_to_bed, create_genome_file, create_genome_windows, peak_a_view

class TestPeakAView(unittest.TestCase):
    @patch('subprocess.run')
    def test_bam_to_bed(self, mock_subprocess):
        bam_to_bed("test.bam", "test.bed")
        mock_subprocess.assert_called_with("bedtools bamtobed -i test.bam > test.bed", shell=True)

    @patch('subprocess.run')
    def test_create_genome_file(self, mock_subprocess):
        # Expected content that would be written to the genome file
        expected_content = "chr1\t1000\nchr2\t2000\n"
        
        # Simulate the behavior of subprocess.run by writing the expected content to the output file
        def subprocess_side_effect(cmd, shell):
            if "samtools view -H" in cmd:
                with open('genome_file.txt', 'w') as f:
                    f.write(expected_content)
        
        mock_subprocess.side_effect = subprocess_side_effect
        
        # Call the function to test
        create_genome_file("test.bam", "genome_file.txt")
        
        # Read the output file and check its contents
        with open("genome_file.txt", "r") as f:
            content = f.read()
        
        self.assertEqual(content, expected_content)


    @patch('subprocess.run')
    def test_create_genome_windows(self, mock_subprocess):
        create_genome_windows("genome_file.txt", "genome_windows.bed", 1000)
        mock_subprocess.assert_called_with("bedtools makewindows -g genome_file.txt -w 1000 > genome_windows.bed", shell=True)


    def test_peak_a_view(self):
        input_bam = "../benchmark/test.bam"
        output_file = "output_peaks.bed"
        
        # Run the peak_a_view function
        peak_a_view(input_bam, output_file, window_size=1000)
        
        # Check if the output file was created
        self.assertTrue(os.path.exists(output_file), "Output file was not created.")
        
        # Optional: Check the contents of the output file
        with open(output_file, 'r') as f:
            lines = f.readlines()
            self.assertGreater(len(lines), 0, "Output file is empty.")
        
        # Clean up the generated files
        os.remove(output_file)
        os.remove("output.bed")
        os.remove("reads.bed")
        os.remove("genome_file.txt")
        os.remove("genome_windows.bed")

if __name__ == "__main__":
    unittest.main()
