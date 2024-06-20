# DNA_Analyzer

The **DNA Analyzer** is a Python tool designed for analyzing and comparing mini-chromosomes from DNA sequences provided in FASTA format. It includes functionalities for validating DNA sequences, identifying start and stop codons, encoding proteins, calculating statistics on encoded proteins, and detecting variants between a reference and subject mini-chromosome.

---

## Table of Contents

- [Features](#features)
- [Usage](#usage)
- [Installation](#installation)
- [Example](#example)
- [Testing](#testing)
- [Contributing](#contributing)

---

## Features

### MiniChromosome Class

- **Validation**: Validates the DNA sequence to ensure it contains only valid nucleotides (A, T, C, G).
- **Scanning**: Identifies the indices of the start codon ('ATG') and the nearest stop codon ('TAA', 'TAG', 'TGA').
- **Encoding**: Translates DNA sequences into amino acid sequences using a predefined codon dictionary.
- **Analysis**: Computes statistics such as the number of encoded proteins, lengths of proteins, number of DNA failures, and non-coding DNA length.

### VariantCaller Class

- **Validation**: Compares the protein structures of a subject mini-chromosome with that of a reference mini-chromosome.
- **Variant Calling**: Identifies variants (differences in amino acids) between the reference and subject mini-chromosomes.
- **Sample Variant Calling**: Processes multiple subject mini-chromosomes and reports variants for each.

### Main Functionality

- **Input**: Reads DNA sequences from a reference FASTA file and one or more subject FASTA files.
- **Processing**: Utilizes `MiniChromosome` to analyze sequences and `VariantCaller` to compare variants.
- **Output**: Displays analysis results including protein statistics and detected variants.

### Unit Testing

- Includes comprehensive unit tests using Python's `unittest` framework (`TestMiniChromosome` and `TestVariantCaller`).
- Tests cover functionalities such as sequence validation, scanning for codons, protein encoding, analysis calculations, and variant calling.

---

## Usage

To run the DNA Analyzer, use the following command format:

```bash
python dna_analyzer.py <reference_file> <subject_file1> [<subject_file2> ...]
```

- Replace `<reference_file>` with the path to the reference mini-chromosome DNA sequence in FASTA format.
- Replace `<subject_file1>`, `<subject_file2>`, etc., with paths to subject mini-chromosome DNA sequences in FASTA format.

---

## Installation

1. Ensure you have Python 3.x installed on your system.
2. Install necessary dependencies using pip:

```bash
pip install biopython
```

- **Biopython** is required for parsing FASTA files.

---

## Example

Suppose you have a reference mini-chromosome in `reference.fasta` and subject mini-chromosomes in `subject1.fasta` and `subject2.fasta`. Run the DNA Analyzer as follows:

```bash
python dna_analyzer.py reference.fasta subject1.fasta subject2.fasta
```

- The results will be printed to the console, including protein statistics and detected variants.

---

## Testing

The DNA Analyzer includes unit tests to ensure functionality and accuracy. To run the tests, use:

```bash
python test_analyzer.py
```

- This command discovers and runs all test cases in the current directory and subdirectories.

---

## Contributing

Contributions to the DNA Analyzer project are welcome. Here are some ways you can contribute:

- Report bugs and issues
- Suggest new features or enhancements
- Improve documentation
- Submit pull requests

---

## Author

- [Mu Li](mailto:mu.li@icahn.mssm.edu)

For questions or feedback, contact [My Email](mailto:mu.li@icahn.mssm.edu).

---
