#!/usr/bin/env python
# Mu Li, 2024-06-20

# Run the main program: python dna_analyzer.py <reference_file> <subject_file> in command line
# pip install biopython
import sys
from Bio import SeqIO


class MiniChromosome:
    """
    Represents a mini-chromosome from a DNA sequence in FASTA format.

    Attributes:
        sequence (str): DNA sequence of the mini-chromosome.
        proteins (list): List of proteins encoded by the mini-chromosome.
        non_coding_dna (int): Total amount of non-coding DNA nucleotides.
        num_encoded_proteins (int): Number of successfully encoded proteins.
        min_protein_len (int): Length of the shortest encoded protein.
        max_protein_len (int): Length of the longest encoded protein.
        total_protein_len (int): Total length of all encoded proteins.
        num_dna_failures (int): Number of failures encountered during analysis.
        start_index (int): Index of the start codon 'ATG' in the sequence.
        stop_index (int): Index of the stop codon in the sequence.
    """

    def __init__(self, sequence):
        """
        Initializes a MiniChromosome object.

        Args:
            sequence (str): DNA sequence of the mini-chromosome.
        """
        self.sequence = sequence
        self.proteins = []
        self.non_coding_dna = 0
        self.num_encoded_proteins = 0
        self.min_protein_len = float('inf')
        self.max_protein_len = 0
        self.total_protein_len = 0
        self.num_dna_failures = 0
        self.start_index = -1
        self.stop_index = -1

    def validate(self):
        """
        Validates the DNA sequence of the mini-chromosome.

        Returns:
            bool: True if the sequence contains only valid nucleotides (A, T, C, G), False otherwise.
        """
        valid_nucleotides = {'A', 'T', 'C', 'G'}
        return all(nuc in valid_nucleotides for nuc in self.sequence)

    def scan(self):
        """
        Scans the DNA sequence to find start and stop codons.

        Returns:
            tuple: Indices of the start and stop codons in the sequence.
        """
        self.start_index = self.sequence.find('ATG')
        stop_indices = [self.sequence.find(codon, self.start_index + 3) for codon in ['TAA', 'TAG', 'TGA']]
        self.stop_index = min((idx for idx in stop_indices if idx != -1), default=-1)
        return self.start_index, self.stop_index

    def encode(self):
        """
        Encodes proteins from the DNA sequence between start and stop codons.

        Returns:
            list: List of proteins translated from codons.
        """
        proteins = []
        if self.start_index != -1 and self.stop_index != -1:
            for i in range(self.start_index + 3, self.stop_index, 3):
                codon_seq = self.sequence[i:i + 3]
                amino_acid = translate_codon(codon_seq)
                if amino_acid == 'STOP':
                    break
                proteins.append(amino_acid)
            self.proteins = proteins  # Store the encoded proteins in self.proteins
        return proteins
    
    # def encode(self):
    #     """
    #     Encodes proteins from the DNA sequence between start and stop codons.

    #     Returns:
    #         str: Concatenated sequence of amino acids translated from codons.
    #     """
    #     if self.start_index != -1 and self.stop_index != -1:
    #         for i in range(self.start_index + 3, self.stop_index, 3):
    #             codon_seq = self.sequence[i:i + 3]
    #             amino_acid = translate_codon(codon_seq)
    #             if amino_acid == 'STOP':
    #                 break
    #             self.proteins.append(amino_acid)
    #         return ''.join(self.proteins)
    #     return ''

    def analyze(self):
        """
        Analyzes the mini-chromosome by validating, scanning for codons, and encoding proteins.

        Updates object attributes based on analysis results.
        """
        if not self.validate():
            self.num_dna_failures += 1
            self.non_coding_dna += len(self.sequence)
            return

        self.scan()
        protein = self.encode()

        if protein:
            self.num_encoded_proteins += 1
            protein_len = len(protein)
            if protein_len == 0:
                self.num_dna_failures += 1
            else:
                self.total_protein_len += protein_len
                self.min_protein_len = min(self.min_protein_len, protein_len)
                self.max_protein_len = max(self.max_protein_len, protein_len)
        else:
            self.num_dna_failures += 1

        if self.start_index != -1:
            self.non_coding_dna += self.start_index
        if self.stop_index != -1:
            self.non_coding_dna += len(self.sequence) - self.stop_index - 3

    def __str__(self):
        """
        Returns a string representation of the mini-chromosome analysis results.

        Returns:
            str: String containing analysis results.
        """
        mean_protein_len = (self.total_protein_len / self.num_encoded_proteins) if self.num_encoded_proteins != 0 else 0
        return (f"Number of encoded proteins: {self.num_encoded_proteins}\n"
                f"Length of shortest protein: {self.min_protein_len if self.min_protein_len != float('inf') else 0} amino acids\n"
                f"Length of longest protein: {self.max_protein_len} amino acids\n"
                f"Mean protein length: {mean_protein_len}\n"
                f"Number of DNA failures: {self.num_dna_failures}\n"
                f"Total amount of non-coding DNA: {self.non_coding_dna} nucleotides")


class VariantCaller:
    """
    A class to compare variants between a reference mini-chromosome and subject mini-chromosomes.

    Attributes:
        reference (MiniChromosome): The reference mini-chromosome for comparison.
    """

    def __init__(self, reference_chromosome):
        """
        Initializes a VariantCaller object with a reference mini-chromosome.

        Args:
            reference_chromosome (MiniChromosome): Reference mini-chromosome object.
        """
        self.reference = reference_chromosome

    def validate(self, subject_chromosome):
        """
        Validates if the subject chromosome has the same number and length of proteins as the reference.

        Args:
            subject_chromosome (MiniChromosome): Subject mini-chromosome object.

        Returns:
            bool: True if the subject chromosome has the same proteins structure, False otherwise.
        """
        if len(self.reference.proteins) != len(subject_chromosome.proteins):
            return False
        for ref_prot, subj_prot in zip(self.reference.proteins, subject_chromosome.proteins):
            if len(ref_prot) != len(subj_prot):
                return False
        return True

    def call(self, subject_chromosome):
        """
        Calls variants between the reference and subject chromosome proteins.

        Args:
            subject_chromosome (MiniChromosome): Subject mini-chromosome object.

        Returns:
            list: List of variant strings in the format 'p.[ref_aa][pos][subj_aa]'.
        """
        variants = []
        try:
            ref_sequence = self.reference.sequence[self.reference.start_index:self.reference.stop_index]
            subj_sequence = subject_chromosome.sequence[subject_chromosome.start_index:subject_chromosome.stop_index]

            for pos, (ref_aa, subj_aa) in enumerate(zip(ref_sequence, subj_sequence)):
                if ref_aa != subj_aa:
                    variant = f"p.{ref_aa}{pos + 1}{subj_aa}"
                    variants.append(variant)
        except Exception as e:
            print(f"Error occurred in call method: {e}")
        
        return variants

    def call_sample(self, subject_chromosomes):
        """
        Calls variants for multiple subject mini-chromosomes.

        Args:
            subject_chromosomes (list): List of subject MiniChromosome objects.

        Returns:
            list: List of variant results for each subject mini-chromosome.
        """
        results = []
        for subj_chrom in subject_chromosomes:
            validation_result = self.validate(subj_chrom)
            print(f"Validation result for subject chromosome: {validation_result}")
            print(f"Reference proteins: {self.reference.proteins}")
            print(f"Subject proteins: {subj_chrom.proteins}")

            if validation_result:
                variants = self.call(subj_chrom)
                print(f"Call returned variants: {variants}")

                if variants:
                    results.append(variants)
                    print(f"Variants for subject chromosome:")
                    for variant in variants:
                        print(variant)
                else:
                    results.append("No variants detected")
                    print("No variants detected for subject chromosome.")
            else:
                results.append("Error: Mismatch in protein lengths or number of proteins")
                print("Error: Mismatch in protein lengths or number of proteins")

        return results



    # def call(self, subject_chromosome):
    #     """
    #     Calls variants between the reference and subject chromosome proteins.

    #     Args:
    #         subject_chromosome (MiniChromosome): Subject mini-chromosome object.

    #     Returns:
    #         list: List of variant strings in the format 'p.[ref_aa][pos][subj_aa]'.
    #     """
    #     variants = []
    #     for i, (ref_prot, subj_prot) in enumerate(zip(self.reference.proteins, subject_chromosome.proteins)):
    #         for j in range(len(ref_prot)):
    #             if ref_prot[j] != subj_prot[j]:
    #                 variant = f"p.{ref_prot[j]}{j + 1}{subj_prot[j]}"
    #                 variants.append(variant)
    #     return variants

    # def call_sample(self, subject_chromosomes):
    #     """
    #     Calls variants for multiple subject mini-chromosomes.

    #     Args:
    #         subject_chromosomes (list): List of subject MiniChromosome objects.

    #     Returns:
    #         list: List of variant results for each subject mini-chromosome.
    #     """
    #     results = []
    #     for subj_chrom in subject_chromosomes:
    #         if self.validate(subj_chrom):
    #             variants = self.call(subj_chrom)
    #             results.append(variants)
    #         else:
    #             results.append("Error: Mismatch in protein lengths or number of proteins")
    #     return results


def translate_codon(codon_seq):
    """
    Translates a DNA codon sequence into an amino acid using a codon dictionary.

    Args:
        codon_seq (str): DNA codon sequence.

    Returns:
        str: Corresponding single-letter amino acid code or 'STOP' for stop codons.
    """
    codon_dict = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', 'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    return codon_dict.get(codon_seq, '')


def read_fasta_file(file_path):
    """
    Reads a FASTA file containing DNA sequences and returns a list of MiniChromosome objects.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        list: List of MiniChromosome objects.
    """
    mini_chromosomes = []
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            mini_chromosome = MiniChromosome(sequence)
            mini_chromosome.analyze()
            mini_chromosomes.append(mini_chromosome)
    return mini_chromosomes


def main(reference_file, subject_files):
    """
    Main function to analyze and compare mini-chromosomes.

    Args:
        reference_file (str): Path to the reference FASTA file.
        subject_files (list): List of paths to subject FASTA files.
    """
    # Reference genome
    reference_chromosomes = read_fasta_file(reference_file)
    if len(reference_chromosomes) != 1:
        print("Error: Reference file should contain exactly one mini-chromosome.")
        return

    reference = reference_chromosomes[0]
    print("Reference Mini-Chromosome Analysis:")
    print(reference)

    # Subject genomes
    subject_chromosomes = []
    for subject_file in subject_files:
        subject_chromosomes.extend(read_fasta_file(subject_file))

    # Call variants
    variant_caller = VariantCaller(reference)
    results = variant_caller.call_sample(subject_chromosomes)

    # Print results
    for i, result in enumerate(results):
        print(f"\nSubject Mini-Chromosome {i + 1} Variants:")
        if isinstance(result, list):
            for variant in result:
                print(variant)
        else:
            print(result)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python dna_analyzer.py <reference_file> <subject_file1> [<subject_file2> ...]")
    else:
        reference_file = sys.argv[1]
        subject_files = sys.argv[2:]
        main(reference_file, subject_files)
