#!/usr/bin/env python
# Mu Li, 2024-06-20
# Run: python test_analyzer.py in command line
# pip install first

import unittest
from dna_analyzer import MiniChromosome, VariantCaller, translate_codon

class TestMiniChromosome(unittest.TestCase):

    def setUp(self):
        # Create a MiniChromosome instance for testing
        self.reference = MiniChromosome("ATGTTTCGCTAA")
        self.subject = MiniChromosome("ATGTTTACTTAA")

    def test_validate(self):
        # Test validate method for reference and subject MiniChromosomes
        self.assertTrue(self.reference.validate())
        self.reference.sequence = "ATGn"
        self.assertFalse(self.reference.validate())

        self.assertTrue(self.subject.validate())
        self.subject.sequence = "ATGn"
        self.assertFalse(self.subject.validate())

    def test_scan(self):
        # Test scan method for reference and subject MiniChromosomes
        start_index_ref, stop_index_ref = self.reference.scan()
        self.assertEqual(start_index_ref, 0)
        self.assertEqual(stop_index_ref, 9)

        start_index_subj, stop_index_subj = self.subject.scan()
        self.assertEqual(start_index_subj, 0)
        self.assertEqual(stop_index_subj, 9)

    def test_encode(self):
        # Test encode method for reference and subject MiniChromosomes
        self.reference.scan()
        proteins_ref = self.reference.encode()
        self.assertEqual(proteins_ref, ['F', 'R'])

        self.subject.scan()
        proteins_subj = self.subject.encode()
        self.assertEqual(proteins_subj, ['F', 'T'])

    def test_analyze(self):
        # Test analyze method for reference and subject MiniChromosomes
        self.reference.analyze()
        self.assertEqual(self.reference.num_encoded_proteins, 1)
        self.assertEqual(self.reference.total_protein_len, 2)
        self.assertEqual(self.reference.min_protein_len, 2)
        self.assertEqual(self.reference.max_protein_len, 2)
        self.assertEqual(self.reference.num_dna_failures, 0)
        self.assertEqual(self.reference.non_coding_dna, 0)

        self.subject.analyze()
        self.assertEqual(self.subject.num_encoded_proteins, 1)
        self.assertEqual(self.subject.total_protein_len, 2)
        self.assertEqual(self.subject.min_protein_len, 2)
        self.assertEqual(self.subject.max_protein_len, 2)
        self.assertEqual(self.subject.num_dna_failures, 0)
        self.assertEqual(self.subject.non_coding_dna, 0)

    def test_str(self):
        # Test __str__ method for reference and subject MiniChromosomes
        self.reference.analyze()
        expected_output_ref = ("Number of encoded proteins: 1\n"
                               "Length of shortest protein: 2 amino acids\n"
                               "Length of longest protein: 2 amino acids\n"
                               "Mean protein length: 2.0\n"
                               "Number of DNA failures: 0\n"
                               "Total amount of non-coding DNA: 0 nucleotides")
        self.assertEqual(str(self.reference), expected_output_ref)

        self.subject.analyze()
        expected_output_subj = ("Number of encoded proteins: 1\n"
                                "Length of shortest protein: 2 amino acids\n"
                                "Length of longest protein: 2 amino acids\n"
                                "Mean protein length: 2.0\n"
                                "Number of DNA failures: 0\n"
                                "Total amount of non-coding DNA: 0 nucleotides")
        self.assertEqual(str(self.subject), expected_output_subj)

# class TestMiniChromosome(unittest.TestCase):
    
#     def setUp(self):
#         # Create a MiniChromosome instance for testing
#         self.mini_chromosome = MiniChromosome("TCCATGTCAAAATGC")

#     def test_validate(self):
#         # Test validate method
#         self.assertTrue(self.mini_chromosome.validate())
#         self.mini_chromosome.sequence = "ATGX"  # Invalid nucleotide 'X'
#         self.assertFalse(self.mini_chromosome.validate())

#     def test_scan(self):
#         # Test scan method
#         start_index, stop_index = self.mini_chromosome.scan()
#         self.assertEqual(start_index, 3)
#         self.assertEqual(stop_index, -1) # No stop codon found

#     def test_encode(self):
#         # Test encode method
#         self.mini_chromosome.scan()
#         proteins = self.mini_chromosome.encode()
#         self.assertEqual(proteins, [])

#     def test_analyze(self):
#         # Test analyze method
#         self.mini_chromosome.analyze()
#         self.assertEqual(self.mini_chromosome.num_encoded_proteins, 0)
#         self.assertEqual(self.mini_chromosome.total_protein_len, 0)
#         self.assertEqual(self.mini_chromosome.min_protein_len, float('inf'))
#         self.assertEqual(self.mini_chromosome.max_protein_len, 0)
#         self.assertEqual(self.mini_chromosome.num_dna_failures, 1)
#         self.assertEqual(self.mini_chromosome.non_coding_dna, 3)

#     def test_str(self):
#         # Test __str__ method
#         self.mini_chromosome.analyze()
#         expected_output = ("Number of encoded proteins: 0\n"
#                            "Length of shortest protein: 0 amino acids\n"
#                            "Length of longest protein: 0 amino acids\n"
#                            "Mean protein length: 0\n"
#                            "Number of DNA failures: 1\n"
#                            "Total amount of non-coding DNA: 3 nucleotides")
#         self.assertEqual(str(self.mini_chromosome), expected_output)

# # Wrong name by using codon and position
# class TestVariantCaller(unittest.TestCase):

#     def setUp(self):
#         # Create reference MiniChromosome instance for testing
#         self.reference = MiniChromosome("ATGTTTACCTAA")
#         self.reference.analyze()  # Important: Set up correctly!!!

#         # Create subject MiniChromosome instance for testing
#         self.subject = MiniChromosome("ATGTTTACTTAA")
#         self.subject.analyze() 

#     def test_validate_same_structure(self):
#         # Test validate method with subject having the same structure as reference
#         variant_caller = VariantCaller(self.reference)
#         self.assertTrue(variant_caller.validate(self.subject))

#     def test_validate_different_structure(self):
#         # Test validate method with subject having different protein structure than reference
#         # Create a subject with different structure
#         self.subject_different = MiniChromosome("ATGCAAACCTTTTAA")
#         self.subject_different.proteins = [['M'], ['L']]
#         self.subject_different.analyze()
        
#         variant_caller = VariantCaller(self.reference)
#         self.assertFalse(variant_caller.validate(self.subject_different))

#     def test_call_no_variants(self):
#         # Test call method when there are no variants between reference and subject
#         self.subject_nov = MiniChromosome("ATGTTTACCTAA")
#         self.subject_nov.analyze()
#         variant_caller = VariantCaller(self.reference)
#         variants = variant_caller.call(self.subject_nov)
#         self.assertEqual(variants, [])

#     def test_call_with_variants(self):
#         # Test call method when there are variants between reference and subject
#         variant_caller = VariantCaller(self.reference)
#         variants = variant_caller.call(self.subject)
#         expected_variants = ['p.C9T']
#         self.assertEqual(variants, expected_variants)

#     def test_call_sample(self):
#         # Test call_sample method with a list of subject MiniChromosome instances
#         subject_chromosomes = [self.subject]
#         variant_caller = VariantCaller(self.reference)
#         results = variant_caller.call_sample(subject_chromosomes)
#         expected_results = [['p.C9T']]
#         self.assertEqual(results, expected_results)

# Correct name by using amino acids and position
class TestVariantCaller(unittest.TestCase):

    def setUp(self):
        # Create reference MiniChromosome instance for testing
        self.reference = MiniChromosome("ATGTTTCCCTAA")
        self.reference.analyze()  # Important: Set up correctly!!!

        # Create subject MiniChromosome instance for testing
        self.subject = MiniChromosome("ATGTTTACTTAA")
        self.subject.analyze() 

    def test_validate_same_structure(self):
        # Test validate method with subject having the same structure as reference
        variant_caller = VariantCaller(self.reference)
        self.assertTrue(variant_caller.validate(self.subject))

    def test_validate_different_structure(self):
        # Test validate method with subject having different protein structure than reference
        # Create a subject with different structure
        self.subject_different = MiniChromosome("ATGCAAACCTTTTAA")
        self.subject_different.proteins = [['M'], ['L']]
        self.subject_different.analyze()
        
        variant_caller = VariantCaller(self.reference)
        self.assertFalse(variant_caller.validate(self.subject_different))

    def test_call_no_variants(self):
        # Test call method when there are no variants between reference and subject
        self.subject_nov = MiniChromosome("ATGTTTCCTTAA")
        self.subject_nov.analyze()
        variant_caller = VariantCaller(self.reference)
        variants = variant_caller.call(self.subject_nov)
        self.assertEqual(variants, [])

    def test_call_with_variants(self):
        # Test call method when there are variants between reference and subject
        variant_caller = VariantCaller(self.reference)
        variants = variant_caller.call(self.subject)
        expected_variants = ['p.P2T']
        self.assertEqual(variants, expected_variants)

    def test_call_sample(self):
        # Test call_sample method with a list of subject MiniChromosome instances
        subject_chromosomes = [self.subject]
        variant_caller = VariantCaller(self.reference)
        results = variant_caller.call_sample(subject_chromosomes)
        expected_results = [['p.P2T']]
        self.assertEqual(results, expected_results)

if __name__ == "__main__":
    unittest.main()
