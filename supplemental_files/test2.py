import unittest
from io import StringIO
from unittest.mock import patch
from dna_analyzer import MiniChromosome, VariantCaller, translate_codon

class TestMiniChromosome(unittest.TestCase):

    def setUp(self):
        self.valid_sequence = "ATGTTGTGCTAA"
        self.invalid_sequence = "ATGTTGNGCTAA"

    def test_validate_valid_sequence(self):
        mini_chromosome = MiniChromosome(self.valid_sequence)
        self.assertTrue(mini_chromosome.validate())

    def test_validate_invalid_sequence(self):
        mini_chromosome = MiniChromosome(self.invalid_sequence)
        self.assertFalse(mini_chromosome.validate())

    def test_scan(self):
        mini_chromosome = MiniChromosome(self.valid_sequence)
        start, stop = mini_chromosome.scan()
        self.assertEqual(start, 0)
        self.assertEqual(stop, 9)

    def test_encode(self):
        mini_chromosome = MiniChromosome(self.valid_sequence)
        mini_chromosome.scan()
        proteins = mini_chromosome.encode()
        self.assertEqual(proteins, ['L', 'C'])

    def test_analyze_valid_sequence(self):
        mini_chromosome = MiniChromosome(self.valid_sequence)
        mini_chromosome.analyze()
        self.assertEqual(mini_chromosome.num_encoded_proteins, 1)
        self.assertEqual(mini_chromosome.min_protein_len, 2)
        self.assertEqual(mini_chromosome.max_protein_len, 2)
        self.assertEqual(mini_chromosome.total_protein_len, 2)

    def test_analyze_invalid_sequence(self):
        mini_chromosome = MiniChromosome(self.invalid_sequence)
        mini_chromosome.analyze()
        self.assertEqual(mini_chromosome.num_dna_failures, 1)
        self.assertEqual(mini_chromosome.non_coding_dna, len(self.invalid_sequence))

class TestVariantCaller(unittest.TestCase):
    def setUp(self):
        reference_sequence = "ATGCGTTAACGT" 
        subject_sequence = "ATGCGTTAAGGT"
        
        self.reference = MiniChromosome(reference_sequence)
        print(f"Reference initial: {self.reference}")
        self.reference.analyze()
        print(f"Reference proteins: {self.reference.proteins}")

        self.subject = MiniChromosome(subject_sequence)
        self.subject.analyze()
        print(f"Subject proteins: {self.subject.proteins}")

    def test_call(self):
        variant_caller = VariantCaller(self.reference)
        print(f"Reference11: {self.reference}")
        variants = variant_caller.call(self.subject)
        print(f"Subject11: {self.subject}")
        print(f"Variants: {variants}")
        self.assertEqual(variants, ['p.A10G'])

if __name__ == "__main__":
    unittest.main()

#     def setUp(self):
#         reference_sequence = "ATGTTGTGCTAATAG"
#         subject_sequence = "ATGTTGCGCTAATAG"
#         self.reference = MiniChromosome(reference_sequence)
#         self.subject = MiniChromosome(subject_sequence)

#     def test_validate_same_structure(self):
#         variant_caller = VariantCaller(self.reference)
#         result = variant_caller.validate(self.subject)
#         self.assertTrue(result)

#     def test_validate_different_structure(self):
#         self.subject.sequence = "ATGTTGCGCTAATAAC"
#         variant_caller = VariantCaller(self.reference)
#         result = variant_caller.validate(self.subject.sequence)
#         self.assertFalse(result)

#     def test_call(self):
#         variant_caller = VariantCaller(self.reference)
#         variants = variant_caller.call(self.subject)
#         self.assertEqual(variants, ['p.T7G'])

#     @patch('sys.stdout', new_callable=StringIO)
#     def test_call_sample(self, mock_stdout):
#         subject_chromosomes = [self.subject]
#         variant_caller = VariantCaller(self.reference)
#         results = variant_caller.call_sample(subject_chromosomes)
#         expected_output = "Validation result for subject chromosome: True\n" \
#                           "Reference proteins: ['M', 'C']\n" \
#                           "Subject proteins: ['M', 'A']\n" \
#                           "Call returned variants: ['p.T4G']\n" \
#                           "Variants for subject chromosome:\n" \
#                           "p.T4G\n"
#         self.assertEqual(mock_stdout.getvalue(), expected_output)

# if __name__ == '__main__':
#     unittest.main()

