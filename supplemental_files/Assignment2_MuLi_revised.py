#!/usr/bin/env python
import sys
######################### Non-coding DNA in failed translation is counted ############
# No matter if DNA translation fails or not, when both START condon and STOP condon are found,
# the non-coding DNA is counted as the area before START codon and after STOP codon,
# when DNA fails, the non-coding DNA is counted as the area before START codon and whole DNA sequence if there is no START codon.
# 2. List of mini chromosomes
# To read file in command line:
# chmod +x trans_trans.py
# ./trans_trans.py mini_chromosomes.txt

def read_mini_chromosomes_from_file(file_path):
    with open(file_path, 'r') as file:
        mini_chromosomes = [line.strip() for line in file.readlines() if line.strip()]
    return mini_chromosomes

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: ./trans_trans.py <mini_chromosomes_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    mini_chromosomes = read_mini_chromosomes_from_file(file_path)

    print('\nMini chromosomes:')
    for chrom in mini_chromosomes:
        print(chrom)

# 3. Codon dictionary: https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#Standard_DNA_codon_table
# ATG is special case START, is not fully represented in dictionary
codon = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', 'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

print('\nCodon dictionary:')
# c is codon, aa is amino acid
for c, aa in codon.items():
    print(f'{c} to {aa}')

# 4. Function of translation of DNA to protein
# Initialize vars
num_encoded_proteins = 0
min_protein_len = float('inf')
max_protein_len = 0
total_protein_len = 0
num_dna_failures = 0
total_non_coding_dna = 0

# Define function: translation of codon to amino acid
def translate_codon(codon_seq):
    # If codon is in dictionary, return amino acid, else return empty string
    return codon.get(codon_seq, '')

# Chromosome iteration
for chro in chr_exp:
    # Find START codon
    start_ind = chro.find('ATG')
    # Find STOP codon
    stop_indices = [chro.find(codon, start_ind + 3) for codon in ['TAA', 'TAG', 'TGA']]
    # Find minimum index of STOP codon, if any, else -1
    stop_ind = min(ind for ind in stop_indices if ind != -1) if any(ind != -1 for ind in stop_indices) else -1

    # Calculate non-coding DNA length
    if start_ind != -1:  # START codon found
        non_coding_before_start = start_ind
    else:
        non_coding_before_start = len(chro)

    if stop_ind != -1:  # STOP codon found
        non_coding_after_stop = len(chro) - stop_ind - 3  # -3 to exclude STOP codon
    else:
        non_coding_after_stop = 0

    # Translation of DNA to protein
    protein = ''
    if start_ind != -1 and stop_ind != -1:  # Both START and STOP codons found
        for i in range(start_ind + 3, stop_ind, 3):  # From START to STOP codon
            codon_seq = chro[i:i+3]
            if codon_seq in codon:
                amino_acid = translate_codon(codon_seq)
                if amino_acid == 'STOP':  # STOP codon found
                    break
                protein += amino_acid
            else:
                num_dna_failures += 1
                break
    
    # Protein properties
    if protein: 
        num_encoded_proteins += 1
        protein_len = len(protein)
        if protein_len == 0:  # Encoded protein length is 0
            num_dna_failures += 1
        else:
            total_protein_len += protein_len
            min_protein_len = min(min_protein_len, protein_len)
            max_protein_len = max(max_protein_len, protein_len)
    else:  # Encodes a protein of length 0
        num_dna_failures += 1

    # Update total non-coding DNA length
    total_non_coding_dna += non_coding_before_start + non_coding_after_stop

# Calculate mean protein length
mean_protein_len = total_protein_len / num_encoded_proteins if num_encoded_proteins != 0 else 0

# Output properties
print('\nProperties of the DNA:')
print(f'Number of encoded proteins: {num_encoded_proteins}')
print(f'len of shortest protein: {min_protein_len} amino acids')
print(f'len of longest protein: {max_protein_len} amino acids')
print(f'Mean protein len: {mean_protein_len:.2f} amino acids')
print(f'Number of DNA failures: {num_dna_failures}')
print(f'Total amount of non-coding DNA: {total_non_coding_dna} nucleotides')

############################################################################################################
# Test cases:
test_case_doc =  (['ATGATTTAA'],
                  '''Test feature: encoding of a protein. 
                  Num encoded proteins: 1 
                  Shortest protein: 1 
                  Longest protein: 1 
                  Mean protein length: 1 
                  Number of DNA failures: 0 
                  Amount of non-coding DNA: 0''')

test_1 = (['ATGGTCTAA', 'CGTGTCCCATAA'],
          '''Test feature: encoding of a protein. 
          Num encoded proteins: 1
          Shortest protein: 1
          Longest protein: 1 
          Mean protein length: 1 
          Number of DNA failures: 1 
          Amount of non-coding DNA: 12 nucleotides''')

test_2 = (['CCCATGGTCTCATGACTTATGCGTGTCCCATAA'],
          '''Test feature: encoding of a protein. 
          Num encoded proteins: 1 
          Shortest protein: 2 
          Longest protein: 2 
          Mean protein length: 2 
          Number of DNA failures: 0 
          Amount of non-coding DNA: 21 nucleotides''')

test_3 = (['CTAGCATCGATCATGGATCGA'],
          '''Test feature: encoding of a protein. 
          Num encoded proteins: 0
          Shortest protein: 0 
          Longest protein: 0
          Mean protein length: 0 
          Number of DNA failures: 1 
          Amount of non-coding DNA: 12''')

test_4 = (['ATGGATCGTACGATCGTACGATCGTAAACG'],
          '''Test feature: encoding of a protein. 
          Num encoded proteins: 1 
          Shortest protein: 7 
          Longest protein: 7 
          Mean protein length: 7 
          Number of DNA failures: 0 
          Amount of non-coding DNA: 3 nucleotides''')

test_5 = (['CCCATGACGCAGATTATGTACGATTAG'],
          '''Test feature: encoding of a protein. 
          Num encoded proteins: 1 
          Shortest protein: 6 
          Longest protein: 6 
          Mean protein length: 6 
          Number of DNA failures: 0 
          Amount of non-coding DNA: 3''')

test_6 = (['ATGTAACAGATTATGTACGATTAGACC', 'CTTATGGCAGAT'],
          '''Test feature: encoding of a protein. 
          Num encoded proteins: 0 
          Shortest protein: 0 
          Longest protein: 0 
          Mean protein length: 0 
          Number of DNA failures: 2 
          Amount of non-coding DNA: 24 nucleotides''')

# def translation(dna, codon):
#     protein = ''
#     start_codon = False
    
#     for i in range(0, len(dna), 3):
#         # 3 codons at a time to translate
#         c = dna[i:i+3]
        
#         if c == 'ATG': # start
#             start_codon = True

#         elif c in ('TAA', 'TAG', 'TGA'): # stop
#             if start_codon:
#                 break

#         if start_codon and c in codon: # translate codon to amino acid
#             aa = codon[c]
#             # add to protein sequence
#             if aa != 'STOP': 
#                 protein += aa
    
#     return protein

# protein_num = 0
# protein_len = []
# failure_num = 0
# noncoding_num = 0

# for chro in chr_exp:
#     protein_seq = translation(chro, codon)
#     if protein_seq:
#         protein_num += 1
#         protein_len.append(len(protein_seq))
#     else:
#         failure_num += 1
#     noncoding_num += (len(chro) - len(protein_seq)) * 3

# # Properties calculation
# shortest_len = min(protein_len)
# longest_len = max(protein_len)
# mean_len = sum(protein_len) / len(protein_len) if protein_len else 0 # avoid dividing by 0

# print('\nDNA properties:')
# print(f'Number of encoded proteins: {protein_num}')
# print(f'len in amino acids of the shortest protein: {shortest_len}')
# print(f'len in amino acids of the longest protein: {longest_len}')
# print(f'Mean protein len in amino acids: {mean_len:.2f}')
# print(f'Number of DNA failures: {failure_num}')
# print(f'Total amount of non-coding DNA in nucleotides: {noncoding_num}')
