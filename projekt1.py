import numpy as np
from Bio.Align import substitution_matrices

MATCH       = +5
MISMATCH    = -4
GAP         = -2

STOP        = 0
DIAG        = 1
UP          = 2
LEFT        = 3

# function for extracting dna seq from fasta lines (only for 1 seq)
def get_seq(fasta_lines):
    dna = ""
    for line in fasta_lines:
        if line[0] != '>':
            dna += line.strip()
    return dna

# function to comapre
def score(a, b):
    if a == b:
        return MATCH
    else:
        return MISMATCH

# function implementing waterman for dna
def waterman_dna(seq1, seq2, score_matrix, traceback_matrix):

    seq1_len = len(seq1)
    seq2_len = len(seq2)

    best_score = 0
    best_y = best_x = 0

    for y in range(1, seq1_len + 1):
        for x in range(1, seq2_len + 1):
            d = score_matrix[y-1 , x-1] + score(seq1[y-1], seq2[x-1])
            u = score_matrix[y-1 , x] + GAP
            l = score_matrix[y , x-1] + GAP

            best = max(0, d, u, l)
            score_matrix[y, x] = best

            if best == 0:
                traceback_matrix[y, x] = STOP
            elif best == d:
                traceback_matrix[y, x] = DIAG
            elif best == u:
                traceback_matrix[y, x] = UP
            else:
                traceback_matrix[y, x] = LEFT

            if best > best_score:
                best_score = best
                best_y, best_x = y, x

    # print("score_matrix:\n")
    # print(score_matrix)

    # print("traceback_matrix:\n")
    # print(traceback_matrix)

    return best_y, best_x

# function implementing waterman for protein
def waterman_protein(seq1, seq2):
    best_score = 0
    best_y = best_x = 0

    for y in range(1, protein1_len + 1):
        for x in range(1, protein2_len + 1):
            d = protein_matrix_score[y-1 , x-1] + pam250[seq1[y-1] , seq2[x-1]]
            u = protein_matrix_score[y-1 , x] + GAP
            l = protein_matrix_score[y , x-1] + GAP

            best = max(0, d, u, l)
            protein_matrix_score[y, x] = best

            if best == 0:
                protein_matrix_traceback[y, x] = STOP
            elif best == d:
                protein_matrix_traceback[y, x] = DIAG
            elif best == u:
                protein_matrix_traceback[y, x] = UP
            else:
                protein_matrix_traceback[y, x] = LEFT

            if best > best_score:
                best_score = best
                best_y, best_x = y, x

    # print("protein_matrix_score:\n")
    # print(protein_matrix_score)

    # print("protein_matrix_traceback:\n")
    # print(protein_matrix_traceback)

    return best_y, best_x

def traceback(seq1, seq2, best_y, best_x, traceback_matrix):
    aligned1 = []
    aligned2 = []

    y, x = best_y, best_x

    while traceback_matrix[y, x] != STOP:
        direction = traceback_matrix[y, x]

        if direction == DIAG:
            aligned1.append(seq1[y-1])
            aligned2.append(seq2[x-1])
            y -= 1
            x -= 1
        elif direction == UP:
            aligned1.append(seq1[y-1])
            aligned2.append("-")
            y -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[x-1])
            x -= 1

    aligned1 = "".join(aligned1[::-1])
    aligned2 = "".join(aligned2[::-1])
    return aligned1, aligned2

def print_aligned(aligned1, aligned2):
    middle = ""
    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i]:
            middle += "|"
        else:
            middle += " "
    
    print(f"aligned1:\t{aligned1}")
    print(f"\t\t{middle}")
    print(f"aligned2:\t{aligned2}")

def main(path_file1, path_file2, dna):
    file1_fasta_lines = open(f"{path_file1}" , "r").readlines()
    file2_fasta_lines = open(f"{path_file2}" , "r").readlines()

    seq1 = get_seq(file1_fasta_lines)
    seq2 = get_seq(file2_fasta_lines)

    seq1_len = len(seq1)
    seq2_len = len(seq2)

    score_matrix = np.zeros((seq1_len + 1, seq2_len + 1))
    traceback_matrix = np.zeros((seq1_len + 1, seq2_len + 1))

    if dna:
        best_y, best_x = waterman_dna(seq1, seq2, score_matrix, traceback_matrix)
    else:
        best_y, best_x = waterman_protein(seq1, seq2)
    
    aligned1, aligned2 = traceback(seq1, seq2, best_y, best_x, traceback_matrix)

    print_aligned(aligned1, aligned2)

main("./dna1.fasta", "./dna2.fasta", 1)

# # dna seq
# dna1_fasta_lines = open("dna1.fasta" , "r").readlines()
# dna2_fasta_lines = open("dna2.fasta" , "r").readlines()

# dna1 = get_seq(dna1_fasta_lines)
# dna2 = get_seq(dna2_fasta_lines)

# dna1_len = len(dna1)
# dna2_len = len(dna2)

# score_matrix = np.zeros((dna1_len + 1, dna2_len + 1))
# traceback_matrix = np.zeros((dna1_len + 1, dna2_len + 1))

# best_score, best_y, best_x = waterman_dna(dna1, dna2)
# # print(best_score, best_y, best_x)

# aligned1, aligned2 = traceback(dna1, dna2, best_y, best_x)

# print_aligned(aligned1, aligned2)

# # protein seq
# protein1_fasta_lines = open("protein1.fasta" , "r").readlines()
# protein2_fasta_lines = open("protein2.fasta" , "r").readlines()

# pam250 = substitution_matrices.load("PAM250")
# # print("pam250:")
# # print(pam250)

# protein1 = get_seq(protein1_fasta_lines)
# protein2 = get_seq(protein2_fasta_lines)

# protein1_len = len(protein1)
# protein2_len = len(protein2)

# protein_matrix_score = np.zeros((protein1_len + 1, protein2_len + 1))
# protein_matrix_traceback = np.zeros((protein1_len + 1, protein2_len + 1))

# best_score, best_y, best_x = waterman_protein(protein1, protein2)
# # print(best_score, best_y, best_x)

# aligned1, aligned2 = traceback(protein1, protein2, best_y, best_x)
# print_aligned(aligned1, aligned2)
