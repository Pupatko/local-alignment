import numpy as np
from Bio.Align import substitution_matrices


MATCH       = +5
MISMATCH    = -4
GAP         = -2

STOP        = 0
DIAG        = 1
UP          = 2
LEFT        = 3

DNA         = 1
PROTEIN     = 0

pam250  = substitution_matrices.load("PAM250")

# function for extracting dna seq from fasta lines (only for 1 seq)
def get_seq(fasta_lines):
    dna = ""
    for line in fasta_lines:
        if line[0] != '>':
            dna += line.strip()
    return dna


# function to comapre the things in DNA
def score_dna(a, b):
    if a == b:
        return MATCH
    else:
        return MISMATCH


# function to comapre the things in PROTEIN
def score_protein(a, b):
    return pam250[a, b]


# function implementing waterman
def waterman(seq1, seq2, score_matrix, traceback_matrix, score):

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

    # print(score_matrix)
    # print(traceback_matrix)
    return best_y, best_x


# function to calculate the aligned seqs
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


# function for pretty print
def write_aligned(aligned1, aligned2, result_file_name):
    file = open(result_file_name, "w")
    middle = ""
    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i]:
            middle += "|"
        else:
            middle += " "
    
    file.write(f"aligned1:\t{aligned1}\n")
    file.write(f"\t\t\t{middle}\n")
    file.write(f"aligned2:\t{aligned2}\n")


# function for pretty print
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


# function to exec alignment
def main(path_file1, path_file2, dna, result_file_name):
    file1_fasta_lines = open(f"{path_file1}" , "r").readlines()
    file2_fasta_lines = open(f"{path_file2}" , "r").readlines()

    seq1 = get_seq(file1_fasta_lines)
    seq2 = get_seq(file2_fasta_lines)

    seq1_len = len(seq1)
    seq2_len = len(seq2)

    score_matrix = np.zeros((seq1_len + 1, seq2_len + 1))
    traceback_matrix = np.zeros((seq1_len + 1, seq2_len + 1))

    if dna:
        best_y, best_x = waterman(seq1, seq2, score_matrix, traceback_matrix, score_dna)
    else:
        best_y, best_x = waterman(seq1, seq2, score_matrix, traceback_matrix, score_protein)
    
    aligned1, aligned2 = traceback(seq1, seq2, best_y, best_x, traceback_matrix)

    print_aligned(aligned1, aligned2)
    write_aligned(aligned1, aligned2, result_file_name)


# filepath_1, filepath_2, DNA/PROTEIN, result_file_name
main("./dna1.fasta", "./dna2.fasta", DNA, "result_dna.txt")
main("./protein1.fasta", "./protein2.fasta", PROTEIN, "result_protein.txt")
