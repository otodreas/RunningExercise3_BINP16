#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 17:22:44 2025

@author: chetanav
"""
"""
Author: Chetana Varada
Date: 24th October, 2025
PART 1-2 (Student 1) of RE - 3 with Yaqi Jiao and Oliver Ekström Todreas

ExamineMSA.py
--------------------
Examines the output files from FastaParser.py - output_fasta_file_mtDNA.fasta output_fasta_file_Y.fasta 
and performs MSA on each individually and produces 2 text files:
  1. One for mtDNA sequences
  2. One for Y chromosome sequences

Description:
    The result fasta files from previous code is used as inputs here, with a weight
    file. The fasta files are read and the sequences are read in aligned form. Pairwise 
    score is generated for all possibilities of nucleotides and multiple sequence
    alignment is performed on each fasta file sequences separately. The output files 
    are written the the given format.
    
Procedure:
    Functions created:
        To read Fasta file into dictionary
        To read a weight matrix
        To calculate identity and alignment score for a pair of sequences
        To write results in a text file
        The main function to perform quality checks and print output statement

Usage:
    python3 ExamineMSA.py output_fasta_file_mtDNA.fasta output_fasta_file_Y.fasta weight_parameters output_MSA_file

Outputs Example:
    output_MSA_file_mtDNA_scores.txt
    output_MSA_file_Y_scores.txt
    
"""

import sys
import itertools

def read_fasta(fasta_file): #To read Fasta file into dictionary
    seq_dict = {} #Initialize empty dictionary to store sequences
    with open(fasta_file, 'r') as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if not line: #To skip empty lines
                continue
            if line.startswith('>'):
                if header:
                    seq_dict[header] = ''.join(seq).upper() #Saves previous sequence
                header = line[1:].strip() #Stores new header
                seq = [] #resets sequence
            else:
                seq.append(line) #To handle multi-line sequences
        if header:
            seq_dict[header] = ''.join(seq).upper() #For the last sequence
    return seq_dict #returns a dictionary mapping sequence names to their aligned sequences

def read_weight_matrix(weight_file): #Reads a weight matrix
    matrix = {} #initialize weights in this dictionary
    with open(weight_file, 'r') as f:
        for line in f:
            if not line.strip(): #Skip empty lines
                continue
            parts = line.strip().split() #nucleotide1, nucleotide2, score
            if len(parts) != 3:
                continue
            a, b, score = parts
            matrix[(a.upper(), b.upper())] = int(score)
    return matrix

def pairwise_score(seq1, seq2, weight_matrix): #To calculate identity and alignment score for a pair of sequences
    assert len(seq1) == len(seq2) #Ensures sequences are same length (aligned)
    total_positions = len(seq1)
    identical = 0
    comparable = 0
    score_sum = 0

    for a, b in zip(seq1, seq2):
        if a in ('?', '-') or b in ('?', '-'): #Treats '?' and '-' as gaps
            continue
        comparable += 1 #Counts positions where both are nucleotides
        if a == b:
            identical += 1 #Counts positions where nucleotides match
        score = weight_matrix.get((a, b)) #Gets score from weight matrix
        if score is None:
            score = weight_matrix.get((b, a), 0)  #Use symmetry if score not found
        score_sum += score #Sum of alignment scores of all comparable positions

    uncertain = total_positions - comparable #Counts positions with ? or - or both
    if total_positions > 0:
        identity_percent = 100 * identical / total_positions #% of matching nucleotides over total length
    else:
        0
    return comparable, uncertain, total_positions, identity_percent, score_sum

def write_results(seq_dict, weight_matrix, output_file): #To write results in a text file
    headers = list(seq_dict.keys()) #Gets sequence names
    results = []

    for h1, h2 in itertools.combinations(headers, 2): #Generates all unique pairs of sequence
        seq1 = seq_dict[h1]
        seq2 = seq_dict[h2]
        comp, uncertain, total, identity, score = pairwise_score(seq1, seq2, weight_matrix)
        results.append((h1, h2, comp, uncertain, total, identity, score))

    with open(output_file, 'w') as out:
        out.write("SampleA | SampleB | Comparable nucleotides | uncertain nucleotides | Total nucleotides | IdentityScore | Score\n")
        for r in results:
            out.write(f"{r[0]} {r[1]} {r[2]} {r[3]} {r[4]} {r[5]:.1f}% {r[6]}\n")

    print(f"Results written to {output_file}")

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 ExamineMSA.py output_fasta_file_mtDNA.fasta output_fasta_file_Y.fasta weight_parameters output_MSA_file")
        sys.exit(1)

    mtDNA_file = sys.argv[1]
    Y_file = sys.argv[2]
    weight_file = sys.argv[3]
    output_prefix = sys.argv[4]

    weight_matrix = read_weight_matrix(weight_file)

    # Process mtDNA
    mtDNA_dict = read_fasta(mtDNA_file)
    mtDNA_output = f"{output_prefix}_mtDNA_scores.txt"
    write_results(mtDNA_dict, weight_matrix, mtDNA_output)

    # Process Y chromosome
    Y_dict = read_fasta(Y_file)
    Y_output = f"{output_prefix}_Y_scores.txt"
    write_results(Y_dict, weight_matrix, Y_output)

if __name__ == "__main__":
    main()
