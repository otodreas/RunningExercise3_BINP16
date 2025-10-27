#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 22:39:02 2025

@author: chetanav
"""

"""
Author: Chetana Varada
Date: 23rd October, 2025
PART 1-1 (Student 1) of RE - 3 with Yaqi Jiao and Oliver Ekström Todreas

FastaParser.py
--------------------
Parses a structured GeneticData.txt file containing DNA sequences
and produces two FASTA files:
  1. One for mtDNA sequences
  2. One for Y chromosome sequences

Description:
    A fasta file -> GeneticData - 4 was selected according to the out group number,
    which had all the ancestoral genetic information of Anastasia and the probable 
    genetic information of real Anastasia.
    This file is then decoded to produce 2 output files - mtDNA and Y.
    
Procedure:
    Functions created:
        To check if the line is a probable sequence, which contains A C T G ? -
        To check for line which have certain keywords
        Parse GeneticData.txt and return two dicts:
          mtDNA_dict and Y_dict
        Handles:
          - multi-line names (e.g., 'Anastasia4 son')
          - weird apostrophe bytes (<92>)
          - both mtDNA and Y chromosome sections
        To write a FASTA file
        The main function to perform quality checks and print output statement

Usage:
    python3 FastaParser.py GeneticData.txt output_fasta_file

Outputs Example:
    output_fasta_file_mtDNA.fasta
    output_fasta_file_Y.fasta
    
"""

import sys
#import os

def is_probable_sequence(line): #To Check if a line looks like a DNA sequence (A/T/C/G/?/-)
    return all(ch in "ATCG?- " for ch in line.strip()) and len(line.strip()) > 0

def is_metadata_line(line): #To check lines that indicate traits by keywords
    keywords = ["carrier", "patient", "not a", "hemophilia", "unknown"]
    low = line.lower() #convert to lower case
    return any(k in low for k in keywords) #match keyword to words in line in input file

def parse_genetic_data(input_file):
    mtDNA_dict = {}
    Y_dict = {}

    with open(input_file, 'r', encoding='utf-8', errors='ignore') as fh: #To read and clean file
        raw_lines = [ln.rstrip() for ln in fh.readlines()]

    lines = [] #clean the file to remove and replace certain characters and make it readable
    for ln in raw_lines:
        fixed = ln.replace("<92>", "'").replace("’", "'").strip()
        if fixed: #Removes empty lines
            lines.append(fixed)

    i = 0 #Current line index
    n = len(lines) #Total number of lines

    while i < n:
        line = lines[i].strip()

        if not line or is_metadata_line(line) or is_probable_sequence(line) or line.lower().startswith(("mtdna", "y chromosome")): 
                #To skip blank or metadata lines
            i += 1
            continue

        person_name = line.lstrip('>').strip() #To detect person name and handle multi-line names

        #Merge consecutive name lines into one if they’re not DNA, metadata, or labels
        while (i + 1 < n and
               not is_probable_sequence(lines[i + 1]) and
               not is_metadata_line(lines[i + 1]) and
               not lines[i + 1].lower().startswith(("mtdna", "y chromosome"))):
            person_name += " " + lines[i + 1].strip()
            i += 1

        i += 1  #move past name line

        while i < n:
            current = lines[i].strip().lower()

            #mtDNA
            if current.startswith("mtdna"):
                if i + 1 < n and is_probable_sequence(lines[i + 1]): #reads next line - sequence
                    seq = lines[i + 1].replace(" ", "").upper()
                    mtDNA_dict[f"{person_name}_mtDNA"] = seq #saves the sequence with person name mtDNA
                    i += 2
                    continue
                else:
                    i += 1
                    continue

            #Y chromosome
            elif current.startswith("y chromosome") or current.startswith("ychr") or current == "y":
                if i + 1 < n and is_probable_sequence(lines[i + 1]): #reads next line - sequence
                    seq = lines[i + 1].replace(" ", "").upper()
                    Y_dict[f"{person_name}_Y"] = seq #saves the sequence with person name Y chromosome
                    i += 2
                    continue
                else:
                    i += 1
                    continue

            if (not is_probable_sequence(lines[i]) 
                and not lines[i].lower().startswith(("mtdna", "y chromosome"))
                and not is_metadata_line(lines[i])): #Break, if next line is a new name (not sequence or label)
                break

            i += 1

    return mtDNA_dict, Y_dict


def write_fasta(output_file, seq_dict): #To write a Fasta file
    with open(output_file, "w") as f:
        for name, seq in seq_dict.items():
            f.write(f">{name}\n{seq}\n")


def main():
    if len(sys.argv) < 3:
        print("Usage: python FastaParser.py GeneticData.txt output_fasta_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_fasta_file = sys.argv[2]

    mtDNA_dict, Y_dict = parse_genetic_data(input_file)

    if not mtDNA_dict and not Y_dict:
        print("Warning: No sequences found in input file.")
        sys.exit(0)

    mt_output = output_fasta_file + "_mtDNA.fasta" #To print file name as mtDNA
    y_output = output_fasta_file + "_Y.fasta" #To print file name as Y

    write_fasta(mt_output, mtDNA_dict)
    write_fasta(y_output, Y_dict)

    print(f"mtDNA sequences written to: {mt_output}")
    print(f"Y chromosome sequences written to: {y_output}")


if __name__ == "__main__":
    main()

