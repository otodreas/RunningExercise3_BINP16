#!/usr/bin/env python3

"""

Author: Yaqi Jiao
Date: 23th October, 2025
PART 2 (Student 2) of RE - 3 with Chetana and Oliver

build_similarity_matrix.py
--------------------------

Description:
This Python script reads a MSA result and computes a genetic distance matrix among samples.
The script performs the following operations:
    1. Reads a tab-separated input file containing pairwise alignment information
    2. Corrects identity scores to account for potential alignment errors.
    3. Calculates pairwise genetic distances using either:
        - Simple P-distance
        - Jukes-Cantor (JC)
    5. Constructs a symmetric genetic distance matrix with all samples and outputs it as a tab-separated file

Input file: build_similarity_matrix.py, input_file
Output file: output_file

Procedure:
    1. Load and prepare input data
    2. Calculate basic genetic distances
    3. Construct genetic distance matrix
    4. Output results

Usage Example:
    python build_similarity_matrix.py output_MSA_file_mtDNA_scores.txt mtDNA_genetic_distance.tsv
"""

# package prepration
import sys
import numpy as np
import pandas as pd


# process MSA data, build dataframe
def load_MSA(input_file):
    rows = []
    with open(input_file, "r") as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split()  # Split each row's elements

            for i, p in enumerate(parts):  # extract names
                if "_" in p:
                    SampleA_idx = i
                    SampleA = " ".join(parts[0 : SampleA_idx + 1])
                    break

            for j in range(SampleA_idx + 1, len(parts)):
                if "_" in parts[j]:
                    SampleB_idx = j
                    SampleB = " ".join(parts[SampleA_idx + 1 : SampleB_idx + 1])
                    break

            numeric_parts = parts[
                SampleB_idx + 1 :
            ]  # Extract the rest part of the line
            comparable = int(numeric_parts[0])
            uncertain = int(numeric_parts[1])
            total = int(numeric_parts[2])
            identity = float(numeric_parts[3].rstrip("%"))
            score = int(numeric_parts[4])

            rows.append(
                [SampleA, SampleB, comparable, uncertain, total, identity, score]
            )

    df = pd.DataFrame(
        rows,
        columns=[
            "SampleA",
            "SampleB",
            "Comparable nucleotides",
            "uncertain nucleotides",
            "Total nucleotides",
            "IdentityScore",
            "Score",
        ],
    )
    return df


# Access arguments passed from the command line written in the terminal
# check we have 3 arguements: the python script, input file and output file
def main():
    if len(sys.argv) == 3:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
    else:
        sys.exit(1)

    # open and read input file

    df = load_MSA(input_file)
    print(df)

    # Identity score correction (optional)
    err = 0.1
    # Remove the percent sign from the identity score and convert it to a decimal
    df["IdentityScore"] = df["IdentityScore"].astype(float) / 100
    # The adjusted identity score (IdentityScore_tr) is computed by reducing the original identity score by
    # a fraction proportional to the error rate, and adding back a compensation term based on the remaining non-identity portion.
    # Specifically, we subtract IdentityScore * err from the original value and add (1 - IdentityScore) * err.
    df["IdentityScore_tr"] = (
        df["IdentityScore"]
        - df["IdentityScore"] * err
        + (1 - df["IdentityScore"]) * err
    )

    # Calculate similarity
    # calculate P-distance
    df["P_distance"] = 1 - df["IdentityScore_tr"]

    # process P-distance
    # Set the threshold of applying JC, then calculate the proportion of "unavailable distance" (resulting in negative log/0)
    JC_LIMIT = 0.75
    EPS = 1e-9
    prop_bad = (df["P_distance"] >= JC_LIMIT).mean()
    p = df["P_distance"]

    # If more than 50% of P-distance is unavailable -> use all original P-distance
    # Otherwise, use JC for the entire group, for the P-distance that may exceed the limit, clip it to JC_LIMIT - EPS first
    if prop_bad > 0.5:
        print(
            f"{prop_bad * 100:.1f}% of P-values >= {JC_LIMIT}, using original P-distance for all."
        )
        basic_distance = p.copy()
    else:
        print(
            f"{prop_bad * 100:.1f}% of P-values >= {JC_LIMIT}, using JC correction for all (with clipping)."
        )
        p_safe = p.clip(upper=JC_LIMIT - EPS)
        basic_distance = -(3 / 4) * np.log(1 - (4 / 3) * p_safe)

        df["genetic_distance"] = basic_distance

    # Construct the final genetic distance matrix
    labels = pd.concat(
        [df["SampleA"], df["SampleB"]]
    ).unique()  # Extrac all elements in data
    gd_df = pd.DataFrame(
        ".", index=labels, columns=labels
    )  # create an empty dataframe to store genetic distances

    # fill the matrix
    for idx, row in df.iterrows():
        s1 = row["SampleA"]
        s2 = row["SampleB"]
        s3 = row["genetic_distance"]
        # optional: Keep s3 to three decimal places
        # s3 = round(s3, 3)

        gd_df.at[s1, s2] = s3
        gd_df.at[s2, s1] = s3

    # Write the matrix to the output file
    gd_df.to_csv(output_file, sep="\t")
    print("Finish")


if __name__ == "__main__":
    main()
