# Import libraries and set random seed
import sys
import os
import random

random.seed(0)


# ==============================
# User inputs and error handling
# ==============================

# Set user-assigned variables
if len(sys.argv) < 3:
    sys.exit(
        "Too few arguments. Please pass at least an input file and a chromosome name"
    )
elif len(sys.argv) < 4:
    input_file, chromo_name = sys.argv[1:3]
    output_file = "output_file.tsv"
elif len(sys.argv) < 5:
    input_file, chromo_name, output_file = sys.argv[1:4]
    if not sys.argv[3].endswith(".tsv"):
        sys.exit("The output file path must be a .tsv file.")
else:
    sys.exis(
        "Too many arguments. Please pass at most an input file, a chromosome name, and an output file"
    )

# Check that input_file exists
if not os.path.isfile(input_file):
    sys.exit("Input file does not exist")

# Get valid chromosome names from the user
valid_chromo_names = ["mtDNA", "Y_chromosome"]

while chromo_name not in valid_chromo_names:
    chromo_name = input(
        f"Invalid chromosome name. Please enter one of the following chromosome names: {valid_chromo_names}, or 'abort' to exit: "
    )
    if chromo_name.lower() == "abort":
        sys.exit("Program aborted.")

# Replace "_" with " " to match data
chromo_name = chromo_name.replace("_", " ")


# Check if the location of the output file is valid.
if "/" in output_file:
    final_slash = len(output_file) - output_file[::-1].find("/")
    if not os.path.exists(output_file[:final_slash]):
        sys.exit("Output path does not exist")
else:
    final_slash = 0


# =============
# Program logic
# =============

### Load data ###


def read_chromosomes(input_file, chromo_name):
    # Open input file as binary
    with open(input_file, "rb") as f:
        # Assign variables
        chromo_line = None
        i = 0
        seqs = []

        # Read lines to strings. Clean and append lines after chromosome header to data.
        while True:
            line_b = f.readline()
            line = str(line_b)

            # The chromosome data will be on the line following the chromosome name.
            if chromo_name in line:
                chromo_line = i + 1
            if i == chromo_line:
                clean_line = line[
                    2 : line.find("\\")
                ]  # since files are read in binary, the first two characters will always be 'b.
                seqs.append(clean_line)
            i += 1
            if not line_b:
                break

        return seqs


seqs = read_chromosomes(input_file, chromo_name)

# Unique values (delete later)
unique_seqs = []
for s in seqs:
    if s not in unique_seqs:
        unique_seqs.append(s)
print("\n", "Number of sequences: ", len(seqs))
print("Number of unique sequences: ", len(unique_seqs), "\n")

# Check that at least one sequence was found
if len(seqs) < 3:
    sys.exit(
        f"Error: too few (<3) chromosomes of '{chromo_name}' found. File corrupted?"
    )

# Get the lengths of each sequence, check that they are equal, assign it to the variable positions.
seq_lens = []
for seq in seqs:
    seq_lens.append(len(seq))
if len(set(seq_lens)) > 1:
    sys.exit(f"Error: sequences for '{chromo_name}' have different lengths.")
else:
    positions = seq_lens[0]


### Perform calculations and write data ###

# Assign variables
header = [
    "Chromosome",
    "Position",
    "Alleles",
    "Major Allele",
    "Minor Allele",
    "Minor Allele Frequency",
]
row_written = False
total_positions = positions * len(seqs)
total_nucleotides = 0

# Loop through positions
for p in range(positions):
    # Assign empty list to nucleotides, where valid nucleotides at position p will be stored.
    nucleotides = []

    # Loop through sequences
    for seq in seqs:
        # Append nucleotide list nucleotides
        if seq[p].upper() in "ACGT":
            nucleotides.append(seq[p].upper())

    # Throw warning if no valid nucleotides were found at position p.
    if len(nucleotides) == 0:
        print(f"Warning: no valid nucleotides found at position {p + 1}")

    # If multiple different nucleotides were found at the position, calculate statistics.
    elif nucleotides.count(nucleotides[0]) < len(nucleotides):
        # get unique nucleotides.
        alleles = set(nucleotides)

        # Store nucleotide frequencies by allele in the dictionary freqs.
        freqs = {}
        for allele in alleles:
            freqs[allele] = nucleotides.count(allele)

        # Since alleles is a set, sort the dictionary for reproducibility
        freqs = dict(sorted(freqs.items()))

        # Find major alleles, minor alleles, and minor frequency of alleles.
        # Ensure random selection if multiple major or minor alleles have the same frequency.
        major_alleles = []
        minor_alleles = []
        for allele in freqs.keys():
            if freqs[allele] == max(list(freqs.values())):
                major_alleles.append(allele)
            if freqs[allele] == min(list(freqs.values())):
                minor_alleles.append(allele)
        major_allele = random.choice(major_alleles)

        if major_allele in minor_alleles:
            minor_alleles.pop(
                minor_alleles.index(major_allele)
            )  # when all allele frequencies are equal, ensure the same allele cannot be both major and minor.
        minor_allele = random.choice(minor_alleles)
        minor_allele_freq = freqs[minor_allele] / sum(
            freqs[allele] for allele in alleles
        )
        row = (
            "\t".join(
                [
                    chromo_name,
                    str(p),
                    "/".join(alleles),
                    str(major_allele),
                    str(minor_allele),
                    str(round(minor_allele_freq, 2)),
                ]
            )
            + "\n"
        )

        # Update total_nucleotides
        total_nucleotides += sum(freqs.values())

        # Write output file, including the header if writing the first row.
        if not row_written:
            with open(output_file, "w") as f:
                header_row = "\t".join(header) + "\n" + row
                f.write(header_row)

                # Print data for inspection and update row_written to avoid re-writing the header.
                print(header_row.strip())
                row_written = True
        else:
            with open(output_file, "a") as f:
                f.write(row)

                # Print data for inspection.
                print(row.strip())

# Print location of saved file.
if final_slash:
    os.chdir(output_file[:final_slash])
print(f"Output saved to {os.path.join(os.getcwd(), output_file[final_slash:])}")

# Warn the user of low-quality data if less than 25% of positions contain valid nucleotides
pct_nucleotides = total_nucleotides / total_positions
if pct_nucleotides < 0.25:
    print(
        f"Warning: {round(pct_nucleotides * 100, 2)}% of positions in the dataset contain valid nucleotides"
    )
