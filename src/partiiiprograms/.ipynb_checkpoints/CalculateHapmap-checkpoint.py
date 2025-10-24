# Import libraries
import sys
import os

# ------------------------------
# User inputs and error handling
# ------------------------------

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
while chromo_name not in ["mtDNA", "Y chromosome"]:
    chromo_name = input(
        "Invalid chromosome name. Please enter one of the following chromosome names: 'mtDNA', 'Y chromosome', or 'abort' to exit: "
    )
    if chromo_name.lower() == "abort":
        sys.exit("Program aborted.")

# Check if the location of the output file is valid.
if "/" in output_file:
    final_slash = len(output_file) - output_file[::-1].find("/")
    if not os.path.exists(output_file[:final_slash]):
        sys.exit("Output path does not exist")
else:
    final_slash = 0


# -------------
# Program logic
# -------------

# Open input file as binary
with open(input_file, "rb") as f:
    # Assign variables
    chromo_line = None
    i = 0
    data = []

    # Read lines to strings. Clean and append lines after chromosome header to data.
    while True:
        line_b = f.readline()
        line = str(line_b)
        if chromo_name in line:
            chromo_line = i + 1
        if i == chromo_line:
            clean_line = line[
                2 : line.find("\\")
            ]  # since files are read in binary, the first two characters will always be 'b.
            data.append(clean_line)
        i += 1
        if not line_b:
            break

# Check that at least one sequence was found
if len(data) < 1:
    sys.exit(f"Error: {chromo_name} not found. File corrupted?")

# Get the lengths of each sequence, check that they are equal, assign it to the variable positions.
seq_lens = []
for seq in data:
    seq_lens.append(len(seq))
if len(set(seq_lens)) > 1:
    sys.exit(f"Error: sequences for {chromo_name} have different lengths.")
else:
    positions = seq_lens[0]

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
total_positions = positions * len(data)
total_nucleotides = 0

# Loop through positions
for p in range(positions):
    # Assign empty list to nucleotides, where valid nucleotides at position p will be stored.
    nucleotides = []

    # Loop through sequences
    for seq in data:
        # Append nucleotide list nucleotides
        if seq[p].upper() in "ACGT":
            nucleotides.append(seq[p].upper())

    # Throw warning if no valid nucleotides were found at position p.
    if len(nucleotides) == 0:
        print(f"Warning: no valid nucleotides found at position {p + 1}")

    # If more than one nucleotide were found at the position, calculate statistics.
    elif nucleotides.count(nucleotides[0]) < len(nucleotides):
        alleles = set(nucleotides)
        freqs = {}

        # Loop through
        for a in alleles:
            freqs[a] = nucleotides.count(a)
        total_nucleotides += sum(freqs.values())
        ma = max(freqs, key=freqs.get)
        mi = min(freqs, key=freqs.get)
        maf = freqs[mi] / sum(freqs[a] for a in alleles)
        row = (
            "\t".join(
                [
                    chromo_name,
                    str(p),
                    "/".join(alleles),
                    str(ma),
                    str(mi),
                    str(round(maf, 2)),
                ]
            )
            + "\n"
        )

        if not row_written:
            with open(output_file, "w") as f:
                header_row = "\t".join(header) + "\n" + row
                f.write(header_row)
                print(header_row.strip())
                row_written = True
        else:
            with open(output_file, "a") as f:
                f.write(row)
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
