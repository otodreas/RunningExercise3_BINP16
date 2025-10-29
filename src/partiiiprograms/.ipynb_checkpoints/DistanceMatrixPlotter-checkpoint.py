# Import libraries
import sys
import os
import math
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt


# ==============================
# User inputs and error handling
# ==============================

# Set user-assigned variables
if len(sys.argv) > 5:
    sys.exis("Too many arguments. Please pass at most four input files.")

# One argument passed:
elif len(sys.argv) < 3:
    # Check if the sole argument is a directory.
    if os.path.isdir(sys.argv[1]):
        # Get .tsv files
        tsv_files = [f for f in os.listdir(sys.argv[1]) if f.endswith(".tsv")]

        # Check that the directory contains no more than 4 .tsv files.
        if len(tsv_files) <= 4:
            # Create empty list input_files to store filepaths in.
            input_files = [""] * len(tsv_files)

            # Loop through the filepaths in the directory passed and add them to input_files
            for i, f in enumerate(tsv_files):
                input_files[i] = os.path.join(sys.argv[1], f)
        else:
            sys.exit(f"Too many files in {sys.argv[1]}")
    # One file was passed to the program:
    else:
        input_files = [sys.argv[1]]

# Multiple files were passed:
else:
    input_files = [""] * (len(sys.argv) - 1)
    for i in range(len(input_files)):
        input_files[i] = sys.argv[i + 1]

# Loop through input_files and check for invalid files.
for f in input_files:
    nonexistent_files = []
    nontsv_files = []
    if not os.path.isfile(f):
        nonexistent_files.append(f)

    # Check that input_file is a .tsv file:
    if not f.endswith(".tsv"):
        nontsv_files.append(f)

# Print errors.
if len(nonexistent_files) > 0:
    sys.exit(
        f"The following files do not exist or are directories: {nonexistent_files}"
    )
if len(nontsv_files) > 0:
    sys.exit(f"The following files are not tsv files: {nontsv_files}")


# =====================
# Classes and functions
# =====================


# User can type "Abort" at any time and abort the program.
def abort(user_input):
    if user_input.lower() == "abort" or user_input.lower() == "exit":
        sys.exit("Program aborted.")
    else:
        return user_input


# Class for list of input files.
class Matrix_File_Metadata:
    def __init__(self, input_files):
        self.input_files = input_files

    # When printing, return input_files
    def __str__(self):
        return str(self.input_files)

    # Getting metadata (chromosome counts, metric counts)
    def get_metadata(self):
        # user input interface
        chromosomes = []
        metrics = []
        metadata = [chromosomes, metrics]
        metadata_input_tags = ["CHROMOSOME", "METRIC"]

        # Loop through files.
        for input_file in self.input_files:
            # Loop through chromosomes and metrics
            for j, tags in enumerate(metadata):
                # Get tag name from user if none have been passed (force input since you cannot iterate through empty list)
                if len(tags) == 0:
                    tags.append(
                        abort(
                            input(
                                f"Please input the {metadata_input_tags[j]} name in the file '{input_file}': "
                            )
                        )
                    )
                # If tags list contains at least 1 entry:
                else:
                    tag_already_labeled = False
                    # Loop through tags
                    for tag in tags:
                        # If tag is already represented in input file, update tag_already_labeled
                        if tag in input_file:
                            tag_already_labeled = True
                    if tag_already_labeled:
                        pass
                    # Prompt user to input tag if it is not represented in list of tags
                    else:
                        tags.append(
                            abort(
                                input(
                                    f"Please input the {metadata_input_tags[j]} name in the file '{input_file}': "
                                )
                            )
                        )

            # Generate dictionaries to return
            chromo_counts = {}
            metric_counts = {}

            # Loop through chromosomes and metrics, initializing dictionaries with 0
            for c in chromosomes:
                chromo_counts[c] = 0
            for m in metrics:
                metric_counts[m] = 0

            # Loop through input files and get the number of occurrences of metadata
            for input_file in self.input_files:
                for c in chromosomes:
                    if c in input_file:
                        chromo_counts[c] += 1
                for m in metrics:
                    if m in input_file:
                        metric_counts[m] += 1

        return (chromo_counts, metric_counts)


# Create object class Matrix_File for each file.
class Matrix_File:
    def __init__(self, filepath):
        self.filepath = filepath

    # Get the chromosome and metric for each file
    def get_name(self, chromosomes, metrics):
        """
        Get the chromosome and metric names for a Matrix_File
        given the dictionaries 'chromosomes' and 'metrics' containing
        counts of each category of metadata in a Matrix_File_Metadata
        object.
        """

        chromosome = []
        metric = []

        # Loop through dictionary keys and grab the correct key
        for c in chromosomes.keys():
            if c in self.filepath:
                chromosome.append(c)
        for m in metrics.keys():
            if m in self.filepath:
                metric.append(m)

        # Check that the file contains exactly one entry of each metadata type.
        if len(chromosome) + len(metric) != 2:
            sys.exit(
                f"File {self.filepath} does not contain exactly one chromosome and one metric name."
            )
        else:
            # Convert metadata from list to string if it passes the check above.
            chromosome = "".join(chromosome)
            metric = "".join(metric)

        self.chromosome = chromosome
        self.metric = metric

        return chromosome, metric

    def get_data(self):
        """
        Extract numerical and categorical data from .tsv file of type Matrix_File.
        Return data in separate numpy arrays. Replace missing numerical values with 0.
        """
        # Open filepath associated with object and store lines as lists in list data
        data = []
        with open(self.filepath, "r") as f:
            while True:
                line = f.readline()
                row = line.strip("\n").split("\t")
                if line:
                    data.append(row)
                else:
                    break

        # Ensure that there are at least 2 columns and that every row has the same number of columns
        row_lengths = {len(row) for row in data}
        if len(row_lengths) > 1 or row_lengths == {1}:
            sys.exit(f"Delimiter issues. Rows read have length(s) {row_lengths}")

        # Set labels variable to the headers of the input file, not including the first empty value
        labels = np.array(data[0][1:])

        # data_dims is assigned to the length of the input data, not including the categorical headers
        data_dims = len(data) - 1

        # create an array of zeroes with dimensions of the input data
        data_array = np.zeros((data_dims, data_dims))
        data_array_debug = np.copy(data_array)

        # Loop through columns of numerical data
        for i in range(data_dims):
            # Check that the vertical label matches the horizontal label. An "upside-down" matrix will not pass.
            if data[i + 1][0] != labels[i]:
                sys.exit("Matrix axes are unaligned. Perhaps the Y axis is inverted?")
            # Loop through rows of numerical data
            non_diagonal_zeros = 0
            for j in range(data_dims):
                # Update value if the input value is numerical
                try:
                    data_array[i][j] = data[i + 1][
                        j + 1
                    ]  # account for the header rows by increasing indecies by 1.
                    # Diagonals must have a score of 0.
                    if data_array[i][j] > 0 and i == j:
                        sys.exit(
                            f"The identity matrix {self.filepath} does not contain a 0 value at diagonal position {i + 1, j + 1}"
                        )

                # Non numerical value handling
                except ValueError:
                    # The only acceptable non-numerical value is "." and "".
                    if data[i + 1][j + 1] != "." and data[i + 1][j + 1] != "":
                        sys.exit(
                            f"The identity matrix {self.filepath} contains the non-numerical value '{data[i + 1][j + 1]}' at position {i + 1, j + 1}"
                        )
                    # Warn the user when non-diagonals have a score of 0. Limit the number of 0 values to 25% of the matrix.
                    if i != j:
                        non_diagonal_zeros += 1
                        if non_diagonal_zeros > data_dims * data_dims / 4:
                            sys.exit(
                                f"The identity matrix {self.filepath} contains more than 25% empty values."
                            )
                        print(
                            f"Warning: the identity matrix {self.filepath} contains an empty value at non-diagonal position {i + 1, j + 1}"
                        )
                    pass  # leave non-numerical values in the diagonal as 0.

        # Access the chromosome of the Matrix_File
        chromosome_tag = self.chromosome

        # Loop through labels, removing chromosome tags if labels are delimited with "_".
        for i, label in enumerate(labels):
            if "_" in label:
                label_terms = label.split("_")
                label_terms_lower = label.lower().split("_")
                while chromosome_tag not in label_terms:
                    chromosome_tag = str(
                        abort(
                            input(
                                f"The matrix row/column labels do not match the chromosome tag '{chromosome_tag}' in the file name '{self.filepath}'. "
                                f"From the matrix data label '{label}', what characters represent the chromosome name? "
                            )
                        )
                    )

                label_terms.pop(label_terms.index(chromosome_tag))
                label_trim = "_".join(label_terms)
                labels[i] = label_trim

        return data_array, labels


# Define function to create dendrogram plots with multiple input files
def custom_dendrogram(input_files):
    # Set subplots based on input files passed
    fig, axs = plt.subplots(
        1, len(input_files), figsize=(10 + 5 * (len(input_files) - 1), 5)
    )
    # Initialize orientations and orientation flipper to mirror plots
    orientations = ["left", "right"]
    orientation_flipper = 0
    # Loop through input files, creating subplots
    for i, input_file in enumerate(input_files):
        # Handle single inputs, when axes objects are not iterable
        ax = axs[i] if len(input_files) > 1 else axs

        # Get custom attributes from input file
        fileobject = Matrix_File(input_file)
        chromosome, metric = fileobject.get_name(chromo_counts, metric_counts)
        data, labels = fileobject.get_data()

        # Flatten 2D input matrix to 1D without updating values using squareform
        c_data = squareform(data)
        # Calculate pairwise euclidean distances between points using the nearest point algorithm
        Z = linkage(c_data)
        # Create dendrogram, update flipper
        dendrogram(
            Z, orientation=orientations[orientation_flipper], labels=labels, ax=ax
        )
        orientation_flipper = not orientation_flipper

        # Color labels by family
        romanovs = np.array(
            ["Romanov", "Olga", "Tatiana", "Marie", "Anastasia", "Alexandra", "Nicolas"]
        )
        # Loop through y tick labels, assign colors
        for j, plot_label in enumerate(ax.get_yticklabels()):
            for k, romanov in enumerate(romanovs):
                # If the label matches any entry in the array romanovs, set the label to blue.
                if romanov.lower() in plot_label.get_text().lower():
                    plot_label.set_color("blue")
                    print(plot_label, "blue")
                    break
                elif k < len(romanovs) - 1:
                    pass
                # If all romanov names have been checked, set the label to red.
                else:
                    plot_label.set_color("red")
                    print(plot_label, "red")

        # Set title of the subplot
        ax.set_title(f"Metric: {metric}")
    # Set the title of the full plot and save.
    fig.suptitle(
        f"Hierarchical Clustering of Genetic Distances (Algorithm: nearest point)\nChromosome: {chromosome}"
    )
    plt.tight_layout()
    plt.savefig(f"dendrogram_{chromosome}.png")


# Define function to create heatplots with multiple input files
def custom_heatmap(input_files):
    # Set subplots based on input files passed
    fig, axs = plt.subplots(
        1, len(input_files), figsize=(10 + 5 * (len(input_files) - 1), 5)
    )

    # Loop through input files, creating subplots
    for i, input_file in enumerate(input_files):
        # Handle single inputs, when axes objects are not iterable
        ax = axs[i] if len(input_files) > 1 else axs

        # Get custom attributes from input file
        fileobject = Matrix_File(input_file)
        chromosome, metric = fileobject.get_name(chromo_counts, metric_counts)
        data, labels = fileobject.get_data()

        # Trim labels if they are too long
        max_char_len = 18
        for j, label in enumerate(labels):
            if len(label) > max_char_len:
                labels[j] = "".join(["...", label[len(label) - max_char_len :]])

        # Flip sign of data so that more closely related pairs have higher values
        data = np.negative(data)
        # Plot image
        im = ax.imshow(data, cmap="nipy_spectral")
        # Create a color bar only once.
        if i == len(input_files) - 1:
            cbar = ax.figure.colorbar(im)

        # Update tick labels
        xticks, yticks = list(labels), list(labels)
        ax.set_xticks(range(len(labels)), xticks, rotation=90)
        ax.set_yticks(range(len(labels)), yticks)

        # Set subtitle
        ax.set_title(f"Metric: {metric}")

    # Set title and save
    fig.suptitle(f"Negative Genetic Distance Heatmap\nChromosome: {chromosome}")
    fig.tight_layout()
    plt.savefig(f"heatmap_{chromosome}.png")


# =============
# Program logic
# =============
    
# Run functions
if __name__ == "__main__":
    # Get metadata
    input_files_metadata = Matrix_File_Metadata(input_files)
    chromo_counts, metric_counts = input_files_metadata.get_metadata()

    # Loop through chromosomes, plotting data that share chromosomes
    for chromosome in chromo_counts.keys():
        chromo_input_files = []
        for input_file in input_files:
            if chromosome in input_file:
                chromo_input_files.append(input_file)

        custom_dendrogram(chromo_input_files)
        custom_heatmap(chromo_input_files)
