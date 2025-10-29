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

elif len(sys.argv) < 3:
    # check if the sole argument is a directory
    if os.path.isdir(sys.argv[1]):
        tsv_files = [f for f in os.listdir(sys.argv[1]) if f.endswith(".tsv")]
        # check that it contains no more than 4 files
        if len(tsv_files) <= 4:
            # create empty list input_files to store filepaths in
            input_files = [""] * len(tsv_files)

            # loop through filepaths in the directory passed and add them to input_files
            for i, f in enumerate(tsv_files):
                input_files[i] = os.path.join(sys.argv[1], f)
        else:
            sys.exit(f"Too many files in {sys.argv[1]}")
    # one file was passed to the program
    else:
        input_files = [sys.argv[1]]

# multiple files were passed.
else:
    input_files = [""] * (len(sys.argv) - 1)
    for i in range(len(input_files)):
        input_files[i] = sys.argv[i + 1]


for f in input_files:
    nonexistent_files = []
    nontsv_files = []
    if not os.path.isfile(f):
        nonexistent_files.append(f)

    # Check that input_file is a .tsv file:
    if not f.endswith(".tsv"):
        nontsv_files.append(f)
if len(nonexistent_files) > 0:
    sys.exit(f"The following files do not exist: {nonexistent_files}")
if len(nontsv_files) > 0:
    sys.exit(f"The following files are not tsv files: {nontsv_files}")


# =====================
# Classes and functions
# =====================

# User can type "Abort" at any time and abort the program
def abort(user_input):
    if user_input.lower() == "abort" or user_input.lower() == "exit":
        sys.exit("Program aborted.")
    else:
        return user_input


class Matrix_File_Metadata:
    def __init__(self, input_files):
        self.input_files = input_files

    def __str__(self):
        return str(self.input_files)

    def get_metadata(self):
        # user input interface
        chromosomes = []
        metrics = []
        # Loop through files
        for i in self.input_files:
            # Get chromosome name if no chromosomes (force input since you cannot iterate through empty list)
            if len(chromosomes) == 0:
                chromosomes.append(
                    abort(
                        input(f"Please input the CHROMOSOME name in the file '{i}': ")
                    )
                )
            else:
                chromosome_already_labeled = False
                # Loop through chromosomes
                for c in chromosomes:
                    if c in i:
                        chromosome_already_labeled = True # update chromosome_already_labeled
                if chromosome_already_labeled:
                    pass
                else:
                    chromosomes.append(
                        abort(
                            input(
                                f"Please input the CHROMOSOME name in the file '{i}': "
                            )
                        )
                    )

            if len(metrics) == 0:
                metrics.append(
                    abort(input(f"Please input the METRIC name in the file '{i}': "))
                )
            else:
                metric_already_labeled = False
                for m in metrics:
                    if m in i:
                        metric_already_labeled = True
                if metric_already_labeled:
                    pass
                else:
                    metrics.append(
                        abort(
                            input(f"Please input the METRIC name in the file '{i}': ")
                        )
                    )

            # Generate dictionaries to return
            chromo_counts = {}
            metric_counts = {}
            for c in chromosomes:
                chromo_counts[c] = 0
            for m in metrics:
                metric_counts[m] = 0
                
            for i in self.input_files:
                for c in chromosomes:
                    if c in i:
                        chromo_counts[c] += 1
                for m in metrics:
                    if m in i:
                        metric_counts[m] += 1
                        

        return (chromo_counts, metric_counts)


# Create object class Matrix_File
class Matrix_File:
    def __init__(self, filepath):
        self.filepath = filepath

    def get_name(self, chromosomes, metrics):
        chromosome = []
        metric = []

        for i in chromosomes.keys():
            if i in self.filepath:
                chromosome.append(i)
        for i in metrics.keys():
            if i in self.filepath:
                metric.append(i)

        if len(chromosome) + len(metric) != 2:
            sys.exit(
                f"File {self.filepath} does not contain exactly one chromosome and one metric name."
            )
        else:
            chromosome = "".join(chromosome)
            metric = "".join(metric)

        self.chromosome = chromosome
        self.metric = metric

        return chromosome, metric

    def get_data(self):
        data = []
        with open(self.filepath, "r") as f:
            while True:
                line = f.readline()
                row = line.strip("\n").split("\t")
                if line:
                    data.append(row)
                else:
                    break

        row_lengths = {len(row) for row in data}
        if len(row_lengths) > 1 or row_lengths == {1}:
            sys.exit(f"Delimiter issues. Rows read have length(s) {row_lengths}")
        labels = np.array(data[0][1:])

        data_dims = len(data) - 1
        data_array = np.zeros((data_dims, data_dims))

        for i in range(data_dims):
            if data[i + 1][0] != labels[i]:
                sys.exit("Matrix axes are unaligned. Perhaps the Y axis is inverted?")
            for j in range(data_dims):
                try:
                    data_array[i][j] = data[i + 1][j + 1]
                except ValueError:
                    data_array[i][j] = (
                        0.0  # this position contains "" in the input data
                    )
        
        chromosome_tag = self.chromosome
        for i, label in enumerate(labels):
            label_terms = label.split("_")
            label_terms_lower = label.lower().split("_")
            while chromosome_tag not in label_terms:
                chromosome_tag = str(
                    abort(
                        input(
                            f"The matrix row/column labels do not match the chromosome tag '{chromosome_tag}' in the file name '{self.filepath}'. From the matrix data label '{label}', what characters represent the chromosome name? "
                        )
                    )
                )

            label_terms.pop(label_terms.index(chromosome_tag))
            label_trim = "_".join(label_terms)
            labels[i] = label_trim

        return data_array, labels


def custom_dendrogram(input_files):
    
    fig, axs = plt.subplots(1, len(input_files), figsize=(10 + 5 * (len(input_files) - 1), 5))
    orientations = ["left", "right"]
    orientation_flipper = 0
    for i, input_file in enumerate(input_files):
        ax = axs[i] if len(input_files) > 1 else axs
        
        fileobject = Matrix_File(input_file)
        chromosome, metric = fileobject.get_name(chromo_counts, metric_counts)
        data, labels = fileobject.get_data()
        
        c_data = squareform(data)
        Z = linkage(c_data)
        dendrogram(Z, orientation=orientations[orientation_flipper], labels=labels, ax=ax)
        orientation_flipper = not orientation_flipper
    
        # Color labels
        romanovs = np.array(
            ["Olga", "Tatiana", "Marie", "Anastasia", "Alexandra", "Nicolas", "Romanov"]
        )
        for j, plot_label in enumerate(ax.get_yticklabels()):
            for k, romanov in enumerate(romanovs):
                if romanov in plot_label.get_text():
                    plot_label.set_color('blue')
                    print(plot_label, "blue")
                    break
                elif k < len(romanovs) - 1:
                    pass
                else:
                    plot_label.set_color('red')
                    print(plot_label, "red")
        
        ax.set_title(f"Metric: {metric}")
    fig.suptitle(f"Hierarchical Clustering of Genetic Distances (Algorithm: nearest point)\nChromosome: {chromosome}")
    plt.tight_layout()
    plt.savefig(f"dendrogram_{chromosome}.png")


def custom_heatmap(input_files):
    

    fig, axs = plt.subplots(1, len(input_files), figsize=(10 + 5 * (len(input_files) - 1), 5))
    for i, input_file in enumerate(input_files):
        ax = axs[i] if len(input_files) > 1 else axs

        fileobject = Matrix_File(input_file)
        chromosome, metric = fileobject.get_name(chromo_counts, metric_counts)
        data, labels = fileobject.get_data()

        max_char_len = 18
        for j, label in enumerate(labels):
            if len(label) > max_char_len:
                labels[j] = "".join(["...", label[len(label) - max_char_len: ]])
    
        data = np.negative(data)
        im = ax.imshow(data, cmap="nipy_spectral")
        if i == len(input_files) - 1:
            cbar = ax.figure.colorbar(im)
        xticks, yticks = list(labels), list(labels)
        ax.set_xticks(range(len(labels)), xticks, rotation=90)
        ax.set_yticks(range(len(labels)), yticks)
        ax.set_title(f"Metric: {metric}")
    fig.suptitle(f"Negative Genetic Distance Heatmap\nChromosome: {chromosome}")
    fig.tight_layout()
    plt.savefig(f"heatmap_{chromosome}.png")


# Run functions
if __name__ == "__main__":
    input_files_metadata = Matrix_File_Metadata(input_files)
    chromo_counts, metric_counts = input_files_metadata.get_metadata()

    for chromosome in chromo_counts.keys():
        chromo_input_files = []
        for input_file in input_files:
            if chromosome in input_file:
                chromo_input_files.append(input_file)

        custom_dendrogram(chromo_input_files)
        custom_heatmap(chromo_input_files)
