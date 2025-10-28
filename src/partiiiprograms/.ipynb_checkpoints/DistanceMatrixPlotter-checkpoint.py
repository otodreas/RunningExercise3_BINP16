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
    print(input_files)
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
        chromosomes = []
        metrics = []
        for i in self.input_files:
            if len(chromosomes) == 0:
                chromosomes.append(
                    abort(
                        input(f"Please input the CHROMOSOME name in the file '{i}': ")
                    )
                )
            else:
                chromosome_already_labeled = False
                for c in chromosomes:
                    if c in i:
                        chromosome_already_labeled = True
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

        return (chromosomes, metrics)


# Create object class Matrix_File
class Matrix_File:
    def __init__(self, filepath):
        self.filepath = filepath

    def get_name(self, chromosomes, metrics):
        chromosome = []
        metric = []

        for i in chromosomes:
            if i in self.filepath:
                chromosome.append(i)
        for i in metrics:
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


def custom_dendrogram(data, labels, chromosome, metric):
    fig, ax = plt.subplots(figsize=(10, 4))
    # c_data = squareform(data)
    Z = linkage(data)
    dendrogram(Z, orientation="left", labels=labels)

    family = np.copy(labels)
    romanovs = np.array(
        ["Olga", "Tatiana", "Marie", "Anastasia", "Alexandra", "Nicolas", "Romanov"]
    )
    for i, ind in enumerate(family):
        for j, romanov in enumerate(romanovs):
            if romanov in ind:
                family[i] = "Romanov"
                break
            elif j < len(romanovs) - 1:
                pass
            else:
                family[i] = (
                    "Other"  # if iteration over romanovs is finished, update family
                )

    ax.set_title(f"Relatedness Dendrogram on Chromosome {chromosome}. Metric: {metric}")
    plt.tight_layout()
    plt.savefig(f"dendrogram_{chromosome}_{metric}.png")


def custom_heatmap(data, labels, chromosome, metric):
    for i in range(len(data)):
        for j in range(len(data)):
            data[i][j] = np.float64(math.log(data[i][j] + 1))

    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(data)
    cbar = ax.figure.colorbar(im)  # , panchor=(0., 0.))
    # cbar.ax.set_ylabel("Log relatedness", rotation=-90)
    # fig.colorbar(pos)
    xticks, yticks = list(labels), list(labels)
    plt.xticks(range(len(labels)), xticks, rotation=90)
    plt.yticks(range(len(labels)), yticks)
    ax.set_title(f"Log Relatedness Heatmap\nChromosome: {chromosome}\nMetric: {metric}")
    fig.tight_layout()
    plt.savefig(f"heatmap_{chromosome}_{metric}.png")


# Run functions
if __name__ == "__main__":
    input_files_metadata = Matrix_File_Metadata(input_files)
    chromosomes, metrics = input_files_metadata.get_metadata()
    
    for input_file in input_files:
        fileobject = Matrix_File(input_file)
        chromosome, metric = fileobject.get_name(chromosomes, metrics)
        data, labels = fileobject.get_data()
        custom_dendrogram(data, labels, chromosome, metric)
        custom_heatmap(data, labels, chromosome, metric)
