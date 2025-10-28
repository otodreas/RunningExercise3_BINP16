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
        # check that it contains no more than 4 files
        if len(os.listdir(sys.argv[1])) <= 4:
            # create empty list input_files to store filepaths in 
            input_files = [""] * len(os.listdir(sys.argv[1]))

            # loop through filepaths in the directory passed and add them to input_files
            for i, f in enumerate(os.listdir(sys.argv[1])):
                input_files[i] = os.path.join(sys.argv[1], f)
        else:
            sys.exit(f"Too many files in {sys.argv[1]}")
    # one file was passed to the program
    else:
        input_files = [sys.argv[1]]

# multiple files were passed. 
else:
    input_files = [""] * len(sys.argv)
    for i in range(len(input_files)):
        input_files[i] = sys.argv[i]
    print(input_files)
    sys.exit()
    

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


# =============
# Program logic
# =============

# Create object class Matrix_File
class Matrix_File:
    def __init__(self, filepath):
        self.filepath = filepath
    def get_name(self):

        # The valid chromosomes
        chromosomes = ["mtDNA", "YDNA"]
        chromosome = None
        # The valid metrics
        metrics = ["alignment", "identity"]
        metric = None
        for i in chromosomes:
            if i in self.filepath:
                chromosome = i
        for i in metrics:
            if i in self.filepath:
                metric = i

        # Check that names are valid
        while chromosome not in chromosomes:
            chromosome = input(f"Chromosome not found. Please input a valid chromosome name {chromosomes} or 'abort' to exit: ")
            if chromosome.lower() == "abort":
                sys.exit("Program aborted.")
        while metric not in metrics:
            metric = input(f"Metric not found. Please input a valid metric name {metrics} or 'abort' to exit: ")
            if metric.lower() == "abort":
                sys.exit("Program aborted.")
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
            if data[i+1][0] != labels[i]:
                sys.exit("Matrix axes are unaligned. Perhaps the Y axis is inverted?")
            for j in range(data_dims):
                try:
                    data_array[i][j] = data[i+1][j+1]
                except ValueError:
                    data_array[i][j] = 0.0  # this position contains "" in the input data
        
        return data_array, labels


def custom_dendrogram(data, labels, chromosome, metric):
    fig, ax = plt.subplots(figsize=(12, 4))
    #c_data = squareform(data)
    Z = linkage(data)
    dendrogram(Z, orientation="left", labels=labels)
    
    family = np.copy(labels)
    romanovs = np.array(["Olga", "Tatiana", "Marie", "Anastasia", "Alexandra", "Nicolas", "Romanov"])
    for i, ind in enumerate(family):
        for j, romanov in enumerate(romanovs):
            if romanov in ind:
                family[i] = "Romanov"
                break
            elif j == len(romanovs) - 1:
                family[i] = "Other" # if iteration over romanovs is finished, update family
            else:
                pass
    
    ax.set_title(f"Relatedness dendrogram on chromosome {chromosome}. Metric: {metric}")
    plt.tight_layout()
    plt.savefig(f"dendrogram_{chromosome}_{metric}.png")

def custom_heatmap(data, labels, chromosome, metric):
    for i in range(len(data)):
        for j in range(len(data)):
            data[i][j] = np.float64(math.log(data[i][j] + 1))
    
    # unique values (delete later)
    print("\n", "Number of pairwise comparisons in matrix: ", int(len(data.flatten()) / 2 - len(data)))
    print("Number of unique values in matrix: ", len(np.unique(data.flatten())), "\n")

    
    # """
    # Create a heatmap from a numpy array and two lists of labels.

    # Parameters
    # ----------
    # data
    #     A 2D numpy array of shape (M, N).
    # row_labels
    #     A list or array of length M with the labels for the rows.
    # col_labels
    #     A list or array of length N with the labels for the columns.
    # ax
    #     A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
    #     not provided, use current Axes or create a new one.  Optional.
    # cbar_kw
    #     A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    # cbarlabel
    #     The label for the colorbar.  Optional.
    # """
        

    # # Plot the heatmap
    # im = ax.imshow(data)
    
    # # Create colorbar
    # cbar = ax.figure.colorbar(im, ax)
    # cbar.ax.set_ylabel("Log relatedness", rotation=-90, va="bottom")

    # xticks, yticks = list(labels), list(labels)
    # plt.xticks(range(len(labels)), xticks, rotation=90)
    # plt.yticks(range(len(labels)), yticks)
    
    # ax.spines[:].set_visible(False)


    fig, ax = plt.subplots(figsize=(8, 8))
    pos = ax.imshow(data)
    fig.colorbar(pos)
    xticks, yticks = list(labels), list(labels)
    plt.xticks(range(len(labels)), xticks, rotation=90)
    plt.yticks(range(len(labels)), yticks)
    ax.set_title(f"Relatedness heatmap on chromosome {chromosome}.\nMetric: {metric}")
    fig.tight_layout()
    plt.savefig(f"heatmap_{chromosome}_{metric}.png")

# Run functions
if __name__ == "__main__":
    for input_file in input_files:
        fileobject = Matrix_File(input_file)
        chromosome, metric = fileobject.get_name()
        data, labels = fileobject.get_data()
        custom_dendrogram(data, labels, chromosome, metric)
        custom_heatmap(data, labels, chromosome, metric)
