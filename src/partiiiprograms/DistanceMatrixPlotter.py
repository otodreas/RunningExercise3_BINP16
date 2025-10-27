import sys
import os
import math
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

# ==============================
# User inputs and error handling
# ==============================

# Set user-assigned variables
if len(sys.argv) < 2:
    sys.exit("Too few arguments. Please pass an input file.")
elif len(sys.argv) < 3:
    input_file = sys.argv[1]
else:
    sys.exis(
        "Too many arguments. Please pass at most an input file, a chromosome name, and an output file"
    )

# Check that input_file exists
if not os.path.isfile(input_file):
    sys.exit("Input file does not exist")

# Check that input_file is a .tsv file:
if not input_file.endswith(".tsv"):
    sys.exit("Input file must be a .tsv file.")


# =============
# Program logic
# =============

def tsv_array_converter(filepath):
    data = []
    with open(filepath, "r") as f:
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
        for j in range(data_dims):
            try:
                data_array[i][j] = data[i+1][j+1]
            except ValueError:
                data_array[i][j] = 0.  # this position contains "" in the input data
    
    return data_array, labels

def custom_dendrogram(data, labels):
    fig, ax = plt.subplots(figsize=(15, 4))
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
    
    plt.tight_layout()
    plt.savefig("dendrogram.png")

def custom_heatmap(data, labels):
    for i in range(len(data)):
        for j in range(len(data)):
            data[i][j] = np.float64(math.log(data[i][j] + 1))
    fig, ax = plt.subplots(figsize=(8, 8))
    pos = ax.imshow(data)#, cmap="hsv")
    fig.colorbar(pos)
    xticks, yticks = list(labels), list(labels)
    plt.xticks(range(len(labels)), xticks, rotation=90)
    plt.yticks(range(len(labels)), yticks)
    plt.tight_layout()
    plt.savefig("heatmap.png")

data, labels = tsv_array_converter(input_file)
custom_dendrogram(data, labels)
custom_heatmap(data, labels)
