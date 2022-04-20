#Imports
from numpy import column_stack
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import seaborn as sns

all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins.csv")

sequence_column = all_encapsulins["Sequence"]

#Functions to apply to columns
def calculate_molecular_weight(seq):
    return(ProteinAnalysis(seq).molecular_weight())

def calculate_pi(seq):
    return(ProteinAnalysis(seq).isoelectric_point())

def get_amino_acid_frequency(seq, residue):
    return(ProteinAnalysis(seq).get_amino_acids_percent()[residue])

def remove_illegal_chars(text):
    text = text.replace("B", "N")
    text = text.replace("Z", "Q")
    text = text.replace("J", "I")
    text = text.replace("X", "I")
    return text

#Clean sequences to remove illegal characters
sequence_column = sequence_column.apply(remove_illegal_chars)

print(all_encapsulins.columns)

#Add sequence feature columns
all_encapsulins["Length"] = sequence_column.str.len()
all_encapsulins["Molecular Weight"] = sequence_column.apply(calculate_molecular_weight)
all_encapsulins["pI"] = sequence_column.apply(calculate_pi)

#Add frequencies of all 20 amino acids
amino_acids = "ACDEFGHIKLMNPQRSTVWY"

for residue in amino_acids:

    column_title = residue + "_frequency"

    all_encapsulins[column_title] = sequence_column.apply(get_amino_acid_frequency, args=(residue,))


#Write dataframe to csv file to avoid having to process each time
output = all_encapsulins.to_csv("all_encapsulins_processed.csv")

print(all_encapsulins.head())