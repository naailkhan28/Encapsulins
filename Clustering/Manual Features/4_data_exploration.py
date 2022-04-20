#Imports
from numpy import amin
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Load previously processed dataframe

all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins_processed.csv")

amino_acids = "ACDEFGHIKLMNPQRSTVWY"

aa_frequency_column_titles = [residue + "_frequency" for residue in amino_acids]

fig, axs = plt.subplots(5, 4)

i = 0

for row in axs:
    for col in row:
        sns.histplot(all_encapsulins[aa_frequency_column_titles[i]], ax=col)
        i += 1

plt.tight_layout()       
plt.show()