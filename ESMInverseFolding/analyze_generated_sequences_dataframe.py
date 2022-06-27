import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scoring_functions import get_instability_index

#Load dataframe
sequences_dataframe = pd.read_csv("ESMInverseFolding/data/BMC_H_sequences_dataframe.csv")
sequences_dataframe = sequences_dataframe.iloc[:, 1:]

# upper_bounds = range(55, 20, -5)
# upper_bounds = [x / 100 for x in upper_bounds]
# lower_bounds = range(50, 15, -5)
# lower_bounds = [y / 100 for y in lower_bounds]

# ids = []

# for upper, lower in zip(upper_bounds, lower_bounds):

#     sequence_subset = sequences_dataframe[sequences_dataframe["% Identity"].between(lower, upper)]

#     sequence_subset = sequence_subset.sort_values(by="Log Likelihood Score", ascending=False)
#     accessions = sequence_subset["id"].to_list()

#     print(sequence_subset)

sns.histplot(data=sequences_dataframe, x="% Identity")
plt.show()