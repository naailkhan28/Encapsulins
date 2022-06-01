import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

#Load dataframe
sequences_dataframe = pd.read_csv("ESMInverseFolding/data/masked_eloop_sequences.csv")
sequences_dataframe = sequences_dataframe.iloc[:, 1:]


#Bin sequences by % identity
upper_bounds = range(500, 300, -25)
upper_bounds = [x / 1000 for x in upper_bounds]
lower_bounds = range(475, 275, -25)
lower_bounds = [y / 1000 for y in lower_bounds]

ids = []

for upper, lower in zip(upper_bounds, lower_bounds):

    sequence_subset = sequences_dataframe[sequences_dataframe["% Identity"].between(lower, upper)]

    sequence_subset = sequence_subset.sort_values(by="Log Likelihood Score", ascending=False)
    accessions = sequence_subset["id"].to_list()

    #Get top 4 and bottom 4 sequences ranked by log likelihood
    accessions = accessions[:4] + accessions[-4:]

    ids += accessions

#Write each sequence to its own FASTA file
# for record in SeqIO.parse("ESMInverseFolding/data/generated_seqs/masked_eloop.fasta", "fasta"):

#     if str(record.id) in ids[:-2]:
#         outfile = SeqIO.write(record, f"ESMInverseFolding/data/structure_predictions/masked_eloop_tests/{str(record.id)}.fasta", "fasta")
#         print(f"{str(record.id)}.fasta")

#Plot % identity against log likelihood with AF2 predicted sequences highlighted
# def get_predicted(id, ids):
#     return id in ids
# sequences_dataframe["Predicted?"] = sequences_dataframe["id"].apply(get_predicted, args=([ids]))

# sns.scatterplot(data=sequences_dataframe, x="% Identity", y="Log Likelihood Score", hue="Predicted?")
# plt.legend([],[], frameon=False)
# plt.show()