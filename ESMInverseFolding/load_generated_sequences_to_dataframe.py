import pickle
import pandas as pd
from Bio import SeqIO
from scoring_functions import get_sequence_identity

#Import likelihood scores
gpus = [0,1,2,3]
all_scores = {}

for gpu in gpus:
    scores = pickle.load(open(f"ESMInverseFolding/data/scores/BMC_H_GPU{gpu}_scores.pkl", "rb"))
    all_scores.update(scores)

#Load sequences and import into a dataframe
filename = f"ESMInverseFolding/data/generated_seqs/BMC_H_generated_sequences.fasta"

all_records = {}
for record in SeqIO.parse(filename, "fasta"):
    all_records[str(record.id)] = str(record.seq)
    output = SeqIO.write(record, f"ESMInverseFolding/data/generated_seqs/BMC_H_seqs_individual/{str(record.id)}.fasta", "fasta")

sequences_dataframe = pd.DataFrame.from_dict(all_records, orient="index", columns=["Sequence"])
sequences_dataframe = sequences_dataframe.reset_index(level=0)
sequences_dataframe = sequences_dataframe.rename(columns={"index": "id", "0": "Sequence"})

#Load log likelihood scores from pickle
def load_score(id):
    return(all_scores[id][0])

sequences_dataframe["Log Likelihood Score"] = sequences_dataframe["id"].apply(load_score)

#Add sequence identity to the seed sequence
seed_sequence = "MGDALGLIETKGLVACIEAADAMCKAANVELIGYENVGSGLVTAMVKGDVGAVKAAVDSGVESAQRIGEVVTSLVIARPHNDISKIVAHYKIAE"

sequences_dataframe["% Identity"] = sequences_dataframe["Sequence"].apply(get_sequence_identity, args=(seed_sequence,))

#Add column for temperature value
sequences_dataframe["Temperature"] = 2

outfile = sequences_dataframe.to_csv("ESMInverseFolding/data/BMC_H_sequences_dataframe.csv")