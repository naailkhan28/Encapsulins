import pickle
import pandas as pd
from Bio import SeqIO
from scoring_functions import get_sequence_identity

#Import likelihood scores
gpus = [0,1,2,3]
all_scores = {}

for gpu in gpus:
    scores = pickle.load(open(f"ESMInverseFolding/data/scores/T4_GPU{gpu}_scores.pkl", "rb"))
    all_scores.update(scores)

#Load sequences and import into a dataframe
filename = f"ESMInverseFolding/data/generated_seqs/T4_generated_seqs.fasta"

all_records = {}
for record in SeqIO.parse(filename, "fasta"):
    all_records[str(record.id)] = str(record.seq)

sequences_dataframe = pd.DataFrame.from_dict(all_records, orient="index", columns=["Sequence"])
sequences_dataframe = sequences_dataframe.reset_index(level=0)
sequences_dataframe = sequences_dataframe.rename(columns={"index": "id", "0": "Sequence"})

#Load log likelihood scores from pickle
def load_score(id):
    return(all_scores[id][0])

sequences_dataframe["Log Likelihood Score"] = sequences_dataframe["id"].apply(load_score)

#Add sequence identity to the seed sequence
seed_sequence = "MNKSQLYPDSPLTDQDFNQLDQTVIEAARRQLVGRRFIELYGPLGRGMQSVFNDIFMESHEAKMDFQGSFDTEVESSRRVNYTIPMLYKDFVLYWRDLEQSKALDIPIDFSVAANAARDVAFLEDQMIFHGSKEFDIPGLMNVKGRLTHLIGNWYESGNAFQDIVEARNKLLEMNHNGPYALVLSPELYSLLHRVHKDTNVLEIEHVRELITAGVFQSPVLKGKSGVIVNTGRNNLDLAISEDFETAYLGEEGMNHPFRVYETVVLRIKRPAAICTLIDPEE"

sequences_dataframe["% Identity"] = sequences_dataframe["Sequence"].apply(get_sequence_identity, args=(seed_sequence,))

#Add column for temperature value
sequences_dataframe["Temperature"] = 1

outfile = sequences_dataframe.to_csv("ESMInverseFolding/data/T4_generated_sequences.csv")