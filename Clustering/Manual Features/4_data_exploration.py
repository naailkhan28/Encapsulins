#Imports
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Load previously processed dataframe
all_encapsulins = pd.read_csv("Clustering/all_encapsulins_processed_manualfeatures.csv")

notfamily1 = all_encapsulins[all_encapsulins["Family"] != "1"]

notfamily1 = notfamily1[~notfamily1["Length"].between(250, 300)]

ids = notfamily1["Enc Uniprot Accession"].to_list()
sequences = notfamily1["Sequence"].to_list()

records = []

for id, sequence in zip(ids, sequences):
    id = f"natural_{id}"
    record = SeqRecord(Seq(str(sequence)), id=id, description="")

    records.append(record)

outfile = SeqIO.write(records, "DiscriminatorModel/data/sequences/validation_natural_sequences.fasta", "fasta")