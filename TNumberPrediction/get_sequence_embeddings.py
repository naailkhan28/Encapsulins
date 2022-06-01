import pandas as pd
import esm
import torch
from sklearn.decomposition import PCA
import numpy as np
import pickle


#Read saved CSV with sequences
encapsulin_data = pd.read_csv("TNumberPrediction/data/all_seqs.csv")

#Prepare sequences for ESM model
accessions = encapsulin_data["Accession"].to_list()
sequences = encapsulin_data["Sequence"].to_list()

#Pick number of sequences to extract representations for - will remove this after testing is complete
protein_data = [(accession, sequence) for accession, sequence in zip(accessions, sequences)]

#Load ESM-1b model
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

#Tokenize sequences
labels, seqs, tokens = batch_converter(protein_data)

#Extract per-residue representations (on CPU)
with torch.no_grad():
    results = model(tokens, repr_layers=[33])
token_representations = results["representations"][33]

#Generate per-sequence representations via averaging
sequence_representations = []
for i, (_, seq) in enumerate(protein_data):
    sequence_representations.append(token_representations[i, 1 : len(seq) + 1].mean(0).numpy())

sequence_representations = np.stack(sequence_representations)

#PCA to reduce size of embeddings
num_pca_components = 120
pca = PCA(num_pca_components)
reduced_embeddings = pca.fit_transform(sequence_representations)

#Write reduced embeddings to a file
embeddings_outfile = pickle.dump(reduced_embeddings, open("encapsulin_t_number_embeddings.pkl", 'wb'))