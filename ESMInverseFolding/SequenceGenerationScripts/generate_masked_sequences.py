import torch
import esm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pickle

#Set up cuda device
device = "cuda:1"

#Initialize model
model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
model = model.to(device)

#Import protein structure
seed_model = "data/7MU1.pdb"
seed_chain = "A"

structure = esm.inverse_folding.util.load_structure(seed_model, seed_chain)
structure = structure[:-31] #:-31 is to remove FMN heteroatoms at the end of the PDB file

coords, seq = esm.inverse_folding.util.extract_coords_from_structure(structure)
coords = torch.tensor(coords)
coords = coords.to(device)

#Mask E-loop in encapsulin structure
coords[44:75, :] = float('inf')

#Set parameters
num_sequences = 10000
T=0.5

#Generate sequences
sampled_sequences = []

for i in range(num_sequences):

    if (i+1) % 10 == 0:
        print(f"Generating sequence {i+1}/{num_sequences}")
        
    sampled_sequence = model.sample(coords, temperature=T, device=device)

    sequence_id = f"masked_T={T}_{i+1}"

    record = SeqRecord(Seq(sampled_sequence), id=sequence_id, description="")
    sampled_sequences.append(record)

outfile = SeqIO.write(sampled_sequences, f"masked_T={T}.fasta", "fasta")

#Score sequences
scores = {}

for i, record in enumerate(sampled_sequences):
    if (i+1) % 10 == 0:
        print(f"Scoring sequence {i+1}/{num_sequences}")
    scores[f"{str(record.id)}"] = esm.inverse_folding.util.score_sequence(model, alphabet, coords, str(record.seq), device)
    
outfile = pickle.dump(scores, open(f"masked_T={T}_scores.pkl", "wb"))