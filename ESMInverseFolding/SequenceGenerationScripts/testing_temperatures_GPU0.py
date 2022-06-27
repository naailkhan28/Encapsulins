import torch
import esm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pickle

#Set up cuda device
device = "cuda:0"

#Initialize model
model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
model = model.to(device)

#Import protein structure
seed_model = "7MU1.pdb"
seed_chain = "A"

structure = esm.inverse_folding.util.load_structure(seed_model, seed_chain)
structure = structure[:-31] #:-31 is to remove FMN heteroatoms at the end of the PDB file

coords, seq = esm.inverse_folding.util.extract_coords_from_structure(structure)
coords = torch.tensor(coords)
coords = coords.to(device)

#Set parameters
num_sequences = 1000
temperatures = [1e-6, 1e-4]

#Generate sequences
sampled_sequences = {}

for T in temperatures:
    print(f"Temperature = {T}")
    sampled_sequences[f"{T}"] = []

    for i in range(num_sequences):

        if (i+1) % 10 == 0:
            print(f"Generating sequence {i+1}/{num_sequences}")
        
        sampled_sequence = model.sample(coords, temperature=T, device=device)

        sequence_id = f"{T}_{i+1}"

        record = SeqRecord(Seq(sampled_sequence), id=sequence_id, description="")
        sampled_sequences[f"{T}"].append(record)

    outfile = SeqIO.write(sampled_sequences[f"{T}"], f"{T}.fasta", "fasta")

#Score sequences
for T, records in sampled_sequences.items():
    scores = {}

    for i, record in enumerate(records):
        if (i+1) % 10 == 0:
            print(f"Scoring sequence {i+1}/{num_sequences}")
        scores[f"str(record.id)"] = esm.inverse_folding.util.score_sequence(model, alphabet, coords, str(record.seq), device)
    
    outfile = pickle.dump(scores, open(f"{T}.pkl", "wb"))