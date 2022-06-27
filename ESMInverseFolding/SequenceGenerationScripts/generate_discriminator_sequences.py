import torch
import esm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np

#Initialize model
rank = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model, _ = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
model = model.to(rank)
    
#Set parameters
seed_models = ["data/7MU1.pdb", "data/6NJ8.pdb", "data/7BOJ.pdb", "data/7OE2.pdb", "data/7S21.pdb"]
seed_chain = "A"
num_sequences = 2235
temperatures = [1e-06, 0.1, 0.5, 1, 2]

#Generate sequences
sampled_sequences = []

for i in range(num_sequences):

    #Randomly choose a seed structure and temperature value from the list
    seed_model = np.random.choice(seed_models)
    print(seed_model)
    T=np.random.choice(temperatures)

    seed_structure = esm.inverse_folding.util.load_structure(seed_model, seed_chain)

    coords, _ = esm.inverse_folding.util.extract_coords_from_structure(seed_structure)
    coords = torch.tensor(coords)
    coords = coords.to(rank)

    if (i+1) % 10 == 0:
        print(f"Generating sequence {i+1}/{num_sequences}")

    sampled_sequence = model.sample(coords, temperature=T, device=rank)

    sequence_id = f"generated_{seed_model}_{i+1}"

    record = SeqRecord(Seq(sampled_sequence), id=sequence_id, description="")
    sampled_sequences.append(record)

outfile = SeqIO.write(sampled_sequences, f"data/discriminator_generated.fasta", "fasta")