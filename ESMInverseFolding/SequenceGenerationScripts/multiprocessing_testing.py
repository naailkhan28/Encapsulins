import torch
import torch.distributed as dist
import torch.multiprocessing as mp
import esm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pickle

def generate_sequences(rank, world_size):
    #Create default process group
    dist.init_process_group("gloo", rank=rank, world_size=world_size)
    
    #Initialize model
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.to(rank)

    #Import protein structure
    seed_model = "data/6NJ8.pdb"
    seed_chain = "A"

    structure = esm.inverse_folding.util.load_structure(seed_model, seed_chain)

    coords, seq = esm.inverse_folding.util.extract_coords_from_structure(structure)
    coords = torch.tensor(coords)
    coords = coords.to(rank)

    #Mask E-loop in encapsulin structure
    coords[44:75, :] = float('inf')

    #Set parameters
    num_sequences = 3333 #Per GPU
    T=0.5

    #Generate sequences
    sampled_sequences = []

    for i in range(num_sequences):

        if (i+1) % 10 == 0:
            print(f"Generating sequence {i+1}/{num_sequences} on GPU {rank}")
            
        sampled_sequence = model.sample(coords, temperature=T, device=rank)

        sequence_id = f"masked_T={T}_{i+1 +(num_sequences * (rank-1))}"

        record = SeqRecord(Seq(sampled_sequence), id=sequence_id, description="")
        sampled_sequences.append(record)

    outfile = SeqIO.write(sampled_sequences, f"masked_T={T}_GPU{rank}.fasta", "fasta")

    #Score sequences
    scores = {}

    for i, record in enumerate(sampled_sequences):
        if (i+1) % 10 == 0:
            print(f"Scoring sequence {i+1}/{num_sequences} on GPU {rank}")
        scores[f"{str(record.id)}"] = esm.inverse_folding.util.score_sequence(model, alphabet, coords, str(record.seq), rank)
        
    outfile = pickle.dump(scores, open(f"masked_T={T}_GPU{rank}_scores.pkl", "wb"))

def main():
    world_size = 4
    mp.spawn(generate_sequences,
        args=(world_size,),
        nprocs=world_size,
        join=True)

if __name__=="__main__":
    main()