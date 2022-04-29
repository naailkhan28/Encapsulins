from scoring_functions import get_average_log_likelihood, get_max_sequence_similarity
from Bio import SeqIO
import esm
import pickle
import torch

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()

model = model.to(device)

sequence_data = {}

for i, record in enumerate(SeqIO.parse("Sequences/generated_1000sequences_pfuriosus_family4_seed.fasta", "fasta")):
    print(f"Iteration {i}")
    sequence = str(record.seq)
    id = int(record.id)
    
    similarity, accession = get_max_sequence_similarity(sequence, "Sequences/all_90percent_clustered/output.fasta")
    likelihood = get_average_log_likelihood(sequence, model, batch_converter, device)

    sequence_data[id] = (likelihood, similarity, accession)

output = pickle.dump(sequence_data, open("generated_encapsulin_scores.pkl", 'wb'))