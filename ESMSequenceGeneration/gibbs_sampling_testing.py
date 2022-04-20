import torch
import torch.nn as nn
import esm
import sequence_generation_tools
import model_prediction
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Use GPU if available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#Initialize model and alphabet, and load batch converter for converting sequences to tokens
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()

#Use parallel model across multuple GPUs if available
if torch.cuda.device_count() > 1:
    print("Using ", torch.cuda.device_count(), "GPUs")
    model = nn.DataParallel(model)

model.to(device)

furiosus = "MSTRGDLIRILGEIEEKMNELKMDGFNPDIILFGREAYNFLSNLLKKEMEEEGPFTHVSNIKIEILEELGGDAVVIDSKVLGLVPGAAKRIKIIK"

furiosus_tokenized = sequence_generation_tools.tokenize(furiosus, batch_converter)
furiosus_tokenized = furiosus_tokenized.to(device)

furiosus_masked = sequence_generation_tools.mask_sequence(furiosus_tokenized, 0.15, alphabet)

records = []

for i in range(1000):
    print("Generating sequence {} / {}".format(i+1, 1000))
    
    furiosus_padded = sequence_generation_tools.pad_sequence(furiosus_masked, 250, alphabet)

    generated = model_prediction.generate_sequence_from_seed(model, alphabet, furiosus_padded, temperature=3)
    filled = model_prediction.fill_unknown_positions(model, alphabet, generated)

    protein_sequence = sequence_generation_tools.detokenize(filled, alphabet)

    seq_record = SeqRecord(Seq(protein_sequence[5:-5]))
    seq_record.id = str(i+1)
    seq_record.description = ""

    records.append(seq_record)

    generated = None
    filled = None
    protein_sequence = None

outfile = SeqIO.write(records, "generated_encapsulins_1000_p_furiosus_family4.fasta", "fasta")