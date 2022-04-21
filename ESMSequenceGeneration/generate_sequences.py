import torch
import torch.nn as nn
import esm
import sequence_generation_tools
import model_prediction
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random

#Use GPU if available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#Initialize model and alphabet, and load batch converter for converting sequences to tokens
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()

#Use parallel model across multiple GPUs if available
if torch.cuda.device_count() > 1:
    print("Using ", torch.cuda.device_count(), "GPUs")
    model = nn.DataParallel(model)

model.to(device)


elongatus = "MTDNAPQLALRDVAARQLANATKTVPQLRTITPRWLVRLLHWTPVEAGIYRVNQVKDASQITVACSERDESELPETFVDYIDNPREYLLSAVNTVVDVHTRISDLYSNPHDQIREQLRLTIEIMKERQESELINSREYGLLNNVAPGQLVHTRNGAPTPDDLDELLIRVWKEPAFFLAHPQAIAAFGRECTRRGVPPATVSLFGSSFITWRGVPLIPSDKVPLENGKTKILLLRVGESRQGVVGLYQPNLPGEQGMGLSVRFMGINRKALASYLVSLYCSLAVLTDDALAVLDNVDVTQYHTYRYNSG"
maritima = "MEFLKRSFAPLTEKQWQEIDNRAREIFKTQLYGRKFVDVEGPYGWEYAAHPLGEVEVLSDENEVVKWGLRKSLPLIELRATFTLDLWELDNLERGKPNVDLSSLEETVRKVAEFEDEVIFRGCEKSGVKGLLSFEERKIECGSTPKDLLEAIVRALSIFSKDGIEGPYTLVINTDRWINFLKEEAGHYPLEKRVEECLRGGKIITTPRIEDALVVSERGGDFKLILGQDLSIGYEDREKDAVRLFITETFTFQVVNPEALILLK"
xanthus = "MPDFLGHAENPLREEEWARLNETVIQVARRSLVGRRILDIYGPLGAGVQTVPYDEFQGVSPGAVDIVGEQETAMVFTDARKFKTIPIIYKDFLLHWRDIEAARTHNMPLDVSAAAGAAALCAQQEDELIFYGDARLGYEGLMTANGRLTVPLGDWTSPGGGFQAIVEATRKLNEQGHFGPYAVVLSPRLYSQLHRIYEKTGVLEIETIRQLASDGVYQSNRLRGESGVVVSTGRENMDLAVSMDMVAAYLGASRMNHPFRVLEALLLRIKHPDAICTLEGAGATERR"
thermotolerans = "MNKSQLYPDSPLTDQDFNQLDQTVIEAARRQLVGRRFIELYGPLGRGMQSVFNDIFMESHEAKMDFQGSFDTEVESSRRVNYTIPMLYKDFVLYWRDLEQSKALDIPIDFSVAANAARDVAFLEDQMIFHGSKEFDIPGLMNVKGRLTHLIGNWYESGNAFQDIVEARNKLLEMNHNGPYALVLSPELYSLLHRVHKDTNVLEIEHVRELITAGVFQSPVLKGKSGVIVNTGRNNLDLAISEDFETAYLGEEGMNHPFRVYETVVLRIKRPAAICTLIDPEE"

seeds = [elongatus, maritima, xanthus, thermotolerans]

records = []
num_sequences = 100000

for i in range(num_sequences):
    print("Generating sequence {} / {}".format(i+1, num_sequences))

    seed = random.choice(seeds)
    seed_tokenized = sequence_generation_tools.tokenize(seed, batch_converter)
    seed_masked = sequence_generation_tools.mask_sequence(seed_tokenized, 0.9, alphabet)

    seed_masked = seed_masked.to(device)

    generated = model_prediction.generate_sequence_from_seed(model, alphabet, seed_masked, temperature=3)
    filled = model_prediction.fill_unknown_positions(model, alphabet, generated)

    protein_sequence = sequence_generation_tools.detokenize(filled, alphabet)

    seq_record = SeqRecord(Seq(protein_sequence[5:-5]))
    seq_record.id = str(i+1)
    seq_record.description = ""

    records.append(seq_record)

    generated = None
    filled = None
    protein_sequence = None

outfile = SeqIO.write(records, "generated_encapsulins_100000_v2.fasta", "fasta")