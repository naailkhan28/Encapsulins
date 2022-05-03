import torch
import esm
import sequence_generation_tools
import model_prediction
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
from scoring_functions import get_average_log_likelihood, get_max_sequence_similarity
import pickle
import pandas as pd

#Use GPU if available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#Initialize model and alphabet, and load batch converter for tokenizing sequences
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()
model.to(device)

#Make a list of seed sequences to choose from when generating
elongatus = "MTDNAPQLALRDVAARQLANATKTVPQLRTITPRWLVRLLHWTPVEAGIYRVNQVKDASQITVACSERDESELPETFVDYIDNPREYLLSAVNTVVDVHTRISDLYSNPHDQIREQLRLTIEIMKERQESELINSREYGLLNNVAPGQLVHTRNGAPTPDDLDELLIRVWKEPAFFLAHPQAIAAFGRECTRRGVPPATVSLFGSSFITWRGVPLIPSDKVPLENGKTKILLLRVGESRQGVVGLYQPNLPGEQGMGLSVRFMGINRKALASYLVSLYCSLAVLTDDALAVLDNVDVTQYHTYRYNSG"
maritima = "MEFLKRSFAPLTEKQWQEIDNRAREIFKTQLYGRKFVDVEGPYGWEYAAHPLGEVEVLSDENEVVKWGLRKSLPLIELRATFTLDLWELDNLERGKPNVDLSSLEETVRKVAEFEDEVIFRGCEKSGVKGLLSFEERKIECGSTPKDLLEAIVRALSIFSKDGIEGPYTLVINTDRWINFLKEEAGHYPLEKRVEECLRGGKIITTPRIEDALVVSERGGDFKLILGQDLSIGYEDREKDAVRLFITETFTFQVVNPEALILLK"
xanthus = "MPDFLGHAENPLREEEWARLNETVIQVARRSLVGRRILDIYGPLGAGVQTVPYDEFQGVSPGAVDIVGEQETAMVFTDARKFKTIPIIYKDFLLHWRDIEAARTHNMPLDVSAAAGAAALCAQQEDELIFYGDARLGYEGLMTANGRLTVPLGDWTSPGGGFQAIVEATRKLNEQGHFGPYAVVLSPRLYSQLHRIYEKTGVLEIETIRQLASDGVYQSNRLRGESGVVVSTGRENMDLAVSMDMVAAYLGASRMNHPFRVLEALLLRIKHPDAICTLEGAGATERR"
thermotolerans = "MNKSQLYPDSPLTDQDFNQLDQTVIEAARRQLVGRRFIELYGPLGRGMQSVFNDIFMESHEAKMDFQGSFDTEVESSRRVNYTIPMLYKDFVLYWRDLEQSKALDIPIDFSVAANAARDVAFLEDQMIFHGSKEFDIPGLMNVKGRLTHLIGNWYESGNAFQDIVEARNKLLEMNHNGPYALVLSPELYSLLHRVHKDTNVLEIEHVRELITAGVFQSPVLKGKSGVIVNTGRNNLDLAISEDFETAYLGEEGMNHPFRVYETVVLRIKRPAAICTLIDPEE"
seeds = [elongatus, maritima, xanthus, thermotolerans]
seed_names = ["elongatus", "maritima", "xanthus", "thermotolerans"]

#List of different masking levels to test - here from 20% to 90%
masking_levels = [x * 0.1 for x in range(2, 4)]
print(f"Masking Levels: {masking_levels}")

#Generate sequences for each masking level
records = []
num_sequences = 1000

for masking_level in masking_levels:
    print(f"Generating sequences with masking level {masking_level}")

    for i in range(num_sequences):
        print("Generating sequence {} / {}".format(i+1, num_sequences))

        seed = random.choice(seeds)
        seed_tokenized = sequence_generation_tools.tokenize(seed, batch_converter)
        seed_masked = sequence_generation_tools.mask_sequence(seed_tokenized, masking_level, alphabet)

        seed_masked = seed_masked.to(device)

        generated = model_prediction.generate_sequence_from_seed(model, alphabet, seed_masked, temperature=3)
        filled = model_prediction.fill_unknown_positions(model, alphabet, generated)

        protein_sequence = sequence_generation_tools.detokenize(filled, alphabet)

        seq_record = SeqRecord(Seq(protein_sequence[5:-5]))
        seq_record.id = f"{i+1}_{seed_names[seeds.index(seed)]}_{masking_level:.1f}"
        seq_record.description = ""

        records.append(seq_record)

#Write all generated sequences to a file - sequences will have a FASTA ID with format number_seed sequence_masking level
sequences_path = "masking_level_testing.fasta"
sequences_outfile = SeqIO.write(records, sequences_path, "fasta")


#Score all generated sequences by sequence similarity and model log likelihood
sequence_data = {}
known_encapsulins_path = "Sequences/all_90percent_clustered/output.fasta"

for i, record in enumerate(SeqIO.parse(sequences_path, "fasta")):
    print(f"Iteration {i}")
    sequence = str(record.seq)
    id = record.id
    
    similarity, accession = get_max_sequence_similarity(sequence, known_encapsulins_path)
    likelihood = get_average_log_likelihood(sequence, model, batch_converter, device)

    sequence_data[id] = (likelihood, similarity, accession)

#Write score data to a pickle file for later reference
scores_path = "masking_levels_testing_scores.pkl"
output = pickle.dump(sequence_data, open(scores_path, 'wb'))

#Load score data to a dataframe for analysis
generated_scores_dataframe = pd.DataFrame.from_dict(sequence_data, orient="index", columns = ["Model Log Likelihood", "Maximum Similarity", "Most Similar Sequence"])

#Take the top 15 and bottom 15 sequences each for similarity and model log likelihood
top15_likelihood = generated_scores_dataframe.sort_values(["Model Log Likelihood"], ascending=False)[:15].index.to_list()
bottom15_likelihood = generated_scores_dataframe.sort_values(["Model Log Likelihood"])[:15].index.to_list()
top15_similarity = generated_scores_dataframe.sort_values(["Maximum Similarity"], ascending=False)[:15].index.to_list()
bottom15_similarity = generated_scores_dataframe.sort_values(["Maximum Similarity"])[:15].index.to_list()

#Find the sequences nearest the median for model log likelihood and similarity
total_sequences = num_sequences * len(masking_levels)
middle_index = int(total_sequences / 2)
median_likelihood = generated_scores_dataframe.sort_values(["Model Log Likelihood"]).index.to_list()[middle_index]
median_similarity = generated_scores_dataframe.sort_values(["Maximum Similarity"]).index.to_list()[middle_index]

#Collect all the desired sequences ranked by scoring into a list
all_sequences = top15_likelihood + bottom15_likelihood + top15_similarity + bottom15_similarity
all_sequences.append(median_likelihood)
all_sequences.append(median_similarity)

print(all_sequences)

for sequence in SeqIO.parse(sequences_path, "fasta"):
    if sequence.id in all_sequences:
        outfile = SeqIO.write(sequence, f"{id}.fasta", "fasta")