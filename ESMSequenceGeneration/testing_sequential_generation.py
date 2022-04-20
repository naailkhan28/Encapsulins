import torch
import torch.nn.functional as F
import esm

#Generate a new protein sequence from a seed sequence using the ESM protein language transformer model
#Recursively add new predicted residues to the end of a seed sequence until a specified maximum length is reached
def generate_sequence(model, batch_converter, residue_tokens, seed_sequence, max_length, temperature=1):

    print(seed_sequence)

    if len(seed_sequence) >= max_length:
        return(seed_sequence)

    model.eval()  # disables dropout for deterministic results

    data = [("seed_sequence", seed_sequence)]
    _, _, tokens = batch_converter(data)

    tokens[0][-1] = 32
    eos = torch.tensor([[2]])
    tokens = torch.cat((tokens, eos), 1)

    results = model(tokens, repr_layers = [33])
    probabilities = F.softmax(results["logits"] / temperature, dim=1).reshape(len(seed_sequence) + 3, -1)

    sample = torch.multinomial(probabilities[-1][4:-9], 1)
    seed_sequence += residue_tokens[sample]

    return(generate_sequence(model, batch_converter, residue_tokens, seed_sequence, max_length))


#NOTES: This approach generates gibberish sequences!
#ESM-1b is a BERT model so uses full self-attention as opposed to masked self-attention
#This means that when producing a probability distribution for a residue it expects to see context either side of that residue
#For generating new proteins sequentially this doesn't work - adding new residues to the end of a sequence means that you have no existing context to the right of that residue

model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()
tokens = alphabet.all_toks[4:-9]
furiosus = "MSTRGDLIRILGEIEEKMNELKMDGFNPDIILFGREAYNFLSNLLKKEMEEEGPFTHVSNIKIEILEELGGDAVVIDSKVLGLVPGAAKRIKIIK"
maritima = "MEFLKRSFAPLTEKQWQEIDNRAREIFKTQLYGRKFVDVEGPYGWEYAAHPLGEVEVLSDENEVVKWGLRKSLPLIELRATFTLDLWELDNLERG"

generate_sequence(model, batch_converter, tokens, maritima, 250)