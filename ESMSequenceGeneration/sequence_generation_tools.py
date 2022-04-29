import torch

#Turn a protein sequence into a sequence of tokens for input into the ESM model
def tokenize(sequence, batch_converter):
    data = [("sequence", sequence)]
    _, _, tokens = batch_converter(data)

    return(tokens)

#Get a protein sequence back from tokens output from the ESM model
def detokenize(sequence_tokens, alphabet):
    all_tokens = alphabet.all_toks
    sequence_tokens = sequence_tokens[0].tolist()

    out_sequence = [all_tokens[token] for token in sequence_tokens]

    return "".join(out_sequence)

#Get the token of a given residue or special character
def get_token_from_residue(residue, alphabet):
    tokens = alphabet.all_toks

    return(tokens.index(residue))

#Randomly replace a given fraction of tokens in a tensor with a masking token
def mask_sequence(sequence_tokens, masked_fraction, alphabet):
    masking_token = get_token_from_residue("<mask>", alphabet)
    start_token = get_token_from_residue("<cls>", alphabet)
    end_token = get_token_from_residue("<eos>", alphabet)

    #Create a random array of floats with equal dimensions to input tensor
    rand = torch.rand(sequence_tokens.shape)

    #Create a masking array, ensuring not to mask out the beginning or end of sequence tokens
    mask_arr = (rand < masked_fraction) * (sequence_tokens != start_token) * (sequence_tokens != end_token)
    mask_arr = torch.flatten(mask_arr[0].nonzero()).tolist()

    #Apply the mask to the input tensor
    for i in mask_arr:
        sequence_tokens[0][i] = masking_token

    return(sequence_tokens)

#Pad a tensor with masking tokens to a maximum length
def pad_sequence(sequence_tokens, max_length, alphabet):
    masking_token = get_token_from_residue("<mask>", alphabet)
    end_token = get_token_from_residue("<eos>", alphabet)

    masking_token_tensor = torch.tensor([[masking_token]])
    end_token_tensor = torch.tensor([[end_token]])

    sequence_tokens[0][-1] = masking_token

    for i in range(len(sequence_tokens[0] -1), max_length + 1):
        sequence_tokens = torch.cat((sequence_tokens, masking_token_tensor), 1)

    sequence_tokens = torch.cat((sequence_tokens, end_token_tensor), 1)

    return(sequence_tokens)

#Remove illegal or rare characters from amino acid sequences
def remove_illegal_chars(text):
    text = text.replace("B", "N")
    text = text.replace("Z", "Q")
    text = text.replace("J", "I")
    text = text.replace("X", "I")
    return text