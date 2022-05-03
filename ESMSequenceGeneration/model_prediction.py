import torch
import torch.nn.functional as F
from sequence_generation_tools import get_token_from_residue
import sampling
import numpy as np

#Get model probabilities for all tokens at a given position in the sequence
def get_probabilities(sequence_tokens, index, model, temperature=1):
    #Disable dropout for deterministic results
    model.eval()

    results = model(sequence_tokens, repr_layers = [33])
    logits = results["logits"] / temperature

    sequence_length = logits.shape[1]
    logits = logits.reshape(sequence_length, -1)

    return(F.softmax(logits[index], dim=0))

#Get model log likelihoods for all positions of a sequence
def get_log_likelihoods(sequence_tokens, model):
    #Disable dropout for deterministic results
    model.eval()

    results = model(sequence_tokens, repr_layers = [33])
    logits = results["logits"]

    sequence_length = logits.shape[1]
    logits = logits.reshape(sequence_length, -1)

    return(logits)

#Return a list of all the masked positions in a sequence of tokens
def get_masked_positions(sequence_tokens, alphabet):

    mask_token = get_token_from_residue("<mask>", alphabet)
    mask_arr = sequence_tokens == mask_token
    mask_arr = torch.flatten(mask_arr[0].nonzero()).tolist()

    return mask_arr

#Return a list of all the unknown positions in a sequence of tokens
def get_unknown_positions(sequence_tokens, alphabet):

    unknown_token = get_token_from_residue("X", alphabet)
    mask_arr = sequence_tokens == unknown_token
    mask_arr = torch.flatten(mask_arr[0].nonzero()).tolist()

    return mask_arr

#Generate a new sequence from a masked input seed sequence
def generate_sequence_from_seed(model, alphabet, sequence_tokens, temperature=1):
    masks = get_masked_positions(sequence_tokens, alphabet)
    iteration = 1

    while masks:
        if iteration % 20 == 0:
            print("Iteration: {}, Number of mask tokens remaining: {}".format(iteration, len(masks)))
        index = np.random.choice(masks)

        probabilities = get_probabilities(sequence_tokens, index, model, temperature=temperature)

        sample = sampling.nucleus_sample(probabilities, alphabet)

        sequence_tokens[0][index] = sample

        masks = get_masked_positions(sequence_tokens, alphabet)
        iteration += 1

    return(sequence_tokens)

#Fill in any unknown "X" residues with new predictions
def fill_unknown_positions(model, alphabet, sequence_tokens):
    unknowns = get_unknown_positions(sequence_tokens, alphabet)
    iteration = 1

    while unknowns:
        index = np.random.choice(unknowns)

        probabilities = get_probabilities(sequence_tokens, index, model)

        sample = sampling.nucleus_sample(probabilities, alphabet)

        sequence_tokens[0][index] = sample

        unknowns = get_unknown_positions(sequence_tokens, alphabet)
        iteration += 1
    
    return(sequence_tokens)