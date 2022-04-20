import torch

#Take the most probable token from the distribution
def greedy_search(probabilities):
    return(torch.argmax(probabilities).item())

#Sample from the distribution of tokens
def simple_sample(probabilities):
    return(torch.multinomial(probabilities, 1))[0]

#Sample from the top K most likely tokens
def topk_sample(probabilities, k):
    values, indices = torch.topk(probabilities, k)
    values = values * (1 / torch.sum(values))
    return(indices[torch.multinomial(values, 1)])[0]

#Nucleus sample from the distribution
def nucleus_sample(probabilities, alphabet, p=0.92):
    sorted_probabilities, indices = torch.sort(probabilities, descending=True)

    cumulative_probability = 0
    max_index = 0

    for i in range(len(sorted_probabilities)):
        if cumulative_probability > 0.92:
            break
        cumulative_probability += sorted_probabilities[i]
        max_index += 1

    sample_nucleus = sorted_probabilities[:max_index+1]
    sample_nucleus = sample_nucleus * (1 / torch.sum(sample_nucleus))

    sample = torch.multinomial(sample_nucleus, 1)

    residue = indices[sample][0].item()

    return(residue)