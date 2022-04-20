import numpy as np

#Make a dictionary with all possible protein n-grams as keys, and their 100d vector representations as values
def get_ngrams(filename):
    ngrams = {}

    with open(filename, 'r') as vectors:
        vectors.readline()
        for line in vectors:
            columns = line.split()

            ngrams[columns[0]] = np.array([float(value) for value in columns[1:]])

    return ngrams

#Represent a protein sequence as a set of ngram vectors
def get_vector_representation(sequence, ngram_dict):

    vectors = []

    for i in range(2, len(sequence)):
        ngram = "".join(sequence[i-2:i+1])
        try:
            vectors.append(ngram_dict[ngram])
        except KeyError:
            vectors.append(ngram_dict["<unk>"])

    return np.concatenate(vectors)

#Pad the vector representation of a protein sequence to a set length
def pad_array(array, length):
    t = length - len(array)
    if t:
        return np.pad(array, pad_width=(0, t), mode='constant')
    return array
