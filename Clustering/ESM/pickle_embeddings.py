import pickle
import torch
import os

embeddings = {}

for file in os.listdir():
    if file.endswith(".pt"):
        embedding = torch.load(file)

    accession = embedding["label"]
    representation = embedding["mean_representations"][33].numpy()
    
    embeddings[accession] = representation


output = pickle.dump(embeddings, open("encapsulin_embeddings.pkl", 'wb'))
