#Imports
import pandas as pd
from get_protvec_representations import get_ngrams, get_vector_representation, pad_array

#Import encapsulin data
all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins.csv")

all_encapsulins = all_encapsulins[["Enc Uniprot Accession", "Sequence"]]

#Get vector representation for each sequence
ngrams = get_ngrams("ProtVec/protVec_100d_3grams.csv")
all_encapsulins["Vector"] = all_encapsulins["Sequence"].apply(get_vector_representation, args=(ngrams,))

#Pad all vector representations with zeros to the length of the longest sequence
max_length = int(all_encapsulins["Sequence"].str.len().max() - 2) * 100
all_encapsulins["Padded Vector"] = all_encapsulins["Vector"].apply(pad_array, args=(max_length,))

#Drop unnecessary columns
all_encapsulins = all_encapsulins.drop(["Sequence", "Vector"], axis=1)

#Split one column with numpy arrays into a set of columns, each with their own value
all_encapsulins_processed = all_encapsulins["Padded Vector"].apply(pd.Series)
all_encapsulins_processed["Accession"] = all_encapsulins["Enc Uniprot Accession"]

#Write the processed DataFrame to CSV for downstream clustering
outfile = all_encapsulins_processed.to_csv("all_encapsulins_processed_vector_representation.csv")