import pandas as pd
from sklearn.cluster import KMeans
from joblib import dump
import pickle

#Load representations from pickle - keys are accession and values are representations as numpy arrays
representations_dict = pickle.load(open("encapsulin_embeddings.pkl", 'rb'))

#Load dict into a dataframe for convenience
embedded_dataframe = pd.DataFrame.from_dict(representations_dict, orient="index")

#Choose features
protein_vectors = embedded_dataframe.iloc[:, 1:-1]

#Fit KMeans model
k=2000
kmeans = KMeans(n_clusters=k).fit(protein_vectors)

#Save model file for later
output = dump(kmeans, "kmeans_2000.joblib")