#Imports
import pandas as pd
from sklearn.cluster import KMeans
from joblib import dump

#Import encapsulin data
print("Reading data...")
all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins_processed_vector_representation.csv")

print(all_encapsulins.head())

#Choose features
protein_vectors = all_encapsulins.iloc[:, 1:-1]

#Fit KMeans model
k=1000
kmeans = KMeans(n_clusters=k).fit(protein_vectors)

output = dump(kmeans, "kmeans_1000.joblib")