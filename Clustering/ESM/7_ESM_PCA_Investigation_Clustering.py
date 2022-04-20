import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
from Bio import SeqIO
import numpy as np

#Load saved PCA data from file
all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins_ESM_PCA.csv")

#Cluster dataset at 90% sequence identity to remove redundancy
accessions = [record.id for record in SeqIO.parse("Sequences/all_90percent_clustered/output.fasta", "fasta")]

all_encapsulins_clustered = all_encapsulins[np.isin(all_encapsulins["Enc Uniprot Accession"].to_numpy(), accessions)]

#Choose features
protein_vectors = all_encapsulins_clustered[[str(i) for i in range(120)]]

#Fit KMeans for  k=6
k=6
kmeans = KMeans(n_clusters=k, random_state=14).fit(protein_vectors)
all_encapsulins_clustered["Cluster"] = kmeans.labels_

#Setup plots
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(projection="3d")

#3D scatterplot of first three principal components
scatter = ax.scatter(all_encapsulins_clustered["0"], all_encapsulins_clustered["1"], all_encapsulins_clustered["2"], marker='.', c=all_encapsulins_clustered["Cluster"])
ax.set_xlabel('Principal Component 1')
ax.set_ylabel('Principal Component 2')
ax.set_zlabel('Principal Component 3')

plt.show()