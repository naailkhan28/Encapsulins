import numpy as np
from joblib import load
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
import pickle

#Load pre-fitted KMeans model
model = load("Clustering Analysis/ESM/kmeans_2000.joblib")

#Load existing encapsulin dataframe
all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins.csv")
all_encapsulins["Cluster"] = model.labels_

#Load representations from pickle - keys are accession and values are representations as numpy arrays
representations_dict = pickle.load(open("ESM/encapsulin_embeddings.pkl", 'rb'))

#Load dict into a dataframe for convenience
embedded_dataframe = pd.DataFrame.from_dict(representations_dict, orient="index")
protein_vectors = embedded_dataframe.iloc[:, 1:-1]

# num_pca_components = 120
# pca = PCA(num_pca_components)
# protein_vectors_PCA = pca.fit_transform(protein_vectors)

# fig_dims = (7, 6)
# fig, ax = plt.subplots(figsize=fig_dims)
# sc = ax.scatter(protein_vectors_PCA[:,0], protein_vectors_PCA[:,1], protein_vectors_PCA[:,2], marker='.')
# ax.set_xlabel('PCA first principal component')
# ax.set_ylabel('PCA second principal component')
# plt.show()

#Plot cluster cardinality
k = 2000

cardinalities = {}
for j in range(k):
    cardinalities[j] = np.sum(model.labels_ == j)

plt.bar(cardinalities.keys(), cardinalities.values())
plt.xlabel("Cluster")
plt.ylabel("Number of Sequences")
#plt.ylim(10)
plt.show()