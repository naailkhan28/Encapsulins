import numpy as np
from joblib import load
import matplotlib.pyplot as plt
import pandas as pd

#Load pre-fitted KMeans model
model = load("Clustering Analysis/ProtVec/kmeans_1000.joblib")

#Load existing encapsulin dataframe
all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins.csv")
all_encapsulins["Cluster"] = model.labels_

#Save all accessions in biggest cluster to file - >400 sequences
# major_cluster_seqs = all_encapsulins[all_encapsulins["Cluster"] == 769]["Enc Uniprot Accession"].tolist()

# with open("major_cluster.txt", "w") as f:
#     for seq in major_cluster_seqs:
#         f.write("Sequences/all/{}.fasta\n".format(seq))

# #Plot cluster cardinality
# k = 1000
# cardinalities = {}

# for j in range(k):
#     cardinalities[j] = np.sum(model.labels_ == j)

# plt.plot(cardinalities.keys(), cardinalities.values())
# plt.xlabel("Cluster")
# plt.ylabel("Number of Sequences")
# plt.ylim(0, 100)
# plt.show()