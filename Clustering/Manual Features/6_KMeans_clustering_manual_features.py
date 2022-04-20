#Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import PowerTransformer, QuantileTransformer
from sklearn.cluster import KMeans

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import euclidean

from collections import defaultdict

#Import encapsulin data
all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins_processed.csv")
all_encapsulins = all_encapsulins.drop(["Unnamed: 0"], axis=1)

#Choose features
quantile_features = ['Length', 'Molecular Weight', 'pI', 'A_frequency',
       'D_frequency', 'E_frequency', 'F_frequency',
       'G_frequency', 'H_frequency', 'I_frequency', 'K_frequency',
       'L_frequency', 'M_frequency', 'N_frequency', 'P_frequency',
       'Q_frequency', 'R_frequency', 'S_frequency', 'T_frequency',
       'V_frequency', 'W_frequency', 'Y_frequency']

log_features = ['C_frequency']

categoric_features = ['Cargo Type', 'Family']


#Set up preprocessor for feature normalization
preprocessor = ColumnTransformer(transformers=[

        ("quantile", QuantileTransformer(output_distribution="normal"), quantile_features),
        ("log", PowerTransformer(), log_features)

])

all_encapsulins_scaled = preprocessor.fit_transform(all_encapsulins)



#Fit KMeans model for k = 100
#Initial assumption of optimal k value for quality checking

k=100

kmeans = KMeans(n_clusters=k).fit(all_encapsulins_scaled)

all_encapsulins["Cluster"] = kmeans.labels_

# #Setup figure and plots

# fig, axs = plt.subplots(1, 3)

# #Plot cluster cardinality
# cardinalities = {}
# for j in range(k):
#     cardinalities[j] = np.sum(kmeans.labels_ == j)

# axs[0].bar(cardinalities.keys(), cardinalities.values())
# axs[0].set_xlabel("Cluster")
# axs[0].set_ylabel("Number of Sequences")

# #Plot cluster magnitude - sum of distances of all observations in each cluster to the centroid

# magnitudes = defaultdict(float)

# for i in range(len(all_encapsulins_scaled)):
#     cluster = kmeans.labels_[i]

#     magnitudes[cluster] += np.abs(euclidean(all_encapsulins_scaled[i], kmeans.cluster_centers_[cluster]))

# axs[1].bar(magnitudes.keys(), magnitudes.values())
# axs[1].set_xlabel("Cluster")
# axs[1].set_ylabel("Magnitude")

# #Plot cluster cardinality versus cluster magnitude
# axs[2].scatter(cardinalities.values(), magnitudes.values())
# axs[2].set_xlabel("Cluster Cardinality")
# axs[2].set_ylabel("Cluster Magnitude")

# plt.show()

plt.hist(all_encapsulins["Length"], edgecolor="black")
plt.xlabel("Length")
plt.ylabel("No. of Sequences")
plt.show()