#Imports
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

import matplotlib.pyplot as plt


#NOTE: This code below was run on a remote server
#Commented out here because it takes forever to run
# #Import encapsulin data
# all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins_processed_vector_representation.csv")
# print("Encapsulin data loaded successfully!")

# #Choose features
# protein_vectors = all_encapsulins.iloc[:, 1:-1]

# #Fit KMeans for various k values and plot sum of squares error against number of clusters
# k_values = [10, 50, 100, 250, 500, 750, 1000, 2000]
# sse = {}

# for k in k_values:
#     print("Running KMeans clustering with {} clusters".format(k))
#     kmeans = KMeans(n_clusters=k).fit(protein_vectors)
#     sse[k] = kmeans.inertia_


# #Write results to file for plotting or analysis
# print("Clustering done! Writing results to file...")

# strings = [f'{key} : {sse[key]}' for key in sse]
# with open('elbow.txt', 'w') as f:
#     [f.write(f'{string}\n') for string in strings]

#This is the code for the elbow plot
#The data comes from a text file generated by the code above
sse = {}

with open("Clustering Analysis/kmeans_protvec_elbowplot.txt") as f:
    for line in f:
        values = line.rstrip().split(":")

        sse[int(values[0])] = float(values[1])

plt.plot(sse.keys(), sse.values())
plt.xlabel("Number of Clusters")
plt.ylabel("Sum of Squares Error")
plt.show()
