import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from Bio import SeqIO

tm_scores = pd.read_csv("ESMInverseFolding/data/BMC_H_TM_scores.csv")
tm_scores = tm_scores.iloc[:, 1:]
hmm_scores = pd.read_csv(f"ESMInverseFolding/data/HMMs/scores/bmc_h_hmm_scores.csv")
hmm_scores["HMM E-Value"] = -1 * np.log10(hmm_scores["HMM E-Value"])

scaler = StandardScaler()
hmm_scores["Scaled Score"] = scaler.fit_transform(hmm_scores[["HMM E-Value"]])

all_scores = pd.merge(left = hmm_scores, right = tm_scores, on="id")

ids = hmm_scores[hmm_scores["Scaled Score"] < -2]["id"].to_list()

sns.scatterplot(data=all_scores, x="Scaled Score", y="TM Score")
plt.hlines(y=0.5, xmin=-5, xmax=3, linestyles="dashed", color="black")
plt.vlines(x=-3, ymin=0, ymax=1, linestyles="dashed", color="black")
plt.title("BMC-H Generated Sequences")

plt.ylabel("TM Score")
plt.show()