import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import more_itertools as mit

# paths = ["top10", "top50", "top100", "top1000","top100_PsiBLAST", "phmmer", "family1"]

# all_ids = []

# for path in paths:
#     dataframe = pd.read_csv(f"ESMInverseFolding/data/HMMs/scores/t4_{path}_hmm_scores.csv")
#     dataframe = dataframe.iloc[:, 1:]
#     dataframe["HMM E-Value"] = -1 * np.log10(dataframe["HMM E-Value"])
#     all_ids.extend(dataframe[dataframe["HMM E-Value"] > 14.5]["id"].to_list())

hmm_dataframe = pd.read_csv("ESMInverseFolding/data/HMMs/scores/bmc_h_hmm_scores.csv")
hmm_dataframe = hmm_dataframe.iloc[:, 1:]
hmm_dataframe["HMM E-Value"] = -1 * np.log10(hmm_dataframe["HMM E-Value"])

ids = hmm_dataframe[hmm_dataframe["HMM E-Value"].between(14, 16)]["id"].to_list()


sns.histplot(data=hmm_dataframe, x="HMM E-Value")
plt.vlines(x=14, ymin=0, ymax=120, linestyles="dashed", color="black")
plt.vlines(x=16, ymin=0, ymax=120, linestyles="dashed", color="black")
plt.show()