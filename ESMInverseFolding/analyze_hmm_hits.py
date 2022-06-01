import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

# paths = ["top10", "top50", "top100", "top1000","top100_PsiBLAST", "phmmer", "family1"]

# all_ids = []

# for path in paths:
#     dataframe = pd.read_csv(f"ESMInverseFolding/data/HMMs/scores/t4_{path}_hmm_scores.csv")
#     dataframe = dataframe.iloc[:, 1:]
#     dataframe["HMM E-Value"] = -1 * np.log10(dataframe["HMM E-Value"])
#     all_ids.extend(dataframe[dataframe["HMM E-Value"] > 14.5]["id"].to_list())

