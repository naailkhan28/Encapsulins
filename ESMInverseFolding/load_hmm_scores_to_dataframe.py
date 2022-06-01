import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

path = "ESMInverseFolding/data/HMMs/scores/q_thermotolerans_t4_generated/t4_phmmer/"
all_files=os.listdir(path)

hmm_scores = {}

for i, file in enumerate(all_files):

    with open(f"{path}/{file}", "r") as f:
        f.readline()
        f.readline()
        f.readline()

        data = f.readline().split()
        try:
            hmm_scores[i] = (data[2], float(data[4]), float(data[5]))
        except IndexError:
            hmm_scores[i] = (file[:-10], 1.0, 0.0)

hmm_dataframe = pd.DataFrame.from_dict(hmm_scores, orient="index", columns=["id", "HMM E-Value", "HMM Score"])

print(hmm_dataframe)
outfile = hmm_dataframe.to_csv("ESMInverseFolding/data/HMMs/scores/t4_phmmer_hmm_scores.csv")

sns.histplot(data=hmm_dataframe, x="HMM E-Value")
plt.show()