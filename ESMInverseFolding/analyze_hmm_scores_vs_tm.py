import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

names = {"top10": "Top 10 BLAST Hits", "top50": "Top 50 BLAST Hits", "top100": "Top 100 BLAST Hits", "top1000": "Top 1000 BLAST Hits",
            "top100_PsiBLAST": "Top 100 Psi-BLAST Hits", "phmmer": "phmmer Search Hits", "family1": "All Family 1 Encapsulins"}


tm_scores_1 = pd.read_csv("ESMInverseFolding/data/hmm_tests_tm_scores_hi.csv")
tm_scores_1 = tm_scores_1.iloc[:, 1:]

tm_scores_2 = pd.read_csv("ESMInverseFolding/data/hmm_tests_tm_scores.csv")
tm_scores_2 = tm_scores_2.iloc[:, 1:]

tm_scores_all = pd.concat([tm_scores_1, tm_scores_2])

fig, axs = plt.subplots(2, 4)

for i, (path, title) in enumerate(names.items()):
    hmm_scores = pd.read_csv(f"ESMInverseFolding/data/HMMs/scores/t4_{path}_hmm_scores.csv")
    hmm_scores = hmm_scores.iloc[:, 1:]


    all_scores = pd.merge(left = hmm_scores, right = tm_scores_all, on="id")

    all_scores["HMM E-Value"] = -1 * np.log10(all_scores["HMM E-Value"])

    sns.scatterplot(data=all_scores, x="HMM E-Value", y="T4", ax=axs.flat[i])
    axs.flat[i].vlines(x=14.5, ymin=0, ymax=1, linestyles="dashed", color="black")
    axs.flat[i].hlines(y=0.5, xmin=0, xmax=20, linestyles="dashed", color="black")
    axs.flat[i].set_title(title)
    axs.flat[i].set_ylabel("TM Score")

plt.tight_layout()
plt.show()