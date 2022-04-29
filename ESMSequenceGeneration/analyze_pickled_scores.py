import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

all_generated_scores = {}

for i in range(4):
    filename = f"ESMSequenceGeneration/Scores/generated_6k_v2_{i}.pkl"

    current_score_dictionary = pickle.load(open(filename, "rb"))

    for key, value in current_score_dictionary.items():
        all_generated_scores[int(key) + (15000 * i)] = value

generated_scores_dataframe = pd.DataFrame.from_dict(all_generated_scores, orient="index", columns = ["Model Log Likelihood", "Maximum Similarity", "Most Similar Sequence"])

top10_likelihood = generated_scores_dataframe.sort_values(["Model Log Likelihood"], ascending=False)[:15].index.to_list()
bottom10_likelihood = generated_scores_dataframe.sort_values(["Model Log Likelihood"])[:15].index.to_list()
top10_similarity = generated_scores_dataframe.sort_values(["Maximum Similarity"], ascending=False)[:15].index.to_list()
bottom10_similarity = generated_scores_dataframe.sort_values(["Maximum Similarity"])[:15].index.to_list()

median_likelihood = generated_scores_dataframe.sort_values(["Model Log Likelihood"]).index.to_list()[30000]
median_similarity = generated_scores_dataframe.sort_values(["Maximum Similarity"]).index.to_list()[30000]

all_sequences = top10_likelihood + bottom10_likelihood + top10_similarity + bottom10_similarity

all_sequences.append(median_likelihood)
all_sequences.append(median_similarity)