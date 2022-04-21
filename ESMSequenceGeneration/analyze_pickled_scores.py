import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

scores = pickle.load(open("ESMSequenceGeneration/generated_encapsulin_scores.pkl", "rb"))

scores_dataframe = pd.DataFrame.from_dict(scores, orient="index", columns = ["Model Log Likelihood", "Maximum Similarity", "Most Similar Sequence"])

likelihood_mean = scores_dataframe["Model Log Likelihood"].mean()
similarity_mean = scores_dataframe["Maximum Similarity"].mean()

likelihood_std = scores_dataframe["Model Log Likelihood"].std()
similarity_std = scores_dataframe["Maximum Similarity"].std()

max_likelihood = scores_dataframe[scores_dataframe["Model Log Likelihood"] > likelihood_mean + likelihood_std]
min_similarity = max_likelihood[max_likelihood["Maximum Similarity"] < similarity_mean - similarity_std]

sequences = min_similarity.index.to_list()

# for record in SeqIO.parse("Sequences/generated_1000sequences_pfuriosus_family4_seed.fasta", "fasta"):
#     if int(record.id) in sequences:
#         seq_record = SeqRecord(record.seq, record.id, description="")

#         outfile = SeqIO.write(seq_record, str(record.id) + ".fasta", "fasta")
#         print(str(record.id) + ".fasta")
fig, ax = plt.subplots(1)
sns.scatterplot(x=scores_dataframe["Model Log Likelihood"], y=scores_dataframe["Maximum Similarity"], ax=ax, s=10)
ax.vlines(x=likelihood_mean + likelihood_std, ymin=scores_dataframe["Maximum Similarity"].min(), ymax=scores_dataframe["Maximum Similarity"].max())
ax.hlines(y=similarity_mean - similarity_std, xmin=scores_dataframe["Model Log Likelihood"].min(), xmax=scores_dataframe["Model Log Likelihood"].max())

plt.show()