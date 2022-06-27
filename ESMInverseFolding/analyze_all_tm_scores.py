import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
from scoring_functions import extract_secondary_structure, get_sequence_identity, get_path_to_ss_file, get_max_identity_parallel

#Load TM score dataframes
tm_dataframe_1 = pd.read_csv("ESMInverseFolding/data/temperature_tests_TM_scores.csv")
tm_dataframe_1 = tm_dataframe_1.iloc[:, 1:]

tm_dataframe_2 = pd.read_csv("ESMInverseFolding/data/masking_tests_TM_scores.csv")
tm_dataframe_2 = tm_dataframe_2.iloc[:, 1:]

tm_dataframe_all = pd.concat([tm_dataframe_1, tm_dataframe_2])

#Load sequence dataframes
temperature_seqs = pd.read_csv("ESMInverseFolding/data/temperature_tests_seqs_dataframe.csv")
temperature_seqs = temperature_seqs.iloc[:, 1:]
masked_seqs = pd.read_csv("ESMInverseFolding/data/masked_eloop_sequences_dataframe.csv")
masked_seqs = masked_seqs.iloc[:, 1:]

sequences_dataframe = pd.concat([temperature_seqs, masked_seqs])

#Merge both dataframes into a single dataframe
ids = tm_dataframe_all["id"].to_list()
sequences_dataframe = sequences_dataframe.query("id in @ids")
tm_dataframe = pd.merge(left=sequences_dataframe, right=tm_dataframe_all, on="id")

#Add secondary structure sequences
tm_dataframe["ss_paths"] = tm_dataframe["id"].apply(get_path_to_ss_file)
tm_dataframe["Secondary Structure"] = tm_dataframe["ss_paths"].apply(extract_secondary_structure)
tm_dataframe = tm_dataframe.drop(["ss_paths"], axis=1)

#Calculate secondary structure accuracy
t1_secondary_structure = extract_secondary_structure("ESMInverseFolding/data/dssp_secondary_structures/7mu1.txt")
tm_dataframe["Secondary Structure Accuracy"] = tm_dataframe["Secondary Structure"].apply(get_sequence_identity, args=(t1_secondary_structure,))
tm_dataframe["E-Loop Accuracy"] = tm_dataframe["Secondary Structure"].astype(str).str[44:75].apply(get_sequence_identity, args=(t1_secondary_structure[44:75],))

#Calculate maximum sequence identity
all_sequences = tm_dataframe["Sequence"].to_list()
max_identities = get_max_identity_parallel(all_sequences, "Sequences/all_encapsulins.fasta", 8)
identity_dataframe = pd.DataFrame.from_dict(max_identities, orient="index", columns=["Max Identity"])
identity_dataframe = identity_dataframe.reset_index()
identity_dataframe = identity_dataframe.rename(columns={"index": "Sequence", "Max Identity": "Max Identity"})
tm_dataframe = pd.merge(left=tm_dataframe, right=identity_dataframe, on="Sequence")

outfile = tm_dataframe.to_csv("ESMInverseFolding/data/masking+temperature_tests_tm_scores_all_data.csv")

sns.scatterplot(data=tm_dataframe, x="Max Identity", y="T1 (maritima)")
plt.show()