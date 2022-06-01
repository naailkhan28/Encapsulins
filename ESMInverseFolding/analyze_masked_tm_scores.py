import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
from scoring_functions import extract_secondary_structure, get_path, get_sequence_identity
from sklearn.linear_model import LinearRegression

#Load dataframes
tm_dataframe = pd.read_csv("ESMInverseFolding/data/scores/af2_masking_tests_TM_scores.csv")
tm_dataframe = tm_dataframe.iloc[:, 1:]
#tm_dataframe = pd.melt(tm_dataframe, id_vars=["id"], value_vars=["T1 (maritima)", "T1 (elongatus)", "T3", "T4"], value_name="TM Score", var_name="Target")

sequences_dataframe = pd.read_csv("ESMInverseFolding/data/masked_eloop_sequences.csv")
sequences_dataframe = sequences_dataframe.iloc[:, 1:]

#Merge both dataframes into a single dataframe
ids = tm_dataframe["id"].to_list()
sequences_dataframe = sequences_dataframe.query("id in @ids")
tm_dataframe = pd.merge(left=sequences_dataframe, right=tm_dataframe, on="id")

#Add secondary structure sequences
tm_dataframe["ss_paths"] = tm_dataframe["id"].apply(get_path)
tm_dataframe["Secondary Structure"] = tm_dataframe["ss_paths"].apply(extract_secondary_structure)
tm_dataframe = tm_dataframe.drop(["ss_paths"], axis=1)

#Calculate secondary structure accuracy
t1_secondary_structure = extract_secondary_structure("ESMInverseFolding/data/dssp_secondary_structures/7mu1.txt")
tm_dataframe["Secondary Structure Accuracy"] = tm_dataframe["Secondary Structure"].apply(get_sequence_identity, args=(t1_secondary_structure,))
tm_dataframe["E-Loop Accuracy"] = tm_dataframe["Secondary Structure"].astype(str).str[44:75].apply(get_sequence_identity, args=(t1_secondary_structure[44:75],))

print(tm_dataframe)