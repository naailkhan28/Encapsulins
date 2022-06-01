import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scoring_functions import extract_secondary_structure, get_path, get_sequence_identity
import numpy as np

#Load dataframes
tm_dataframe_1 = pd.read_csv("ESMInverseFolding/data/scores/af2_temperature_tests_TM_scores.csv")
tm_dataframe_1 = tm_dataframe_1.iloc[:, 1:]
#tm_dataframe_1 = pd.melt(tm_dataframe_1, id_vars=["id"], value_vars=["T1 (maritima)", "T1 (elongatus)", "T3", "T4"], value_name="TM Score", var_name="Target")

tm_dataframe_2 = pd.read_csv("ESMInverseFolding/data/scores/af2_temperature_tests_TM_scores_lowestID.csv")
tm_dataframe_2 = tm_dataframe_2.iloc[:, 1:]
#tm_dataframe_2 = pd.melt(tm_dataframe_2, id_vars=["id"], value_vars=["T1 (maritima)", "T1 (elongatus)", "T3", "T4"], value_name="TM Score", var_name="Target")

tm_dataframe_3 = pd.read_csv("ESMInverseFolding/data/scores/af2_masking_tests_TM_scores.csv")
tm_dataframe_3 = tm_dataframe_3.iloc[:, 1:]
#tm_dataframe_3 = pd.melt(tm_dataframe_3, id_vars=["id"], value_vars=["T1 (maritima)", "T1 (elongatus)", "T3", "T4"], value_name="TM Score", var_name="Target")

tm_dataframe_all = pd.concat([tm_dataframe_1, tm_dataframe_2, tm_dataframe_3])

sequences_dataframe = pd.read_csv("ESMInverseFolding/data/temperature_tests.csv")
sequences_dataframe = sequences_dataframe.iloc[:, 1:]

hmm_dataframe = pd.read_csv("ESMInverseFolding/data/HMMs/scores/family1_hmm_scores_dataframe.csv")
hmm_dataframe = hmm_dataframe.iloc[:, 1:]
sequences_dataframe = pd.merge(left=sequences_dataframe, right=hmm_dataframe, on="id")
sequences_dataframe = pd.merge(left=sequences_dataframe, right=tm_dataframe_all, on="id")

sequences_dataframe["HMM E-Value"] = -1 * np.log10(sequences_dataframe["HMM E-Value"])

sns.scatterplot(data=sequences_dataframe, x="HMM E-Value", y="T1 (maritima)", hue="Temperature")
plt.vlines(x=12, ymin=-0, ymax=1.0, linestyles="dashed", color="black")
plt.show()

# #Merge both dataframes into a single dataframe
# ids = tm_dataframe_all["id"].to_list()
# sequences_dataframe = sequences_dataframe.query("id in @ids")
# tm_dataframe = pd.merge(left=sequences_dataframe, right=tm_dataframe_all, on="id")

# #Fit linear regression to Log Likelihood Score vs TM-Score
# X = tm_dataframe["Log Likelihood Score"].to_numpy().reshape(-1, 1)
# y = tm_dataframe["T1 (maritima)"].to_numpy()
# regr = LinearRegression().fit(X, y)
# y_pred = regr.predict(X)

# #Add secondary structure sequences
# tm_dataframe["ss_paths"] = tm_dataframe["id"].apply(get_path)
# tm_dataframe["Secondary Structure"] = tm_dataframe["ss_paths"].apply(extract_secondary_structure)
# tm_dataframe = tm_dataframe.drop(["ss_paths"], axis=1)

# #Calculate secondary structure accuracy
# t1_secondary_structure = extract_secondary_structure("ESMInverseFolding/data/dssp_secondary_structures/7mu1.txt")
# tm_dataframe["Secondary Structure Accuracy"] = tm_dataframe["Secondary Structure"].apply(get_sequence_identity, args=(t1_secondary_structure,))
# tm_dataframe["E-Loop Accuracy"] = tm_dataframe["Secondary Structure"].astype(str).str[44:75].apply(get_sequence_identity, args=(t1_secondary_structure[44:75],))

# tm_dataframe["HMM E-Value"] = -1 * np.log10(tm_dataframe["HMM E-Value"])

# maritima_e_value = -1 * np.log10(2e-97)
# maritima_score = 311.7

# sns.scatterplot(data=tm_dataframe, x="HMM Score", y="T1 (maritima)")
# plt.vlines(x=maritima_score, ymin=0, ymax=1.0, linestyles="dashed", color="red")
# plt.xlabel("HMM Score")
# plt.xlim(0, 45)
# plt.show()