import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

#Load sequence and scores dataframes
masked_sequences_dataframe = pd.read_csv("ESMInverseFolding/data/masked_eloop_sequences.csv")
masked_sequences_dataframe = masked_sequences_dataframe.iloc[:, 1:]
temperature_sequences_dataframe = pd.read_csv("ESMInverseFolding/data/temperature_tests.csv")
temperature_sequences_dataframe = temperature_sequences_dataframe.iloc[:, 1:]

#Load TM score dataframes
masked_tm_dataframe = pd.read_csv("ESMInverseFolding/data/scores/af2_masking_tests_TM_scores.csv")
masked_tm_dataframe = masked_tm_dataframe.iloc[:, 1:]
temperature_tm_dataframe = pd.read_csv("ESMInverseFolding/data/scores/af2_temperature_tests_TM_scores.csv")
temperature_tm_dataframe = temperature_tm_dataframe.iloc[:, 1:]
temperature_lowID_tm_dataframe = pd.read_csv("ESMInverseFolding/data/scores/af2_temperature_tests_TM_scores_lowestID.csv")
temperature_lowID_tm_dataframe = temperature_lowID_tm_dataframe.iloc[:, 1:]

#Load discriminator scores dataframe
discriminator_scores_dataframe = pd.read_csv("ESMInverseFolding/data/scores/discriminator_scores.csv", header=0, names=["id", "Discriminator Score"])

#Merge all TM score and sequence frames
all_tm_dataframes = pd.concat([masked_tm_dataframe, temperature_tm_dataframe, temperature_lowID_tm_dataframe])
all_sequences_dataframe = pd.concat([masked_sequences_dataframe, temperature_sequences_dataframe])

#Merge discrminator scores dataframe with sequences dataframe
all_generated_ids = all_sequences_dataframe["id"].to_list()
discriminator_scores_dataframe = discriminator_scores_dataframe.query("id in @all_generated_ids")
all_scores_dataframe = pd.merge(left=all_sequences_dataframe, right=discriminator_scores_dataframe, on="id")

sns.scatterplot(data=all_scores_dataframe, x="Log Likelihood Score", y="Discriminator Score")
plt.show()