import pandas as pd
import re

#Load TM scores from file
scores_path = "ESMInverseFolding/data/scores/t4_hmm_test_tm_scores_hi.txt"
tm_scores = {}

#Regex to match the 4x TM Scores written to the text file - four match groups corresponding to T4, T1 (elongatus), T1 (maritima), and T3 structures
id_pattern = re.compile(r'models\/(.*)_ranked')
scores_pattern = re.compile(r'6NJ8\.pdbTM-score= (0\.\d+).*6X8T\.pdbTM-score= (0\.\d+).*7MU1\.pdbTM-score= (0\.\d+).*7S2T\.pdbTM-score= (0\.\d+)')

#Load all text data into an object
with open(scores_path, "r") as f:
    all_lines = f.readlines()

#Split the text into chunks of 13 lines - each of these chunks is the data for one AlphaFold2 model
chunks = [all_lines[i:i+13] for i in range(0, len(all_lines)-13, 13)]
chunks = ["".join(line.rstrip() for line in chunk) for chunk in chunks]

#For each chunk of 13 lines, get the ID of the protein and the 4 TM scores against each template model
for i, chunk in enumerate(chunks):
    id = re.search(id_pattern, chunk).groups()[0]

    scores = [float(group) for group in re.search(scores_pattern, chunk).groups()]

    scores.insert(0, id)

    tm_scores[i] = scores

#Load this into a dataframe and save

tm_dataframe = pd.DataFrame.from_dict(tm_scores, orient="index", columns=["id", "T4", "T1 (elongatus)", "T1 (maritima)", "T3"])
outfile = tm_dataframe.to_csv("ESMInverseFolding/data/scores/hmm_tests_tm_scores_hi.csv")