import os
from matplotlib import pyplot as plt

mainpath = "Alphafold2Predictions/ESMInverseFolding/maskingtesting/"

folders = [folder for folder in os.listdir(mainpath)]

#Get all filenames in directory, and make a list of all filenames ending in ".pdb"

pdb_filenames = []

for directory in folders:
    if directory == ".DS_Store":
        continue
    all_files = os.listdir(mainpath + directory)
    for file in all_files:
        if file == "ranked_0.pdb":
            pdb_filenames.append(f"{mainpath}{directory}/{file}")


#Residue numbers and pLDDT scores are stored in a dictionary
#FILENAME_residues is a key whose value is a list of every residue number
#FILENAME_pLDDT is a key whose value is a list of every pLDDT score from the pdb file
all_values = {}
average_pLDDT = {}

for model in pdb_filenames:
    all_values["{}_residues".format(model)] = []
    all_values["{}_pLDDT".format(model)] = []

    with open(model, "r") as file:
        header = file.readline()

        all_lines = file.readlines()

    columns = []
    
    for line in all_lines:
        columns = line.split()    
        if columns[0] == "TER":
            break
        all_values["{}_residues".format(model)].append(float(columns[5]))
        all_values["{}_pLDDT".format(model)].append(float(columns[10]))

    average_pLDDT[model] = sum(all_values["{}_pLDDT".format(model)]) / len(all_values["{}_pLDDT".format(model)])   


#Set up figure in matplotlib
fig, axs = plt.subplots(1, 2)

#For every pdb file in this folder, plot its pLDDT scores against residue numbers
plots = [axs[0].plot(all_values["{}_residues".format(filename)], all_values["{}_pLDDT".format(filename)], label=filename[-16:-19]) for filename in pdb_filenames]
#plt.xlim(92, 170)
axs[0].set_xlabel("Residue")
axs[0].set_ylabel("pLDDT")
#plt.legend(loc="upper right")

#Plot average pLDDT scores for every pdb file
axs[1].bar([key[0:-4] for key in average_pLDDT], [average_pLDDT[key] for key in average_pLDDT], color=plt.rcParams['axes.prop_cycle'].by_key()['color'], edgecolor="black")
axs[1].set_xticklabels([filename[-16:-19] for filename in pdb_filenames])
axs[1].set_ylabel("Average pLDDT")
axs[1].set_ylim(50, 100)

plt.show()