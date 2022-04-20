#Imports
import pandas as pd
from Bio import SeqIO

#Import encapsulin Excel spreadsheet (Andreas and Giessen, 2021) as a dictionary of 4 dataframes
#Each dataframe represents one family of encapsulins
encapsulin_data = pd.read_excel("Clustering Analysis/giesen_et_al_encapsulins.xlsx", sheet_name=[0,1,2,3], engine="openpyxl")

#Add a consistent "Family" column to each dataframe
encapsulin_data[0]["Family"] = pd.Series([1 for x in range(len(encapsulin_data[0].index))])
encapsulin_data[1]["Family"] = encapsulin_data[1]["Family 2 (A or B?)"].str.slice(-2,)
encapsulin_data[2]["Family"] = pd.Series([3 for x in range(len(encapsulin_data[0].index))])
encapsulin_data[3] = encapsulin_data[3].rename(columns={'Family Type':"Family"})

#Rename the "Natural Product Cluster Type" column to "Cargo Type" to be consistent with all the other families
encapsulin_data[2] = encapsulin_data[2].rename(columns={'Natural Product Cluster Type:':"Cargo Type"})

#Drop unneeded columns from the dataframes
#Leave behind only Accession, Phylum, Kingdom, Cargo Type, and Family
encapsulin_data[0]= encapsulin_data[0].drop(columns=['Gene Name', 'Organism',
       'Taxonomy ID', 'Cargo Uniprot Accession',
       'Manually Curated Cargo Accession',
       'Genome Neighbors (Genome neighborhood index compared to Encapsulin (Negative values are upstream)| Uniprot Accession | Description| Pfam Family)'])
encapsulin_data[1] = encapsulin_data[1].drop(columns=['Gene Name', 'Organism',
       'Taxonomy ID', 'Cysteine Desulfurase Accessions', 'Terpene Cyclase Accessions',
       'Polyprenyl Transferase Accession', 'Xylulose Kinase Accession',
       'Split Sequence?','Family 2 (A or B?)', 'Double Encapsulin?', 'Double Enc Operon Position',
       'Enc_Neighbor_Accession',
       'Genome Neighbors (Genome neighborhood index compared to Encapsulin (Negative values are upstream)| Uniprot Accession | Description| Pfam Family)'])
encapsulin_data[2] = encapsulin_data[2].drop(columns=['Gene Name', 'Organism',
       'Taxonomy ID',
       'Genome Neighbors (Genome neighborhood index compared to Encapsulin (Negative values are upstream)| Uniprot Accession | Description| Pfam Family)'])
encapsulin_data[3] = encapsulin_data[3].drop(columns=['Gene Name', 'Organism',
       'Taxonomy ID',
       'Cargo Accession',
       'Genome Neighbors (Genome neighborhood index compared to Encapsulin (Negative values are upstream)| Uniprot Accession | Description| Pfam Family)'])

#Merge all 4 dataframes into one big dataframe
all_encapsulins = pd.concat([df for keys, df in encapsulin_data.items()], ignore_index=True)

#Drop rows with all missing values
all_encapsulins = all_encapsulins.dropna(axis=0, how="all")

#Replace missing phyla with "None"
all_encapsulins["Phylum"] = all_encapsulins["Phylum"].fillna("None")

#Convert all values in "Family" column to a string
all_encapsulins["Family"] = all_encapsulins["Family"].astype(str)

#Load sequence data from previously downloaded FASTA files
filenames = ["Sequences/all/" + accession + ".fasta" for accession in all_encapsulins["Enc Uniprot Accession"].tolist()]

sequences_dictionary = {}

for filename in filenames:
       try:
              seq_record = SeqIO.read(filename, "fasta")
              sequences_dictionary[filename.removeprefix("Sequences/all/")[:-6]] = str(seq_record.seq)
       except ValueError:
              continue

#Import sequences into dataframe
all_encapsulins["Sequence"] = all_encapsulins["Enc Uniprot Accession"].map(sequences_dictionary)

#Drop rows with missing sequence values
all_encapsulins = all_encapsulins.dropna(axis=0)

#Write dataframe to a file for faster loading in future
output = all_encapsulins.to_csv('Clustering Analysis/all_encapsulins.csv', index=False)