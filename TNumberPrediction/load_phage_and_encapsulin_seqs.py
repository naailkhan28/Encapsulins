import pandas as pd
from Bio import SeqIO

#Load phage T-number data from XLSX file
encapsulin_data = pd.read_excel("TNumberPrediction/data/phage_tnumber_data.xlsx", engine="openpyxl")

#Drop all columns except accession, sequence, and T number
encapsulin_data = encapsulin_data[["NCBI_GENPEPT_PROTEIN_ID", "TRANSLATION", "Tnearest"]]
encapsulin_data = encapsulin_data.rename(columns={"NCBI_GENPEPT_PROTEIN_ID": "Accession", "TRANSLATION": "Sequence", "Tnearest": "T"})

#Dictionary of known encapsulin structure accessions with their T numbers
known_encapsulins = {
    "Q9WZP2": 1,
    "A0R4H0": 1,
    "K5BEG2": 1,
    "D0LZ74": 1,
    "A0A3G6X0H1": 1,
    "Q1D6H4": 3,
    "I6U7J4": 3,
    "A0A0F5HPP7": 4
}

#Remove sequences above maximum length - ESM model can't handle sequences more than 1024 residues
encapsulin_data["Length"] = encapsulin_data["Sequence"].apply(len)
encapsulin_data = encapsulin_data[encapsulin_data["Length"] < 1024]

#Add known encapsulin structures with T numbers to the DataFrame
for accession, T in known_encapsulins.items():
    record = SeqIO.read(f"Sequences/all/{accession}.fasta", "fasta")

    dict = {"Accession": accession, "Sequence": str(record.seq), "T": T}

    encapsulin_data = encapsulin_data.append(dict, ignore_index=True)

#Write dataframe to CSV file
outfile = encapsulin_data.to_csv("TNumberPrediction/all_seqs.csv")