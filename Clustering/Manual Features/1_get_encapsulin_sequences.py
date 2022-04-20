from os import access
import urllib.request
import pandas as pd
from Bio import SeqIO

base_url="https://www.uniprot.org/uniprot/"

encapsulin_data = pd.read_excel("Clustering Analysis/giesen_et_al_encapsulins.xlsx", sheet_name=[0,1,2,3], engine="openpyxl")

accession_data = []

for dataframe in encapsulin_data.values():
    accession_data.append(dataframe["Enc Uniprot Accession"].tolist())

# for accession in accession_data[0]:
#     link = base_url + accession + ".fasta"
#     filename = "family1/" + accession + ".fasta"
#     urllib.request.urlretrieve(link, filename)
    
# for accession in accession_data[1]:
#     link = base_url + str(accession) + ".fasta"
#     filename = "family2/" + str(accession) + ".fasta"
#     urllib.request.urlretrieve(link, filename)

for accession in accession_data[2]:
    link = base_url + str(accession) + ".fasta"
    filename = "family3/" + str(accession) + ".fasta"
    urllib.request.urlretrieve(link, filename)

for accession in accession_data[3]:
    link = base_url + str(accession) + ".fasta"
    filename = "family4/" + str(accession) + ".fasta"
    urllib.request.urlretrieve(link, filename)