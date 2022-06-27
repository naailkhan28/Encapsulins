from Bio import AlignIO
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import nwalign3 as nw
import time

#Calculate the sequence identity between two equal length sequences
def get_sequence_identity(sequence1, sequence2):

    num_matches =  sum([a == b for a, b, in zip(sequence1, sequence2)])

    length = len(sequence1)

    return num_matches / length

#Get the maximum sequence identity between a query sequence and a database of sequences
def get_maxid(query_sequence, database):
    max_id = 0

    for record in SeqIO.parse(database, "fasta"):
        alignment = nw.global_align(s1=query_sequence, s2=str(record.seq), matrix="BLOSUM62")
        identity = get_sequence_identity(alignment[0], alignment[1])
        if identity > max_id:
            max_id = identity
        
    return(max_id)

#Get total counts of all 21 amino acids in each column of an MSA
def get_residue_counts(alignmentfile, format):
    all_residues = defaultdict(list)
    alignment = AlignIO.read(alignmentfile, format)

    for record in alignment:
        for i, residue in enumerate(str(record.seq)):
            all_residues[i+1].append(residue)

    return(all_residues)

#Calculate Shannon Entropy for each column in an MSA
def calculate_variability(residue_counts, num_sequences):

    num_sequences = 1000
    alphabet = list("ACDEFGHIKLMNPQRSTVWY-")
    variability = {}

    for position, residues in residue_counts.items():

        frequencies = [residues.count(amino_acid) / num_sequences for amino_acid in alphabet]

        entropy = -1 * sum([(n * np.log(n)) if n else 0 for n in frequencies])

        variability[position] = entropy

    return(variability)

#Calculate the difference between per-column Shannon Entropies for two MSAs
def calculate_variability_difference(list1, list2):

    min_length = min(len(list1), len(list2))

    array1 = np.array(list1[:min_length])
    array2 = np.array(list2[:min_length])

    difference = abs(array1 - array2)

    positions = np.array([x for x in range(1, min_length + 1)])

    return(positions, difference)

#Extract a secondary structure sequence from a DSSP file
def extract_secondary_structure(infile):
    
    ss= []

    with open(infile, "r") as f:
        for _ in range(28):
            f.readline()
        
        for line in f:
            ss_element = line[16]
            if ss_element == " ":
                ss.append("L")
            else:
                ss.append(ss_element)

    return("".join(ss))

t4_generated = "TNWQEPTQDFNYNDLYSNICQGYIDLDTAAGFIMSKEALGAFIDHVNAVAVQNIREVRESNNDLARVFNLSEEQDSYSIHESITAWVLSLRSNDLHSQMNRPTPIQNAKSSASELQKVQMQDMIGGQSLTADNGTIALSGVLLWEEGKPDVTEEFLPEYTQALHSLASHNHQGNYVMWLSPPWCAHYKRSRIDSNVPVVAAAKDIINGGSIAVRKIYQATGIMTLHHIGLRQLKEMLKLVVICLQERGLDKPQWIINTIIVHHEFYGTHVHCQFNDH" 
path = "data/all_encapsulins.fasta"

start = time.time()
output = get_maxid(t4_generated, path)
end = time.time()
print("Took %f ms" % ((end - start) * 1000.0))