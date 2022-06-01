from Bio import SeqIO
import nwalign3 as nw
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#Calculate the sequence identity between two equal length sequences
def get_sequence_identity(sequence1, sequence2):

    num_matches =  sum([a == b for a, b, in zip(sequence1, sequence2)])

    length = len(sequence1)

    return num_matches / length

#Get the maximum sequence identity between a query sequence and a database of sequences
def get_max_identity(query_sequence, database):
    max_id = 0

    for record in SeqIO.parse(database, "fasta"):
        alignment = nw.global_align(s1=query_sequence, s2=str(record.seq), matrix="BLOSUM62")
        identity = get_sequence_identity(alignment[0], alignment[1])
        if identity > max_id:
            max_id = identity
        
    return(max_id)

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

#Calcualte the relative solvent exposed surface area for each residue in a DSSP file
def extract_rsa_values(infile):

    #Surface areas for each residue from Tien et al, 2013 (theoretical)
    surface_areas = {"A":129.0,"R":274.0,"N":195.0,"D":193.0,"C":167.0,
                    "E":223.0,"Q":225.0, "G":104.0,"H":224.0,"I":197.0,
                    "L":201.0,"K":236.0, "M":224.0,"F":240.0,"P":159.0,
                    "S":155.0, "T":172.0,"W":285.0,"Y":263.0,"V":174.0}
    
    relative_surface_areas = {}

    with open(infile, "r") as f:
        for _ in range(28):
            f.readline()
        
        for i, line in enumerate(f):
            residue = line[13]
            exposed_surface_area = float(line[35:38])
            total_surface_area = surface_areas[residue]

            relative_surface_areas[i+1] = exposed_surface_area / total_surface_area
    
    return(relative_surface_areas)

#Return hydrophobicity values for a residue
def get_hydrophobicity_scores(residue):

    #Hydrophobicity scale from Wimley and White, 1996
    hydrophobicities = {"E":-1.23,"D":-2.02,"N":-0.42,"Q":-0.58,"K":-0.99,
                        "R":-0.81,"H":-0.96,"G":-0.01,"P":-0.45,"S":-0.13,
                        "T":-0.14,"C":0.24,"M":0.23,"A":-0.17,"V":-0.07,
                        "I":0.31,"L":0.56,"F":1.13,"W":1.85,"Y":0.94}
    
    return(hydrophobicities[residue])

#Remove illegal characters from amino acid sequences
def remove_illegal_chars(sequence):
    sequence = sequence.replace("B", "N")
    sequence = sequence.replace("Z", "Q")
    sequence = sequence.replace("J", "I")
    sequence = sequence.replace("X", "I")
    sequence = sequence.replace("<eos>", "")
    sequence = sequence.replace("<mask>", "")
    sequence = sequence.replace("<null_1>", "")
    sequence = sequence.replace("<af2>", "")
    return sequence

#Get instability index of a protein sequence
def get_instability_index(sequence):
    analysis_object = ProteinAnalysis(remove_illegal_chars(sequence))

    return(analysis_object.instability_index())