from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO
from model_prediction import get_log_likelihoods
from sequence_generation_tools import tokenize, remove_illegal_chars

#Given a query sequence, find the sequence with the highest similarity in a FASTA file
#Similarity is calculated using pairwise alignment with BLOSUM62 scoring, gap open and extension penalties can be defined too
#Returns a tuple of the highest score and the accession of the sequence with that similarity score to the query
def get_max_sequence_similarity(query_sequence, database_file, gap_open_penalty=-0.5, gap_extend_penalty=-0.1, end_gap_score=0):

    query_sequence = remove_illegal_chars(query_sequence)

    #Initilize aligner object with parameters
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "global"
    aligner.open_gap_score = gap_open_penalty
    aligner.extend_gap_score = gap_extend_penalty
    aligner.target_end_gap_score = end_gap_score
    aligner.query_end_gap_score = end_gap_score

    scores = {}

    #Iterate through the database file and calculate similarity scores for each sequence
    for record in SeqIO.parse(database_file, "fasta"):
        target_sequence = remove_illegal_chars(str(record.seq))

        if "<eos>" in query_sequence:
            continue

        score = aligner.score(query_sequence, target_sequence)

        normalization_length = max(len(query_sequence), len(target_sequence))

        scores[score / normalization_length] = record.id

    #Return the maximum score and the sequence that produces it
    max_score = max(scores.keys())

    return(max_score, scores[max_score])

#Calculate average per-token log likelihood for a generated sequence
def get_average_log_likelihood(query_sequence, model, batch_converter, device):

    sequence_tokens = tokenize(query_sequence, batch_converter)
    sequence_tokens = sequence_tokens.to(device)
    all_likelihoods = get_log_likelihoods(sequence_tokens, model)

    total_log_likelihood = 0.0

    for position, likelihoods in enumerate(all_likelihoods[1:-1]):
        actual_residue = sequence_tokens[0][position + 1]

        total_log_likelihood += likelihoods[actual_residue]

    return((total_log_likelihood / len(query_sequence)).item())