from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

cleaned_records = []

for old_record in SeqIO.parse("ESM Sequence Generation/generated_encapsulins_1000_p_furiosus_family4.fasta", "fasta"):

    id = old_record.id

    old_seq = str(old_record.seq)

    new_seq = Seq(old_seq[4:-4])

    new_record = SeqRecord(new_seq, id=id, description="")

    cleaned_records.append(new_record)

outfile = SeqIO.write(cleaned_records, "ESM Sequence Generation/processed_generated_1000sequences_pfuriosus_family4_seed.fasta", "fasta")