from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

all_records = []


for i, record in enumerate(SeqIO.parse("DiscriminatorModel/data/discriminator_generated.fasta", "fasta")):
    new_record = SeqRecord(record.seq, id=str(i+1), description="")

    all_records.append(new_record)

outfile = SeqIO.write(all_records, "DiscriminatorModel/data/generated_sequences_fixed_ids.fasta", "fasta")