from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sequence_list = ['148_maritima_0.2', '829_maritima_0.2', '979_maritima_0.2', '966_maritima_0.2', '110_maritima_0.2', '465_maritima_0.2', '482_maritima_0.2', '111_maritima_0.2', '641_maritima_0.2', '368_maritima_0.2', '464_maritima_0.2', '992_maritima_0.2', '159_maritima_0.2', '599_maritima_0.2', '904_maritima_0.2', '205_elongatus_0.3', '773_elongatus_0.2', '662_thermotolerans_0.3', '10_elongatus_0.3', '426_elongatus_0.3', '812_elongatus_0.3', '50_thermotolerans_0.3', '231_thermotolerans_0.3', '679_elongatus_0.3', '90_thermotolerans_0.3', '629_thermotolerans_0.2', '193_elongatus_0.3', '668_elongatus_0.3', '922_elongatus_0.3', '791_maritima_0.3', '439_thermotolerans_0.2', '965_thermotolerans_0.2', '803_thermotolerans_0.2', '718_elongatus_0.2', '609_thermotolerans_0.2', '376_elongatus_0.2', '153_thermotolerans_0.2', '639_thermotolerans_0.2', '932_thermotolerans_0.2', '831_thermotolerans_0.2', '939_thermotolerans_0.2', '185_thermotolerans_0.2', '303_thermotolerans_0.2', '532_xanthus_0.2', '632_thermotolerans_0.2', '278_maritima_0.3', '760_xanthus_0.3', '658_thermotolerans_0.3', '99_maritima_0.3', '7_maritima_0.3', '314_maritima_0.3', '286_elongatus_0.3', '187_elongatus_0.3', '646_maritima_0.3', '894_maritima_0.3', '631_elongatus_0.3', '755_maritima_0.3', '270_maritima_0.3', '161_maritima_0.3', '156_xanthus_0.3', '972_elongatus_0.3', '900_elongatus_0.2']

for record in SeqIO.parse("Sequences/Generated/masking_level_testing.fasta", "fasta"):
    if str(record.id) in sequence_list:
        new_id = str(record.id).replace(".", "_")

        new_record = SeqRecord(record.seq, id=new_id, description="")

        outfile = SeqIO.write(new_record, f"{new_id}.fasta", "fasta")
