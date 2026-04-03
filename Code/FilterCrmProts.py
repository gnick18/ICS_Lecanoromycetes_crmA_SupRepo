from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os 
import pdb

# Load full sequence dictionary from fasta
output_crmTreeFiles = '/Volumes/T7/MiscProjects/Singh/April_FinalWork/crmTree'
crmFastaPath = os.path.join(output_crmTreeFiles, 'crm_outgroups.faa')
crm = SeqIO.to_dict(SeqIO.parse(crmFastaPath, 'fasta'))

# Parse tblout file and collect domain ranges for DIT1_PvcA or PF05141*
domain_hits = {}
hmmPath = os.path.join(output_crmTreeFiles, 'hmmsearch.out')

# Define column headers as per HMMER documentation
columns = [
    "target_name", "target_accession", "tlen", "query_name", "query_accession", "qlen",
    "E_value", "score", "bias", "domain_num", "domain_total", "c_Evalue", "i_Evalue",
    "domain_score", "domain_bias", "hmm_from", "hmm_to", "ali_from", "ali_to",
    "env_from", "env_to", "acc", "description_of_target"
]

# Load data into a DataFrame
df = pd.read_csv(hmmPath, delim_whitespace=True, header=None, names=columns, comment='#')

# Filter for rows matching PF05141* or DIT1_PvcA
ditDomains = df[df['query_name'] == 'DIT1_PvcA']
# if there are duplicates, keep the row with the highest score
ditDomains = ditDomains.loc[ditDomains.groupby('target_name')['domain_score'].idxmax()]

# Extract domain regions and build new SeqRecord list
ditDomains_SeqIO = []
for index, row in ditDomains.iterrows():
    seq_id = row['target_name']
    start = int(row['ali_from']) - 1  # convert to 0-based indexing
    end = int(row['ali_to'])
    if seq_id in crm:
        original_record = crm[seq_id]
        domain_seq = original_record.seq[start:end]
        domain_record = SeqRecord(domain_seq, id=original_record.id + "_ICS", description="DIT1_PvcA domain")
        ditDomains_SeqIO.append(domain_record)
    else:
        print(f"Warning: {seq_id} not found in the CRM dictionary.")
        pdb.set_trace()

## Write the new sequences to a FASTA file
outputName = os.path.join(output_crmTreeFiles, 'DIT1_PvcA_crm_domains.faa')
SeqIO.write(ditDomains_SeqIO, outputName, 'fasta')

