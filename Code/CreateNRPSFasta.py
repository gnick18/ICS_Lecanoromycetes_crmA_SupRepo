from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os 
import pdb
import shutil
import copy
from Bio.SeqFeature import FeatureLocation


# Define column headers as per HMMER documentation
columns = [
"target_name", "target_accession", "tlen", "query_name", "query_accession", "qlen",
"E_value", "score", "bias", "domain_num", "domain_total", "c_Evalue", "i_Evalue",
"domain_score", "domain_bias", "hmm_from", "hmm_to", "ali_from", "ali_to",
"env_from", "env_to", "acc", "description_of_target"
]
# Load data into a DataFrame
crm_hmmer_search = pd.read_csv('/Volumes/T7/MiscProjects/Singh/April_FinalWork/crmTree/hmmsearch.out', comment="#",delim_whitespace=True, header=None, names=columns)

crmTable = pd.read_csv("/Volumes/T7/MiscProjects/Singh/April_FinalWork/GetNRPS_crm/NRPS_CRM_FinalTable.tsv", sep='\t')

nrpsSeqs = []
proteinsAdded = []

fused = crmTable[crmTable['NRPS_attached_ICS?'] == "Yes"]
split = crmTable[crmTable['NRPS_attached_ICS?'] == "No"]

## For split we just grab the nrps seq and make a seqrecord
for i, row in split.iterrows():
	nrpsProtein = row['nrps_protein']
	if nrpsProtein in proteinsAdded:
		continue
	nrpsSeqs.append(SeqRecord(row['nrps_proteinSeq'], id=row['nrps_protein'] + "_split", description=""))
	proteinsAdded.append(nrpsProtein)

## For the fused we want to remove the ics domain, and add the rest.
for i, row in fused.iterrows():
	icsSeq = row['crm_proteinSeq']
	fullcrmLength = len(icsSeq)
	icsName = row['crm_protein']
	if icsName in proteinsAdded:
		continue
	accession = row['Accession']
	species = row['Species']
	## make a seqrecord of the ics seq
	icsSeqRecord = SeqRecord(icsSeq, id=species+"_"+icsName, description=accession)
	icsHmmerRow = crm_hmmer_search[(crm_hmmer_search['target_name'].str.contains(icsName)) & (crm_hmmer_search['query_name'] == 'DIT1_PvcA')]
	if len(icsHmmerRow) == 0:
		print("ICS protein not found in hmmsearch output")
		continue
	if len(icsHmmerRow) > 1:
		## grab the first one
		icsHmmerRow = icsHmmerRow.iloc[0]
	seq_id = icsHmmerRow['target_name']
	start = int(icsHmmerRow['ali_from']) - 1  # convert to 0-based indexing
	end = int(icsHmmerRow['ali_to'])
	### check if start is closer to 0 than the end is to the length of the seq
	if start < (fullcrmLength - end):
		beginingORend = "begining"
	else:
		beginingORend = "end"

	nrpsSeqRecord = copy.deepcopy(icsSeqRecord)
	if beginingORend == "begining":
		nrpsSeqRecord.seq = icsSeq[end:]	
	elif beginingORend == "end":
		nrpsSeqRecord.seq = icsSeq[:start]
	#
	nrpsSeqRecord.id = species+"_"+icsName+"_nrps"
	nrpsSeqRecord.description = accession
	nrpsSeqs.append(nrpsSeqRecord)
	proteinsAdded.append(icsName)


## Write the new sequences to a FASTA file
outputFolder = "/Volumes/T7/MiscProjects/Singh/April_FinalWork/nrpsTree"
outputName = os.path.join(outputFolder, 'nrps_crm_domains.faa')
#
SeqIO.write(nrpsSeqs, outputName, 'fasta')
