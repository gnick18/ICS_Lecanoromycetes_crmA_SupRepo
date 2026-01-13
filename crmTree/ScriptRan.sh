###### Trimming the alignments with different trimal parameters
echo "Trimming the alignments with different trimal parameters"
trimal -in ./DIT1_PvcA_crm_domains_aln.faa -out ./DIT1_PvcA_crm_domains_gappyout.fas -gappyout
trimal -in ./DIT1_PvcA_crm_domains_aln.faa -out ./DIT1_PvcA_crm_domains_gt0-1.fas -gt 0.1
trimal -in ./DIT1_PvcA_crm_domains_aln.faa -out ./DIT1_PvcA_crm_domains_gt0-9.fas -gt 0.9
echo "Trimming complete."
###### IQ-Tree run
# ~Model selected on previous run with ModelFinder 
# Gappyout
echo "Running IQ-Tree on the trimmed alignments"
if [ ! -d "Gappyout" ]; then
    mkdir Gappyout
	# run the IQ-Tree
	echo "Running IQ-Tree on the Gappyout trimmed alignment"
	iqtree2 -s ./DIT1_PvcA_crm_domains_gappyout.fas -m Q.plant+I+G4 -T 8 --ufboot 1000 -safe 
	mv ./DIT1_PvcA_crm_domains_gappyout* ./Gappyout/
	echo "IQ-Tree run complete for Gappyout."
fi
# GT 0.1
if [ ! -d "GT0-1" ]; then
    mkdir GT0-1
	# run the IQ-Tree
	echo "Running IQ-Tree on the GT0-1 trimmed alignment"
	iqtree2 -s ./DIT1_PvcA_crm_domains_gt0-1.fas -m Q.plant+I+G4 -T 8 --ufboot 1000 -safe 
	mv ./DIT1_PvcA_crm_domains_gt0-1* ./GT0-1/
	echo "IQ-Tree run complete for GT0-1."
fi
# GT 0.9
if [ ! -d "GT0-9" ]; then
    mkdir GT0-9
	# run the IQ-Tree
	echo "Running IQ-Tree on the GT0-9 trimmed alignment"
	iqtree2 -s ./DIT1_PvcA_crm_domains_gt0-9.fas -m Q.plant+I+G4 -T 8 --ufboot 1000 -safe 
	mv ./DIT1_PvcA_crm_domains_gt0-9* ./GT0-9/
	echo "IQ-Tree run complete for GT0-9."
fi

echo "IQ-Tree runs complete."