############### Third step

# Steps required to run ARACNe

# 1. Calculate a threshold for Mutual Information
# 2. Run ARACNe on bootstraps of the input matrix
# 3. Consolidate, i.e. combine the bootstraps into a final network file

echo -e "Step2: consolidation of candidate microRNAs by network analysis \n"
echo -e "Network reconstruction using ARACNE \n"


for d in ../results/boot/*;
	do
	
	eval "rna=`ls $d/*expr*`"
	eval "mirna=`ls $d/*mirna*`"
	

	#suffix1=_expr_new.txt
	#suffix2=_mirna_new.txt
	#eval "rna=$d$suffix1"
	#eval "mirna=$d$suffix2"

# 1. Calculate threshold
	
	java -Xmx5G -jar ./ARACNe-AP/dist/aracne.jar -e $rna  -o $d --tfs $mirna --pvalue 1E-8 --seed 1 --calculateThreshold

# 2. Run ARACNe on bootstraps of the input matrix

	for i in {1..100}
	do
		java -Xmx5G -jar ./ARACNe-AP/dist/aracne.jar -e $rna  -o $d --tfs $mirna --pvalue 1E-8 --seed $i 
	done
	
# 3. Consolidate, i.e. combine the bootstraps into a final network file

	java -Xmx5G -jar ./ARACNe-AP/dist/aracne.jar -o $d --consolidate

done; 

