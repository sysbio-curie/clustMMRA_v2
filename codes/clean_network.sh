for d in ../results/boot/*;
	do
	eval "net=`ls $d/network*`"
	awk -F"\t" '{if ($2 !~ /miR/) {print }}' $net > cleaned_network.txt
	mv cleaned_network.txt $d
	
done;
