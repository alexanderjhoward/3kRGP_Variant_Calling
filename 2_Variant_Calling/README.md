# Step 2: Call sample genome variants relative to reference genome

With our genomic locations determined, we can now go through the process of downloading sample genomes, saving genomic regions of interest, and calling variants of the sample genome relative to the Nipponbare genome.

This is all done through the following script.

```{bash}

	screen -S PSY_run # Setting a session name for the screen helps distinguish it

	while read i; do
        	sbatch Scripts/CallOsPSY.sh $i
		sleep 30s
	done < ../1_Reference_Genome/Output/OsPSY_locs.txt

	## press Ctrl+a d to detach screen and let it run in the background

	screen -ls # lists all active screens, should see one named (string of numbers).PSY_run 
	screen -r 2771512 # The string of numbers for my session, returns my screen

```
