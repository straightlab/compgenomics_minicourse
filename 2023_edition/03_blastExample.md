# Taxonomic analysis

We start with a fastq file. We want to identify where the reads are coming from. We will start with some simple blast. 


## Blast

Blast takes a fasta as an input, so we need to convert the fastq into fasta.

```bash
#recall your folder name
export me=$GROUP_SCRATCH/biochem_minicourse_2023/<your_dir>

#move to the folder for the dataset
cd $me/data/external

awk '(NR%4==1){print ">"$0; getline; print $0}' SRR13403380_subset.fastq > SRR13403380_subset.fasta
```

Let's look at the first 10 reads first
```
head -n 20 SRR13403380_subset.fasta > SRR13403380.head.fasta
```

We have multipe databases installed on sherlock, which correspond to some of those available on ncbi blast online. Here, because we know what type of samples we are dealing with we have downloaded just a few databases and combined them into one blastdb. The path of this database is  
```
/scratch/groups/astraigh/biochem_minicourse_2023/straightlab/blastdb/combined_db
```

Let's run blast using this database. We will put the blast output in a `blast` folder

```bash
#create folder for output
mkdir -p blast

#we can set a variable for this blast database so we don't have to type out the full path everytime
myblastdb=/scratch/groups/astraigh/biochem_minicourse_2023/straightlab/blastdb/combined_db

#run blast
blastn -query SRR13403380.head.fasta -db $myblastdb -num_threads 2 > blast/SRR13403380_rawoutput.first10.txt
```

Taking a look at this file with `less` shows that this command returns the blast hits in a format similar to what we obtain when running blast online. 

Note we can use "process substitution" with `<()` to get the first 10 reads on the fly without having to create an intermediate file

```
blastn -query <(head -n 20 SRR13403380_subset.fasta) -db $myblastdb -num_threads 2 > blast/SRR13403380_raw.output.first10.txt
```

Taking a look at this file with `less` shows that this command returns the blast hits in a format similar to the format we obtain when running blast online. 

We can instead ask blast to produce a tabulated output, which is more computationally friendy, using the `--outfmt 6` flag.

```
blastn -query <(head -n 20 SRR13403380_subset.fasta) -db $myblastdb -num_threads 2 -subject_besthit -outfmt "6 stitle score" > blast/SRR13403380_tabulated.output.first10.txt
```

