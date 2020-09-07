# Taxonomic analysis

We start with a fastq file. We want to identify where the reads are coming from. We will start with some simple blast. 


## Blast

Blast takes a fasta as an input, so we need to convert the fastq into fasta.

```
awk '(NR%4==1){print ">"$0; getline; print $0}' raw/SRR8351023.fastq >raw/SRR8351023.fasta
```

Let's look at the first 10 reads first
```
head -n 20 raw/SRR8351023.fasta >raw/SRR8351023.head.fasta
```

We have multipe databases installed on sherlock, which correspond more or less to the ones available on ncbi blast online. Here, because we know we are dealing with a bacterial sample, we will use the database containing "representative" refseq prokaryote genomes. The path of this database is  
```
/scratch/groups/astraigh/kraken/ref_prok_rep_genomes/blast/ref_prok_rep_genomes
```

Let's run blast using this database. We will put the blast output in a `blast` folder

```
#create folder for output
mkdir -p blast

#this is just a technical thing for blast to run well
export BLASTDB=/scratch/groups/astraigh/kraken

#run blast
blastn -query raw/SRR8351023.head.fasta -db /scratch/groups/astraigh/kraken/ref_prok_rep_genomes/blast/ref_prok_rep_genomes >blast/rawoutput.first10.txt
```


Note we can use "process substitution" with `<()` to get the first 10 reads on the fly without having to create an intermediate file

```
blastn -query <(head -n 20 raw/SRR8351023.fasta) -db /scratch/groups/astraigh/kraken/ref_prok_rep_genomes/blast/ref_prok_rep_genomes >blast/raw.output.first10.txt
```

Taking a look at this file with `less` shows that this command returns the blast hits in a format similar to that of we ontain when running blast online. 

We can instead ask blast to produce a tabulated output, which is more computationally friendy, using the `--outfmt 6` flag.

```
blastn -query <(head -n 40 raw/SRR8351023.fasta) -db /scratch/groups/astraigh/kraken/ref_prok_rep_genomes/blast/ref_prok_rep_genomes -outfmt "6 qseqid sseqid pident length mismatches gapopen evalue bitscore ssciname" > >blast/tabulated.output.first10.txt
```

As you can see, blast takes a very long time to run, and there is an additionnal complexity which is that it returns more than one hit per read.

## Running Kraken
Rather than running blasts, taxonomic classification can be done using dedicated tools. Kraken2 (https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) is a commonly used taxonomic sequence classifier, which examines the k-mers within a query sequence and compare these k-mers with those present in a taxonomic database. 

We have prebuild a taxonomic kmer database for kranken2 on Sherlock. This database contains RefSeq complete genomes for all archea, bacteria, virus, and the human genome. 


```
kraken2 -db /scratch/groups/astraigh/kraken/standard --threads 2 --report kraken/standard.kreport --use-names --output kranken/standard.out.txt raw/SRR8351023.fastq
```

Take a look at the kraken report

```
less kraken/standard.kreport
```

Now we want to extract species level abundance and sort by the most abundant first


```
cat kraken/standard.kreport | sed 's/ \+ /\t/g' | sed 's/ *$//' | tr -s "\t" | cut -c 2- | awk -F $'\t' '($4=="S")' | sort -k1,1nr >kraken/standard.species.txt
```

## Plotting abundance

Now we have some data ready to load in R for plotting. Follow the instructions in notebook 04 to open a remote R session. We will do some analysis within the R notebook.


