# Mystery Sample Analysis

Now we are going to process and analyze the sequencing data we generated earlier. 

## Basecall Fast5s & Demultiplex samples

Oxford Nanopore sequencers output the raw data in fast5 format. Fast5s contain the electrical conductivity signal from each nanopore. Each set of bases that are fed through the pore create a characteristic change in electrical conductivity which can be deciphered into the sequence of bases. To convert from this raw electrical signal to bases, we used a basecaller program called Guppy which outputs the data as a fastq (this is the same format of the data we downloaded from GEO). 

```bash
#recall your folder name
export me=$GROUP_SCRATCH/biochem_minicourse_2023/<your_dir>

#Go to your data folder we created earlier & create subfolders that the fastq files will go in. Change sampleX to your sample number.

cd $me/data
mkdir sampleX
```
Copy the data for your sample from the straightlab foler. Change sampleX to your sameple number (ie. sample2). Repeat this for each sample number for your group
```
cp $GROUP_SCRATCH/biochem_minicourse_2023/straightlab/data/samples/sampleX.fastq $me/data/samples
```

## BasicQC
As we did earlier for our downloaded data, we'll take a look at the quality of our sequencing run.

```
cd $me/data/samples
```
Create histogram and save it into a file. 
Replace sampleX with the sample# for your group in both the input fastq file & the output readlength histogram file
```
cat $me/data/samples/sampleX.fastq | awk '(NR%4==2){print length($0)}' | sort | uniq -c | sort -k1,1nr > sampleX_readlength.hist.txt

less sampleX_readlength.hist.txt
```


##Using Nanostat
We used FastQC earlier to get some qc metrics, but we're now working with long read data, and thus we'll only use NanoStat for metrics.
```
NanoStat --fastq $me/data/samples/sampleX.fastq -n sampleX_nanostat.summary
```
Transfer the output file to the Downloads folder on your computer. Run this command on your local terminal.
```
rsync -ah --progress <username>@dtn.sherlock.stanford.edu:/scratch/groups/astraigh/biochem_minicourse_2023/<your_dir>/data/samples/sampleX_nanostat.summary ~/Downloads
```
Open the file using your text edit application.


## Blast

Blast takes a fasta as an input, so we need to convert the fastq into fasta.

```bash
#move to the folder for this dataset
cd $me/data/samples

awk '(NR%4==1){print ">"$0; getline; print $0}' sampleX.fastq > sampleX.fasta
```

Let's look at the first 50 reads first
```
head -n 100 sampleX.fasta > sampleX.head.fasta
```

We have multipe databases installed on sherlock, which correspond to some of those available on ncbi blast online. Here, because we know we what type of samples we are dealing with we have downloaded just a few databases and combined them into one blastdb. The path of this database is  
```
/scratch/groups/astraigh/biochem_minicourse_2023/straightlab/blastdb/combined_db

#We can set a variable for this blast database so we don't have to type out the full path eachtime
export myblastdb=/scratch/groups/astraigh/biochem_minicourse_2023/straightlab/blastdb/combined_db

```

Let's run blast using this database. We will put the blast output in a `blast` folder

```bash
#create folder for output
cd $me/data/samples
mkdir blast

#run blast
blastn -query <(head -n 100 sampleX.fasta) -db $myblastdb -num_threads 2 > blast/sampleX.rawoutput.first50.txt
```

Taking a look at this file with `less` shows that this command returns the blast hits in a format similar to that of we ontain when running blast online. 
```
less blast/sampleX.rawoutput.first50.txt
```


We can instead ask blast to produce a tabulated output, which is more computationally friendy, using the `--outfmt 6` flag.

```
blastn -query <(head -n 100 sampleX.fasta) -db $myblastdb -num_threads 2 -subject_besthit -culling_limit 2 -outfmt "6 stitle score mismatches" > blast/sampleX_tabulated.output.first50.txt
```
Let's sort the output and print the uniq entries
```
cat blast/sampleX_tabulated.output.first50.txt | awk '{OFS=" "; print $1,$2,"\t",$(NF)}' | sort -k3,3 -nr > blast/sampleX_sort.output.first50.tsv
```
Based on the blast output, which organism do you think each DNA sample came from?


There are several other output format options that you can use to customize blast's output. More info here: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/

We can run a simple python script to summarize the blast output by species

```
python /scratch/groups/astraigh/biochem_minicourse_2023/straightlab/scripts/summarize_species.py blast/sampleX_sort.output.first50.tsv blast/sampleX_species_summary.tsv
```
