# Mystery Sample Analysis

Now we are going to process and analyze the sequencing data we generated earlier. 

## Basecall Fast5s & Demultiplex samples

Oxford Nanopore sequencers output the raw data in fast5 format. Fast5s contain the electrical conductivity signal from each nanopore. Each set of bases that are fed through the pore create a characteristic change in electrical conductivity which can be deciphered into the sequence of bases. To convert from this raw electrical signal to bases, we will use a basecaller program called Guppy which outputs the data as a fastq (this is the same format of the data we downloaded from GEO).

```bash
#recall your folder name
export me="teamStraight"

#Go to your data folder we created earlier & create subfolders that the basecalled data will go in. Change sampleX to your sample number.

cd $GROUP_SCRATCH/biochem_minicourse_2021/$me/data
mkdir sampleX
mkdir sampleX/basecalled

#Now we want to define a variable that points to your sample's fast5s. 
#We have uploaded the data from each sample to the following directory. Swap sampleX for your sample number.

fast5s=$GROUP_SCRATCH/biochem_minicourse_2021/shared/class_dat/sampleX
```
**Basecall**

Guppy takes an input_path, save_path, and config as inputs. Input_path is the directory where your fast5 files are located (in shared/class_dat/sampleX). Save_path is the directory where guppy will output the resulting fastq files. Config refers to the specific settings we input into guppy. We don't need to change any of these settings to we'll leave that as is.

```
/home/groups/astraigh/software/ont-guppy-cpu/bin/guppy_basecaller --input_path $fast5s --save_path "$GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/basecalled" --config /home/groups/astraigh/software/ont-guppy-cpu/data/dna_r9.4.1_450bps_fast.cfg --num_callers 2
```
**Demultiplex**

Because we multiplexed our samples to sequence multiple samples in one run, we will now run a program that computationally separates our samples based on the barcode sequence.
We'll start by making a directory for the split data:
```
mkdir -p $GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/basecalled/split

/home/groups/astraigh/software/ont-guppy-cpu/bin/guppy_barcoder --input_path "$GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/basecalled/pass" --save_path "$GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/basecalled/split" --config /home/groups/astraigh/software/ont-guppy-cpu/data/barcoding/configuration.cfg -t 2
```
We can now change to the output directory & list the contents to see what guppy output
```
cd $GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/basecalled/split
ls 
```
To combine the output fastqs, we'll first create a merged folder & then cat the multiple fastqs together. 
```
mkdir -p /split/barcode*/merged
cd /barcode20/merged
cat ./*.fastq > "$GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/basecalled/split/barcode20/merged/final.fastq"
```

## BasicQC
As we did earlier for our downloaded data, we'll take a look at the quality of our sequencing run.

make a directory for the qc analysis:
```
cd $GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/
mkdir -p qc
```
Create histogram and save it into a file. 
Replace my_fastq with the output file name from guppy.
```
$GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/basecalled/guppy/split/barcode20/merged/final.fastq | awk '(NR%4==2){print length($0)}' | sort | uniq -c | sort -k1,1nr > qc/readlength.hist.txt

less qc/readlength.hist.txt
```
Using FastQC
```
fastqc basecalled/guppy/split/barcode20/merged/final.fastq -o qc/ -t 2 
```
Transfer the output file to the Downloads folder on your computer. Be sure to change the directory to match your fastqc output.
```
rsync -ah --progress <username>@dtn.sherlock.stanford.edu:/scratch/groups/astraigh/biochem_minicourse_2021/teamStraight/data/sampleX/qc/my_fastqc.html ~/Downloads
```
Run on your local terminal:
```
open ~/Dowloads/my_fastqc.html
```
#Using Nanostat
```
NanoStat --fastq basecalled/guppy/split/barcode20/merged/final.fastq -o qc -n nanostat.summary
```

## Blast

Blast takes a fasta as an input, so we need to convert the fastq into fasta.

```bash
#recall your folder name
export me="teamStraight"

#move to the folder for this dataset
cd $GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX/

awk '(NR%4==1){print ">"$0; getline; print $0}' basecalled/guppy/split/barcode20/merged/final.fastq > basecalled/guppy/split/barcode20/merged/mysample.fasta
```

Let's look at the first 50 reads first
```
head -n 100 basecalled/guppy/split/barcode20/merged/mysample.fasta > basecalled/guppy/split/barcode20/merged/mysample.head.fasta
```

We have multipe databases installed on sherlock, which correspond to some of those available on ncbi blast online. Here, because we know we what type of samples we are dealing with we have downloaded just a few databases and combined them into one blastdb. The path of this database is  
```
/scratch/groups/astraigh/biochem_minicourse_2021/shared/blastdb/full_db.nal
```

Let's run blast using this database. We will put the blast output in a `blast` folder

```bash
#create folder for output
cd $GROUP_SCRATCH/biochem_minicourse_2021/$me/data/sampleX
mkdir -p blast

#We can set a variable for this blast database so we don't have to type out the full path eachtime
export BLASTDB=/scratch/groups/astraigh/biochem_minicourse_2021/shared/blastdb/full_db.nal

#run blast
blastn -query basecalled/guppy/split/barcode20/merged/mysample.head.fasta -db BLASTDB -num_threads 4 > blast/mysample.rawoutput.first50.txt
```

Taking a look at this file with `less` shows that this command returns the blast hits in a format similar to that of we ontain when running blast online. 
```
less blast/mysample.rawoutput.first50.txt
```

Note we can use "process substitution" with `<()` to get the first 50 reads on the fly without having to create an intermediate file

```
blastn -query <(head -n 100 basecalled/guppy/split/barcode20/merged/mysample.fasta) -db BLASTDB -num_threads 4 > blast/mysample.rawoutput.first50.txt
``` 

We can instead ask blast to produce a tabulated output, which is more computationally friendy, using the `--outfmt 6` flag.

```
blastn -query <(head -n 100 basecalled/guppy/split/barcode20/merged/mysample.fasta) -db BLASTDB -num_threads 4 -subject_besthit -outfmt "6 qseqid sseqid pident length mismatches gapopen evalue bitscore ssciname staxid sskingdom" > blast/mysample_tabulated.output.first50.txt
```


