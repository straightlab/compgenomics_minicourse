# Intro to long read data

## Goals
In this notebook, we're going to download some sequencing data from the Gene Expression Omnibus https://www.ncbi.nlm.nih.gov/geo/, which is a public genomics data repository . For the purpose of familiarizing ourself with some simple long read sequencing data and metagenomics analysis, we will use some of the dataset reported in the Sevin et al. paper (Shotgun metagenome data of a defined mock community using Oxford Nanopore, PacBio and Illumina technologies). They did MinIon sequencing of synthetic microbial community `BMock12` consisting of 12 known 12 bacterial strains. Our goal is to recover the identity of these strains starting from the (almost) raw data. The dataset accession number for the ONT reads without size selection is SRX5161985 (https://www.ncbi.nlm.nih.gov/sra/SRX5161985) 

Covered tools and concepts in this notebook:
- fastq files
- use basic commands to view files `cat`,`head` and `less`
- use `wc` to count lines and very sample `awk` scripts to perform line by line operations on a file
- use `sort` and `uniq` to create histograms
- chain linux commands with `|` and redirect output to a file with `>`
- run a command line program such as `fastqc` to process the data in some manner
- download data back from Sherlock to your local machine with `rsync`

## GEO data download 
Data stored on GEO can be downloaded using `fasterq-dump` tool.
We'll download the raw data within the subfolder `raw` within the folder `data/woyke_mockcommunity` we've created earlier. Always keep things organized to avoid headaches down the road!

```bash
#move to the folder for this dataset
cd $PI_SCRATCH/bootcamp2020/teamCKO/data/woyke_mockcommunity

#make a new directory for the raw data
mkdir raw

# download the data
fasterq-dump SRR8351023 --progress --threads 2 --temp ./tmp/ --outdir ./raw
```

Let's do some quick inspection of the data we downloaded.
```
ls -lh raw
```
We can see the previous command just downloaded a file called SRR8351023.fastq wich is 1.8Gb in size.

This is a fastq formatted file, which is the standard format for sequencing data. We can figure out how the file is organized by looking at the first 10 lines.

```
head raw/SRR8351023.fastq
```

This is of the form 
```text
<@READID (first read) and some other read info>
<DNA sequence>
<+repeat of the READID line>
<Base call quality>
<@READID (next read) and some other read info>
etc..
```
So each read takes up 4 lines, with the sequence for read 1 at line2, sequence for read2 at line 4+2=6, etc... 

## Basic qc 
### The quick manual way (optional)
A lot of basic analysis of sequencing data can be done without relying on external software. It is extremely useful to know a few linux built-in commands that will come in handy in many contexts and especially for parsing raw sequencing data. 4 commands that are used ofter are `wc`, `sort`, `awk`, `uniq` (if you want to learn commands more you should also look at `cut`, `join` and `paste`). Let's show and example of real world application for these commands by doing the following qc analysis:
- how many reads are present in our dataset
- what is the read length distribution
- 
We can see how many reads are in this file by counting the number of lines and dividing by 4

```
wc -l raw/SRR8351023.fastq
# 579904 SRR8351023.fastq
# 579904/4 = 144,976 reads
```

What is the length of each read? We need to count the number of character for each sequence line. First let's  extract the sequences from the file. `awk` is a great unix program to manipulate text files line by line. Here we just filter lines whit a number modulo 4 = 2. Let's extract the first 3 sequences. The pipe operator `|` allows us to chain commands. The parenthesis in awk serves as an if statement. $0 is a variable containing the whole line.

```
head -n12 raw/SRR8351023.fastq | awk '(NR%4==2){print $0}'
```

Now to get the length of these first three sequences, we just print the lengh of the line instead of the line itself

```
head -n12 raw/SRR8351023.fastq | awk '(NR%4==2){print length($0)}'
```

Finally, we want to but a histogram of the how many times each read length is represented. For that, we need two useful commands: sort and uniq. Sort will just sort the lines, and uniq (wich requires sorted input), will count how many times each unique line occurs. We're ready to that for the whole file rather than the first three sequences so let's replace `head` with `cat`, which just reads through the whole file

```
cat raw/SRR8351023.fastq | awk '(NR%4==2){print length($0)} | sort | uniq -c 
```

We can redirect the output to a file with the > operator. Let's put call this file `readlength.hist.txt` and put it in a dedicated `qc` folder One last command we may want to chain in an other sort command so that the histogram is sorted so that the mode of the distribution comes first (the length that shows up the most often). This is achieved with 

```bash
# make a directory for the qc analysis
mkdir -p qc

#create histogram and save it into a file
cat raw/SRR8351023.fastq | awk '(NR%4==2){print length($0)} | sort | uniq -c | sort -k1,1nr > qc/readlength.hist.txt
```

We can look at this file in a scrollable way with `less`.
``` 
less qc/readlength.hist.txt
```

### QC with a dedicated tool fastqc and downloading data back
The `fastqc` tool (preinstalled on the lab partition) can be used to get some other qc metrics, in particular about sequencing quality and overrepresentated sequences.

```
fastqc raw/SRR8351023.fastq -o qc/ -t 2 
```

This produced an html file, which we need to download to our computer to look at. File transfer operations to and from Sherlock are best done using `rsync`

Open a new terminal tab in your computer, and then run 
```
rsync -avh --progress dtn.sherlock.stanford.edu:/scratch/groups/astraigh/bootcamp2020/teamCKO/external_data/woyke_mockcommunity/qc/*.html ~/Downloads

```
The command has the form `rsync <option flags> [source path] [destination path]`.  The optional flags -vh and --progress are just to tune the behavior of rsync and tell it to display progress in a nice way (-vh --prgress) and -a is to preserve timestamps on files. 

The source in on sherlock, and more specifically for file transfer we want to use a dedicated file transfer node on sherlock `dtn`. The targt path is just your local Download directory.

Now let's take a look. On you local terminal, run
```
open ~/Dowloads/SRR8351023_fastqc.html
```

### NanoStat
Another simple tool which produces qc metrics more relevant to long read data is `NanoStat`

```bash
NanoStat --fastq raw/SRR8351023.fastq -o qc -n nanostat.summary
```

Questions:
- What is the N50 for this dataset?
- What is the median Q score
- what is the mean error rate?

