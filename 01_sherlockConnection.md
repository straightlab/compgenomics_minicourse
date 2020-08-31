# Connection to Sherlock and familiarization with the linux shell 

We will do all of our bioinformatics analysis on the Stanford High Performance Computer Cluster Sherlock. This requires basic knowledge of the linux bash shell. In this part, we will 
- connect to the cluster
- setup up our working directory
- download some external sequencing dataset
- do some very basic qc on this dataset and learn some of the most useful linux commands along the way! 

At the end of this notebook you will have had a glimpse at how to :
- navigate the linux filesystem with `ls, cd, mkdir`
- use basic commands to view files `cat`,`head` and `less`
- use `wc` to count lines and very sample `awk` scripts to perform line by line operations on a file
- use `sort` and `uniq` to create histograms
- run a command line program such as `fastqc` to process the data in some manner
- download data back from Sherlock to your local machine with `rsync`


## Connection to Sherlock
In your browser go to `login.sherlock.stanford.edu`, then click on `>Shell`
Alternatively, you can ssh into sherlock using a terminal app on your computer. On a Mac, you can use the native terminal app Term. Open up Term (Command-space term), then
`ssh <username>@login.sherlock.stanford.edu`

## Setting up our workspace for the project.
Go to the bootcamp directory
```
cd $SCRATCH/bootcamp2020
```
Note that $SCRATCH is a bash variable, which contains the path to a default "scratch" folder for our lab. You can see the content of this variable with
```
echo $SCRATCH
```

Make a new folder for your team, and move to that folder. For example for team CKO
```
mkdir teamCKO
cd teamCKO
```

This folder is currently empty. You can list the content of the folder with 
```
ls -lah teamCKO
```

## Downloading a sequencing dataset
Now we're going to download some sequencing data from the Gene Expression Omnibus https://www.ncbi.nlm.nih.gov/geo/, which is a public genomics data repository . For the purpose of familiarizing ourself with some simple long read sequencing data and metagenomics analysis, we will use some of the dataset reported in the Sevin et al. paper (Shotgun metagenome data of a defined mock community using Oxford Nanopore, PacBio and Illumina technologies). They did MinIon sequencing of synthetic microbial community `BMock12` consisting of 12 known 12 bacterial strains. Our goal is going to be to recover the identity of these strains starting from the (almost) raw data. The dataset accession number for the ONT reads without size selection is SRX5161985 (https://www.ncbi.nlm.nih.gov/sra/SRX5161985) 

```
mkdir -p external/woyke_mockcommunity/raw
fasterq-dump SRR8351023 --progress --threads 6 --temp ./tmp/ --outdir ./external/woyke_mockcommunity/raw
```

Let's do some quick inspection of the data we downloaded.
```
cd /external/woyke_mockcommunity/raw
ls -lah
```
We can see the previous command just downloaded a file called SRR8351023.fastq wich is ~ Gb in size.

This is a fastq formatted file, which is the standard format for sequencing data. We can figure out how the file is organized by looking at the first 10 lines.

```
head SRR8351023.fastq
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
### The quick manual way
A lot of basic analysis of sequencing data can be done without relying on external sofware. If is extremely useful to know a few linux built-in commands that will come in handy in many context and especially for parsing raw sequencing data. Arguably the top 4 for a bioinformatician are `wc`, `sort`, `awk`, `uniq` (if you want to learn commands more you should also look at `cut`, `join` and `paste`). Let's show and example of real world application for these commands by doing the following qc analysis:
- how many reads are present in our dataset
- what is the read length distribution


We can see how many reads are in this file by counting the number of lines and dividing by 4

```
wc -l SRR8351023.fastq
# 579904 SRR8351023.fastq
# 579904/4 = 144,976 reads
```

What is the length of each read? We need to count the number of character for each sequence line. First let's  extract the sequences from the file. `awk` is a great unix program to manipulate text files line by line. Here we just filter lines whit a number modulo 4 = 2. Let's extract the first 3 sequences. The pipe operator allows us to chain commands. The parenthesis in awk serves as an if statement. $0 is a variable containing the whole line.

```
head -n12 SRR8351023.fastq | awk '(NR%4==2){print $0}'
```

Now to get the length of these first three sequences, we just print the lengh of the line instead of the line itself

```
head -n12 SRR8351023.fastq | awk '(NR%4==2){print length($0)}'
```

Finally, we want to but a histogram of the how many times each read length is represented. For that, we need two useful commands: sort and uniq. Sort will just sort the lines, and uniq (wich requires sorted input), will count how many times each unique line occurs. We're ready to that for the whole file rather than the first three sequences so let's replace `head` with `cat`, which just reads through the whole file

```
cat SRR8351023.fastq | awk '(NR%4==2){print length($0)} | sort | uniq -c 
```

We can redirect the output to a file with the > operator. Let's put call this file `readlength.hist.txt` and put it in a dedicated `qc` folder One last command we may want to chain in an other sort command so that the histogram is sorted so that the mode of the distribution comes first (the length that shows up the most often). This is achieved with 

```
cd /scratch/groups/astraigh/bootcamp2020/teamCKO/external_data/woyke_mockcommunity

mkdir -p qc

cat raw/SRR8351023.fastq | awk '(NR%4==2){print length($0)} | sort | uniq -c | sort -k1,1nr > qc/readlength.hist.txt
```

We can look at this file in a scrollable way with `less`.
``` 
less readlength.hist.txt
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
The command has the form `rsync <option flags> [source path] [destination path]`.  The optionnal flags -vh and --progress are just to tune the behavior of rsync and tell it to display progress in a nice way (-vh --prgress) and -a is to preserve timestamps on files. 

The source in on sherlock, and more specifically for file transfer we want to use a dedicated file transfer node on sherlock `dtn`. The targt path is just your local Download directory.

Now let's take a look. On you local terminal, run
```
open ~/Dowloads/SRR8351023_fastqc.html
```



