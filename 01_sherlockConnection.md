# Connection to Sherlock and familiarization with the linux shell

We will do all of our bioinformatics analysis on the Stanford High Performance Computer Cluster Sherlock. This requires basic knowledge of the linux bash shell. In this part, we will connect to the cluster and download some external data to familiarize ourself with the linux shell. 

In your browser go to `login.sherlock.stanford.edu`, then click on `>Shell`
Alternatively, you can ssh into sherlock using a terminal app on your computer. On a Mac, you can use the native terminal app Term. Open up Term (Command-space term), then
`ssh <username>@login.sherlock.stanford.edu`

Go to the bootcamp directory
`cd $SCRATCH/bootcamp2020/`
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
ls -la teamCKO
```

Now we're going to download some sequencing data from the Gene Expression Omnibus https://www.ncbi.nlm.nih.gov/geo/, which is a public genomics data repository . For the purpose of familiarizing ourself with some simple long read sequencing data and metagenomics analysis, we will use some of the dataset reported in the Sevin et al. paper (Shotgun metagenome data of a defined mock community using Oxford Nanopore, PacBio and Illumina technologies). They did MinIon sequencing of synthetic microbial community `BMock12` consisting of 12 known 12 bacterial strains. Our goal is going to be to recover the identity of these strains starting from the (almost) raw data. The dataset accession number for the ONT reads without size selection is SRX5161985 (https://www.ncbi.nlm.nih.gov/sra/SRX5161985) 

```
mkdir -p external/woyke_mockcommunity/raw
fasterq-dump SRR8351023 --progress --threads 6 --temp ./tmp/ --outdir ./external/woyke_mockcommunity/raw

cd /scratch/groups/astraigh/bootcamp2020/teamCKO/external_data/woyke_mockcommunity
```

# Running Kraken
We have prebuild a full taxonomic kmer db for kranken2 on Sherlock. 

```
/home/groups/astraigh/software/KRAKEN2/kraken2 -db /scratch/groups/astraigh/kraken/standard --threads 6 --report kraken/standard.kreport --use-names --output - raw/SRR8351023.fastq
```


