# Taxonomic analysis

## Blast


## Running Kraken
Rather than running blasts, taxonomic classification can be done using dedicated tools. Kraken2 (https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) is a commonly used taxonomic sequence classifier, which examines the k-mers within a query sequence and compare these k-mers with those present in a taxonomic database. 

We have prebuild a taxonomic kmer database for kranken2 on Sherlock. This database contains RefSeq complete genomes for all archea, bacteria, virus, and the human genome. 

```
/home/groups/astraigh/software/KRAKEN2/kraken2 -db /scratch/groups/astraigh/kraken/standard --threads 6 --report kraken/standard.kreport --use-names --output - raw/SRR8351023.fastq
```

Take a look at the kraken report

```
less kraken/standard.kreport
```

Now we want to extract species level abundance and sort by the most abundant first

```
cut -f1,4,6 kraken/standard.kreport | awk '($2=="S")' | sort -k1,1nr | sed 's/  */ /g' > kraken/standard.species.txt
```

## Plotting abundance

In R

## Plotting normalized abundance


using Bracken


## Alignment

### Download genome
First we need to download the Muricauda_ruestringensis genome.


We get useful info here
https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/refseq_category:representative

https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastdbcmd_application_opti/


Only 1 team does that.
```bash
cd  $PI_SCRATCH/boootcamp2020/genomes

wget ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Muricauda_ruestringensis/representative/GCF_000224085.1_ASM22408v1/GCF_000224085.1_ASM22408v1_genomic.fna.gz

```

Unzip that file
```
gunzip GCF_000224085.1_ASM22408v1_genomic.fna.gz
```


### Alignement

minimap2 -x map-ont -d GCF_000224085.1_ASM22408v1_genomic_muricauda-rues.ont.mmi GCF_000224085.1_ASM22408v1_genomic.fna

minimap2 -a /scratch/groups/astraigh/bootcamp2020/genomes/GCF_000224085.1_ASM22408v1_genomic_muricauda-rues.ont.mmi raw/SRR8351023.fastq > alignments/ont.sam

head alignments/ont.sam

How many reads align?
samtools view -F 4 -F 256 -c alignments/ont.sam

How many reads don't align
samtools view -F 4 -c alignments/ont.sam

--> 23%

samtools view -b alignments/ont.sam >alignments/ont.bam

samtools sort -o alignments/ont.sorted.bam alignments/ont.bam

samtools index alignments/ont.sorted.bam

samtools view alignments/ont.sorted.bam | head

bedtools genomecov -bga -ibam alignments/ont.sorted.bam >alignments/ont.bedgraph

Question:
- take a quick look at the bedgraph file and guesstimate the genome coverage

Should we have them visualize that in IGV?



#
mkdir canu
canu -p canu -d canu -fast genomeSize=30m useGrid=False minThreads=2 maxThreads=2 maxMemory=30 -nanopore-raw raw/



ownloading nt database from server... done.
Uncompressing nt database...done.
Parsing nt FASTA file...done.
Masking low-complexity regions of downloaded library...srun: Force Terminated job 6792746
###
No 
readlength distribution --> manual


