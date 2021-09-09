# Alignment

In this notebook, we will align our reads to the genomes of the top 19 identified species.

## Download genome files

Firs we need to download the relevant genomes.

Information about the genomes in the kraken database are found in a file called `assembly_summary.txt`. For example, for the bacterial library (a subset of the kraken database), this is `/scratch/groups/astraigh/kraken/standard/library/bacteria/assembly_summary.txt`


To find info on a specific genome in that file, for example cohaesibacter, we just look for the word "cohaesibacter". `grep` is a linux tool for searching plain-text data sets for lines that contain a query pattern (here the word coahesibacter). The `-i` flag allows for case-insensitive search. 

```
grep -i "cohaesibacter" /scratch/groups/astraigh/kraken/standard/library/bacteria/assembly_summary.txt
```

We simplity this file to retain only the following columns
- 6: taxid (can be either species level or strain)
- 7: species level taxid
- 8: taxon name 
- 20: path to file on ncbi ftp
- 21: wether it's excluded from refseq or not


```
grep -v "#" /scratch/groups/astraigh/kraken/standard/library/bacteria//assembly_summary.txt | cut -f6,7,8,20,21 > /scratch/groups/astraigh/kraken/standard/library/bacteria/genomes.wget.txt
```

Now we join this file with the top 19 species to get the path  download their genome

```
cd /scratch/groups/astraigh/genomics_minicourse/shared/genomes/mock12

awk -F $'\t' 'BEGIN{OFS=FS}(NR==FNR){x[$5]=$6; next;}{if($2 in x){split($4, y, "/"); print $1, $2, $3, $4"/"y[length(y)]"_genomic.fna.gz", $5, x[$2];}}' <(head -n 19 kraken/standard.species.txt) /scratch/groups/astraigh/kraken/standard/library/bacteria/genomes.wget.txt  >top19.genomes.txt
```

Finally we can download all the genomes using wget
```
module load system parallel

parallel -j1 'wget {1} -P /scratch/groups/astraigh/genomics_minicourse/shared/genomes/mock12' :::: <(cut -f4,4 top19.genomes.txt)
```

And concatenate them in a single fasta
```
cd /scratch/groups/astraigh/genomics_minicourse/shared/genomes/mock12
zcat ./*.fna.gz >>/scratch/groups/astraigh/genomics_minicourse/shared/genomes/mock12/all_genomic.fna
```

So all our genomes are in `/scratch/groups/astraigh/genomics_minicourse/shared/genomes/mock12/all_genomic.fna`

Note this fasta is annotted with accession numbers rather than taxonomic id (and some species have multiple scaffolds), so we need a conversion file to link accession number to taxid.

```
awk -F $'\t' '{split($1,x,"|"); print(x[3],x[2]);}' /scratch/groups/astraigh/kraken/standard/library/bacteria/seqid2taxid.map >/scratch/groups/astraigh/kraken/standard/library/bacteria/scaffold2taxid.map
```


## Alignement

We will align our reads with the minimap2 aligner. First let's make an index for the aligner

```bash
cd /scratch/groups/astraigh/genomics_minicourse/shared/data/mock12

#mk a dir for the output
mkdir -p alignments

#create index
minimap2 -a ../../genomes/mock12/all_genomic.ont.mmi raw/SRR8351023.fastq >alignments/ont.sam
```

Then align
```
minimap2 -t 6 -a ../../genomes/mock12/all_genomic.ont.mmi raw/SRR8351023.fastq >alignments/ont.sam
```

We can quickly check :
- How many reads align?
```
samtools view -F 4 -F 256 -c alignments/ont.sam
```

- How many reads don't align
```
samtools view -F 4 -c alignments/ont.sam
```

Finally we convert the sam into bam (binary format) and sort
```bash 
samtools view -b alignments/ont.sam >alignments/ont.bam
samtools sort -o alignments/ont.sorted.bam alignments/ont.bam
samtools index alignments/ont.sorted.bam
```


## Coverage

We can get genome coverage with `bedtools`

```bash
bedtools genomecov -ibam alignments/ont.sorted.bam >alignments/ont.coverage.txt
```


