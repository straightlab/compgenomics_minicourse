conda activate /share/PI/astraigh/miniconda3/envs/bootcamp

ROOTDIR="/scratch/groups/astraigh/biochem_minicourse_2021/practice_run" 

cd /scratch/groups/astraigh/
chmod -R g+rwx biochem_minicourse_2021
setfacl -m g:astraigh:rwX biochem_minicourse_2021
setfacl -d -m g:astraigh:rwX biochem_minicourse_2021

mkdir -p "${ROOTDIR}/dataset1/guppy"

## data download
fasterq-dump SRR8351023 --progress --threads 2 --temp "${ROOTDIR}/tmp/" --outdir "${ROOTDIR}/raw/"

##
fast5s = "/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial5/full_fast_gpu/fast5_split/fast5_eg"

/home/groups/astraigh/software/ont-guppy-cpu/bin/guppy_basecaller --input_path $fast5s --save_path "${ROOTDIR}/dataset1/guppy" --config /home/groups/astraigh/software/ont-guppy-cpu/data/dna_r9.4.1_450bps_fast.cfg --num_callers 2

## demultiplex
mkdir -p "${ROOTDIR}/dataset1/guppy/split"

/home/groups/astraigh/software/ont-guppy-cpu/bin/guppy_barcoder --input_path "${ROOTDIR}/dataset1/guppy/pass" --save_path "${ROOTDIR}/dataset1/guppy/split" --config /home/groups/astraigh/software/ont-guppy-cpu/data/barcoding/configuration.cfg -t 2

mkdir -p "${ROOTDIR}/dataset1/guppy/split/barcode20/merged"

cd "${ROOTDIR}/dataset1/guppy/split/barcode20"

cat ./*.fastq >"${ROOTDIR}/dataset1/guppy/split/barcode20/merged/final.fastq"

cd merged
NanoStat --fastq final.fastq -o merged -n nanostat.summary

fastqc final.fastq -o . -t 2

## Make a blastdb


cd /oak/stanford/groups/astraigh/T2T/PROseq

/home/groups/astraigh/software/ncbi-blast-2.10.1+/bin/makeblastdb -in chm13.draft_v1.0.fasta -dbtype nucl -parse_seqids

blastdb_aliastool -dblist "aiptasia bacteria b_impatiens c_hoffmanni human mouse racoon salmon shrimp symbio xenopus" -dbtype nucl -title Fulldb_nochm13 -out full_final

## Make bowtie2 db
ml biology bowtie2
bowtie2-build -f --threads 16 chm13.draft_v1.0.fasta chm13.draft_v1.0
