ROOTDIR = "/scratch/groups/astraigh/biochem_minicourse_2021/practice_run" 

mkdir -p "${ROOTDIR}/dataset1/guppy"

## data download
fasterq-dump SRR8351023 --progress --threads 2 --temp "${ROOTDIR}/tmp/" --outdir "${ROOTDIR}/raw/"

##
fast5s = "/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial5/full_fast_gpu/fast5_split/fast5_eg"

/home/groups/astraigh/software/ont-guppy-cpu/bin/guppy_basecaller --input_path $fast5s --save_path "${ROOTDIR}/dataset1/guppy" --config /home/groups/astraigh/software/ont-guppy-cpu/data/dna_r9.4.1_450bps_fast.cfg --num_callers 2