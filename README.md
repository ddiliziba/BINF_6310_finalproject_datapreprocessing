# 6310 RNA processing (QC → Trim → Align) with Docker

## Prereqs
- Docker Desktop
- On Apple Silicon: Docker > Settings > Resources → Memory 24–32 GB (for STAR index)

## Layout
GroupProject/
  raw/    # put fastq.gz here (R1/R2)
  ref/    # put GRCh38.fasta and refseq.annotation.gtf here
  work/   # intermediates
  results/# reports/outputs

## Quickstart
make build                  # build the image
make qc                     # FastQC + MultiQC
make trim                   # cutadapt (produces work/NA12878/R*.trim.fq.gz)
make star-index             # one-time STAR genome index
make align                  # STAR → BAM + index
make where                  # show key outputs

### (Optional) Download GM12878 ENCODE pair
cd raw
curl -sL 'https://www.encodeproject.org/experiments/ENCSR000AED/?format=json' -o encsr.json
python3 ../scripts/select_fastqs.py > urls.txt
while read -r u; do curl -L --retry 5 --retry-connrefused --retry-delay 5 -C - -O "$u"; done < urls.txt
# rename to NA12878_GM12878_rep1_R1.fastq.gz / R2.fastq.gz

### Apple Silicon
Use: make build PLAT=--platform=linux/amd64
