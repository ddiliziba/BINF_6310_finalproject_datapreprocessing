
if docker run --help 2>/dev/null | grep -q -- '--platform'; then
  PLATFORM="--platform=linux/amd64"
else
  PLATFORM=""
fi

alias d_fastqc='docker run --rm $PLATFORM -v "$PWD":/data -w /data quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 fastqc'
alias d_multiqc='docker run --rm $PLATFORM -v "$PWD":/data -w /data quay.io/biocontainers/multiqc:1.24.1--pyhdfd78af_0 multiqc'
alias d_cutadapt='docker run --rm $PLATFORM -v "$PWD":/data -w /data quay.io/biocontainers/cutadapt:5.4--py311h3f08180_0 cutadapt'
alias d_star='docker run --rm $PLATFORM -v "$PWD":/data -w /data quay.io/biocontainers/star:2.7.11b--h43eeafb_2 STAR'
alias d_samtools='docker run --rm $PLATFORM -v "$PWD":/data -w /data staphb/samtools:1.20 samtools'
alias d_bcftools='docker run --rm $PLATFORM -v "$PWD":/data -w /data staphb/bcftools:1.20 bcftools'
alias d_bwa='docker run --rm $PLATFORM -v "$PWD":/data -w /data biocontainers/bwa:v0.7.17_cv1 bwa'
alias d_sra='docker run --rm $PLATFORM -v "$PWD":/data -w /data ncbi/sra-tools fasterq-dump'
alias d_gatk='docker run --rm $PLATFORM -v "$PWD":/data -w /data broadinstitute/gatk:4.6.2.0 gatk'
