PLAT ?= --platform=linux/amd64
IMG  ?= hap-phase:latest
PWD  := $(shell pwd)

build:
	@docker build $(PLAT) -t $(IMG) .

shell:
	@docker run --rm -it $(PLAT) -v "$(PWD)":/data -w /data $(IMG)

qc:
	@docker run --rm $(PLAT) -v "$(PWD)":/data -w /data $(IMG) \
	'mkdir -p results/rawQC && fastqc -t 8 -o results/rawQC raw/*.fastq.gz && multiqc -o results results/rawQC'

trim:
	@docker run --rm $(PLAT) -v "$(PWD)":/data -w /data $(IMG) \
	'mkdir -p work/NA12878 && cutadapt -j 8 -q 20,20 -m 30 \
	 -o work/NA12878/R1.trim.fq.gz -p work/NA12878/R2.trim.fq.gz \
	 raw/NA12878_GM12878_rep1_R1.fastq.gz raw/NA12878_GM12878_rep1_R2.fastq.gz'

star-index:
	@docker run --rm $(PLAT) -v "$(PWD)":/data -w /data $(IMG) \
	'mkdir -p ref/STAR_GRCh38 && samtools faidx ref/GRCh38.fasta && \
	 STAR --runThreadN 8 --runMode genomeGenerate \
	     --genomeDir ref/STAR_GRCh38 \
	     --genomeFastaFiles ref/GRCh38.fasta \
	     --sjdbGTFfile ref/refseq.annotation.gtf \
	     --limitGenomeGenerateRAM 32000000000 --genomeSAindexNbases 12'

align:
	@docker run --rm $(PLAT) -v "$(PWD)":/data -w /data $(IMG) \
	'STAR --runThreadN 8 --genomeDir ref/STAR_GRCh38 \
	      --readFilesIn work/NA12878/R1.trim.fq.gz work/NA12878/R2.trim.fq.gz \
	      --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
	      --twopassMode Basic --sjdbGTFfile ref/refseq.annotation.gtf \
	      --outSAMattrRGline ID=NA12878_RNA SM=NA12878 PL=ILLUMINA LB=LIB1 \
	      --outFileNamePrefix work/NA12878/rna_ \
	 && samtools index work/NA12878/rna_Aligned.sortedByCoord.out.bam'

rna: qc trim align

where:
	@echo "FastQC/MultiQC: results/rawQC/, results/multiqc_report.html"
	@echo "Trimmed FASTQs: work/NA12878/R1.trim.fq.gz, R2.trim.fq.gz"
	@echo "RNA BAM:        work/NA12878/rna_Aligned.sortedByCoord.out.bam(.bai)"
