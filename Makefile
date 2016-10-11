#
# Helpers
#

mkdir = mkdir -p '$(dir $@)'

bsub = scripts/bsub -K

keep = $(foreach i,$2,$(if $(findstring $1,$i),$i))

raw-library-files = $(shell echo $(addsuffix *_L001_*,$(addprefix raw/,$(shell grep $1 raw/samples.csv | cut -d, -f2 | tr _ -))))
library-files = $(subst _L001_,_merged_,$(addprefix $2,$(notdir $(call raw-library-files,$1))))
read-files = $(foreach i,R1 R2,$(subst /mapped/,/trimmed/,$(subst _001,_${i}_001,${1:.bam=.fastq.gz})))

#
# Helper variables
#

sequences := ../shared/reference
index-dir := ../shared/index
apis-reference = ${sequences}/apis_mellifera/Apis_mellifera.GCA_000002195.1.dna.toplevel.fa.gz
homo-reference = ${sequences}/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
viral-reference = ${sequences}/bee-viruses/bee_viruses.fasta
viral-index = ${index-dir}/bee-viruses/Genome
index = ${index-dir}/bee-contamination/Genome
viral-references = $(foreach i,$(shell cut -d, -f5 supporting/bee-virus-list.csv),${sequences}/bee-viruses/$i.fasta)
raw-reads = $(shell ls raw/*.fastq.gz)
merged-files = $(addprefix data/merged/,$(notdir $(subst L001,merged,$(call keep,L001,${raw-reads}))))
.PRECIOUS: ${merged-files}
long-trimmed-libraries = $(call library-files,long,data/trimmed/long/)
short-trimmed-libraries = $(call library-files,short,data/trimmed/short/)
.PRECIOUS: ${long-trimmed-libraries}
.PRECIOUS: ${short-trimmed-libraries}
fastqc-files = $(patsubst %.fastq.gz,%_fastqc.zip,$(call library-files,long,data/qc/long/) $(call library-files,short,data/qc/short/))

mapped-reads = $(subst .fastq.gz,.bam,$(subst _R1_,_,$(call keep,_R1_,$(subst /trimmed/,/mapped/,${long-trimmed-libraries}, ${short-trimmed-libraries}))))
.PRECIOUS: ${mapped-reads}

#
# Download and/or build the various reference genomes
#

${apis-reference}:
	@$(mkdir)
	wget ftp://ftp.ensemblgenomes.org/pub/release-32/metazoa/fasta/apis_mellifera/dna/$(notdir $@) -O '$@'

${sequences}/bee-viruses/%.fasta:
	@$(mkdir)
	./scripts/efetch-fasta '$*' '$@'

.PHONY: viral-reference
## Build the viral reference metagenome
viral-reference: ${viral-reference}
	
${viral-reference}: ${viral-references}
	cat ${sequences}/bee-viruses/*.fasta > '$@'

${homo-reference}:
	@$(mkdir)
	wget ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/$(notdir $@) -O '$@'

#
# Build mapping indices
#

${viral-index}: ${viral-reference}
	@$(mkdir)
	STAR --runMode genomeGenerate --genomeSAindexNbases 3 \
		--genomeDir '$(dir $@)' --genomeFastaFiles '$+'
	rm Log.out

${index}: ${apis-reference} ${viral-reference} ${homo-reference}
	@$(mkdir)
	${bsub} -M24000 -R'select[mem>24000] rusage[mem=24000]' \
		"STAR --runMode genomeGenerate --genomeDir '$(dir $@)' \
		--genomeFastaFiles '$+'"
	rm Log.out

#
# First step in the input processing: merge files across lanes.
#

.SECONDEXPANSION:

lane-files = $(foreach i,L001 L002,raw/$(notdir $(subst merged,$i,$1)))

merged-files: ${merged-files}
data/merged/%.fastq.gz: $$(call lane-files,$$@)
	@$(mkdir)
	${bsub} "./scripts/merge-lanes $+ '$@'"

#
# Preliminary QC report on the unmapped input files
#

.PHONY: trim-long
## Trim long RNA-seq libraries
trim-long: ${long-trimmed-libraries}

data/trimmed/long/%R1_001.fastq.gz: data/merged/%R1_001.fastq.gz data/merged/%R2_001.fastq.gz
	@$(mkdir)
	${bsub} "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
		-o '$@' -p '$(subst _R1_,_R2_,$@)' '$(firstword $+)' '$(lastword $+)'"

data/trimmed/long/%R2_001.fastq.gz: data/merged/%R1_001.fastq.gz data/merged/%R2_001.fastq.gz
	# $@ already created by the preceding rule.

.PHONY: trim-short
## Trim short RNA-seq libraries
trim-short: ${short-trimmed-libraries}

data/trimmed/short/%R1_001.fastq.gz: data/merged/%R1_001.fastq.gz data/merged/%R2_001.fastq.gz
	@$(mkdir)
	${bsub} "cutadapt --cut 4 -a NNNNTGGAATTCTCGGGTGCCAAGG -A NNNNGATCGTCGGACTGTAGAACTCTGAAC \
		-o '$@' -p '$(subst _R1_,_R2_,$@)' '$(firstword $+)' '$(lastword $+)'"

data/trimmed/short/%R2_001.fastq.gz: data/merged/%R1_001.fastq.gz data/merged/%R2_001.fastq.gz
	# $@ already created by the preceding rule.

data/qc/%_fastqc.zip: data/trimmed/%.fastq.gz
	@$(mkdir)
	${bsub} -M4000 -R'select[mem>4000] rusage[mem=4000]' \
		"fastqc --outdir '$(dir $@)' '$<'"
	rm '${@:.zip=.html}'

.PHONY: qc-report
## Quality control report
qc-report: data/qc/multiqc_report.html

data/qc/multiqc_report.html: ${fastqc-files}
	multiqc --force --outdir data/qc data/qc

#
# Remove viral contamination
#

read-files = $(foreach i,xyz,raw/$(notdir ${1:.bam=$i.fastq.gz}))

data/viral-mapped/%.bam: $$(call read-files,$$@) ${viral-index}
	@$(mkdir)
	${bsub} -n 12 -M24000 -R'select[mem>24000] rusage[mem=24000]' \
		"STAR --runThreadN 12 --genomeDir '$(dir $(lastword $^))' \
		--runMode alignReads --alignEndsType Local \
		--outFilterMismatchNoverLmax 0.15 --outFilterMultimapNmax 1000 \
		--readFilesIn $(call read-files,%) --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM Unsorted --outFileNamePrefix '$(basename $@)'"
	mv "$(basename $@)Alignment.out.bam" "$(basename $@).bam"

# FIXME: Align reads against viral/human contaminants with soft-clipping and 85% identity for match, allowing multi-mapping reads.
# Map against the human genome *a lot* more stringently because of conservation.
# Remove anything that maps from database, map the rest against apis.

.DEFAULT_GOAL := show-help
# See <https://gist.github.com/klmr/575726c7e05d8780505a> for explanation.
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)";echo;sed -ne"/^## /{h;s/.*//;:d" -e"H;n;s/^## //;td" -e"s/:.*//;G;s/\\n## /---/;s/\\n/ /g;p;}" ${MAKEFILE_LIST}|LC_ALL='C' sort -f|awk -F --- -v n=$$(tput cols) -v i=19 -v a="$$(tput setaf 6)" -v z="$$(tput sgr0)" '{printf"%s%*s%s ",a,-i,$$1,z;m=split($$2,w," ");l=n-i;for(j=1;j<=m;j++){l-=length(w[j])+1;if(l<= 0){l=n-i-length(w[j])-1;printf"\n%*s ",-i," ";}printf"%s ",w[j];}printf"\n";}'|more $(shell test $(shell uname) == Darwin && echo '-Xr')
