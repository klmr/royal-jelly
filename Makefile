#
# Helpers
#

mkdir = mkdir -p '$(dir $@)'

bsub = scripts/bsub -K

keep = $(foreach i,$2,$(if $(findstring $1,$i),$i))

mem = -M$1 -R'select[mem>$1] rusage[mem=$1]'

raw-library-files = $(shell echo $(addsuffix *_L001_*,$(addprefix raw/,$(shell grep $1 raw/samples.csv | grep '^1' | cut -d, -f2 | tr _ -))))
library-files = $(subst _L001_,_merged_,$(addprefix $2,$(notdir $(call raw-library-files,$1))))
read-files = $(foreach i,R1 R2,$(shell sed 's,data/[^/]*/,data/trimmed/,' <<< '$(subst _001,_${i}_001,${1:.bam=.fastq.gz})'))

include binaries.make

#
# Helper variables
#

sequence-dir := ../shared/reference
index-dir := ../shared/index
annotation-dir = ../shared/annotation
apis-reference = ${sequence-dir}/apis_mellifera/Apis_mellifera.GCA_000002195.1.dna.toplevel.fa
homo-reference = ${sequence-dir}/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa
viral-reference = ${sequence-dir}/bee-viruses/bee_viruses.fasta

apis-annotation = ${annotation-dir}/apis_mellifera/Apis_mellifera.GCA_000002195.1.32.gtf
viral-annotation = ${annotation-dir}/bee-viruses/bee_viruses.gtf
apis-viral-annotation = ${annotation-dir}/$(notdir $(abspath ${apis-annotation}/..))-viruses/$(notdir ${apis-annotation:%.gtf=%-viruses.gtf})

apis-index = ${index-dir}/$(notdir $(basename ${apis-reference}))/Genome
homo-index = ${index-dir}/$(notdir $(basename ${homo-reference}))/Genome
viral-index = ${index-dir}/bee-viruses/Genome
apis-viral-index = ${index-dir}/$(notdir $(abspath ${apis-index}/..))-viruses/Genome

viral-references = $(foreach i,$(shell cut -d, -f5 supporting/bee-virus-list.csv),${sequence-dir}/bee-viruses/$i.fasta)
raw-reads = $(shell ls raw/*.fastq.gz)
merged-files = $(addprefix data/merged/,$(notdir $(subst L001,merged,$(call keep,L001,${raw-reads}))))
.PRECIOUS: ${merged-files}

long-trimmed-libraries = $(call library-files,long,data/trimmed/long/)
short-trimmed-libraries = $(call library-files,short,data/trimmed/short/)
trimmed-libraries = ${long-trimmed-libraries} ${short-trimmed-libraries}
.PRECIOUS: ${trimmed-libraries}

fastqc-files = $(patsubst %.fastq.gz,%_fastqc.zip,$(call library-files,long,data/qc/long/) $(call library-files,short,data/qc/short/))
read-lengths = $(subst .fastq.gz,.txt,$(subst /trimmed/,/qc/read-lengths/,${trimmed-libraries}))

mapped-reads = $(subst .fastq.gz,.bam,$(subst _R1_,_,$(call keep,_R1_,$(subst /trimmed/,/mapped/,${trimmed-libraries}))))
.PRECIOUS: ${mapped-reads}

mapped-indexed-reads = ${mapped-reads:.bam=.bam.bai}
.PRECIOUS: ${mapped-indexed-reads}

read-coverage = $(patsubst %.bam,%.genomecov,$(subst /mapped/,/coverage/,${mapped-reads}))
.PRECIOUS: ${read-coverage}

bigwig-files = $(patsubst %.bam,%.bw,$(subst /mapped/,/coverage/,${mapped-reads}))
.PRECIOUS: ${bigwig-files}

homo-mapped-reads = $(subst /mapped/,/human-mapped/,${mapped-reads})
.PRECIOUS: ${homo-mapped-reads}

feature-counts = $(patsubst %.bam,%.tsv,$(subst /mapped/,/quant/,${mapped-reads}))
.PRECIOUS: ${feature-counts}

insert-sizes = $(patsubst %.bam,%-insert-size.tsv,$(subst /mapped/,/isize/,${mapped-reads}))
.PRECIOUS: ${insert-sizes}

#
# Download and/or build the various reference genomes and annotations
#

${apis-reference}:
	@$(mkdir)
	wget ftp://ftp.ensemblgenomes.org/pub/release-32/metazoa/fasta/apis_mellifera/dna/$(notdir $@).gz -O '$@.gz'
	gunzip $@

${sequence-dir}/bee-viruses/%.fasta:
	@$(mkdir)
	./scripts/efetch-fasta '$*' '$@'

.PHONY: viral-reference
## Build the viral reference metagenome
viral-reference: ${viral-reference}
	
${viral-reference}: ${viral-references}
	cat $+ > '$@'

${homo-reference}:
	@$(mkdir)
	wget ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/$(notdir $@).gz -O '$@.gz'
	gunzip $@

${apis-annotation}:
	@$(mkdir)
	wget ftp://ftp.ensemblgenomes.org/pub/release-32/metazoa/gtf/apis_mellifera/$(notdir $@).gz -O '$@.gz'
	gunzip $@

${viral-annotation}: ${viral-references}
	@$(mkdir)
	./scripts/fasta-to-gtf -o $@ $+

.PHONY: annotation
## Merge the bee and bee virus annotation file
annotation: ${apis-viral-annotation}

${apis-viral-annotation}: ${apis-annotation} ${viral-annotation}
	@$(mkdir)
	cat $+ > '$@'

#
# Build mapping indices
#

define build-index
@$(mkdir)
${bsub} $(call mem,${MEM}) \
	"STAR --runMode genomeGenerate --genomeDir '$(dir $@)' \
	--genomeFastaFiles $+"; rm -rf _STARtmp Log.out
endef

${apis-index}: ${apis-reference}
	$(eval MEM=128000)
	$(build-index)

${homo-index}: ${homo-reference}
	$(eval MEM=48000)
	$(build-index)

${apis-viral-index}: ${apis-reference} ${viral-reference}
	$(eval MEM=128000)
	$(build-index)

${viral-index}: ${viral-reference}
	@$(mkdir)
	STAR --runMode genomeGenerate --genomeSAindexNbases 3 \
		--genomeDir '$(dir $@)' --genomeFastaFiles '$+'
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
# Trimming
#

.PHONY: trim-long
## Trim long RNA-seq libraries
trim-long: ${long-trimmed-libraries}

data/trimmed/long/%R1_001.fastq.gz: data/merged/%R1_001.fastq.gz data/merged/%R2_001.fastq.gz
	@$(mkdir)
	${bsub} -n3 "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
		--minimum-length 10 \
		-o '$@' -p '$(subst _R1_,_R2_,$@)' '$(firstword $+)' '$(lastword $+)' \
		> '${@:.fastq.gz=.log}'"

data/trimmed/long/%R2_001.fastq.gz: data/merged/%R1_001.fastq.gz data/merged/%R2_001.fastq.gz
	# $@ already created by the preceding rule.

.PHONY: trim-short
## Trim short RNA-seq libraries
trim-short: ${short-trimmed-libraries}

data/trimmed/short/%R1_001.fastq.gz: data/merged/%R1_001.fastq.gz data/merged/%R2_001.fastq.gz
	@$(mkdir)
	${bsub} -n3 "cutadapt -a NNNNTGGAATTCTCGGGTGCCAAGG -A NNNNGATCGTCGGACTGTAGAACTCTGAAC \
		--minimum-length 10 \
		-o '$@' -p '$(subst _R1_,_R2_,$@)' '$(firstword $+)' '$(lastword $+)' \
		> '${@:.fastq.gz=.log}'"

data/trimmed/short/%R2_001.fastq.gz: data/merged/%R1_001.fastq.gz data/merged/%R2_001.fastq.gz
	# $@ already created by the preceding rule.

#
# Map potential human contamination
#

.PHONY: map-to-human
## Map reads to human reference
map-to-human: ${homo-mapped-reads}

data/human-mapped/%.bam: $$(call read-files,$$@) ${homo-index}
	@$(mkdir)
	${bsub} -n 12 $(call mem,24000) \
		"STAR --runThreadN 12 --genomeDir '$(dir $(lastword $^))' \
		--runMode alignReads --alignEndsType Local \
		--outFilterMismatchNoverLmax 0.15 --outFilterMultimapNmax 1000 \
		--readFilesIn $(call read-files,$@) --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
		--outFileNamePrefix '$(basename $@)'"
	mv "$(basename $@)Aligned.sortedByCoord.out.bam" "$(basename $@).bam"

# FIXME: Align reads against viral/human contaminants with soft-clipping and 85% identity for match, allowing multi-mapping reads.
# Map against the human genome *a lot* more stringently because of conservation.
# Remove anything that maps from database, map the rest against apis.

#
# Read mapping
#

.PHONY: map-reads
## Map reads to reference of bee genome and viral genomes
map-reads: ${mapped-reads}

data/mapped/long/%.bam: $$(call read-files,$$@) ${apis-viral-index}
	@$(mkdir)
	${bsub} -n 12 $(call mem,24000) \
		"STAR --runThreadN 12 --genomeDir '$(dir $(lastword $^))' \
		--runMode alignReads --alignEndsType Local \
		--outFilterMismatchNoverLmax 0.15 --outFilterMultimapNmax 10 \
		--readFilesIn $(call read-files,$@) --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
		--outFileNamePrefix '$(basename $@)'"
	mv "$(basename $@)Aligned.sortedByCoord.out.bam" "$(basename $@).bam"

data/mapped/short/%.bam: $$(call read-files,$$@) ${apis-viral-index}
	@$(mkdir)
	${bsub} -n 12 $(call mem,24000) \
		"STAR --runThreadN 12 --genomeDir '$(dir $(lastword $^))' \
		--runMode alignReads --alignEndsType Local \
		--outFilterMatchNmin 18 \
		--outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 10000 \
		--readFilesIn $(call read-files,$@) --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
		--outFileNamePrefix '$(basename $@)'"
	mv "$(basename $@)Aligned.sortedByCoord.out.bam" "$(basename $@).bam"

#
# QC report
#

data/qc/%_fastqc.zip: data/trimmed/%.fastq.gz
	@$(mkdir)
	${bsub} $(call mem,4000) \
		"fastqc --outdir '$(dir $@)' '$<'"
	rm '${@:.zip=.html}'

.PHONY: qc-report
## Quality control report
qc-report: data/qc/multiqc_report.html

data/qc/multiqc_report.html: ${fastqc-files} ${trimmed-libraries} ${mapped-reads} ${feature-counts} ${insert-sizes}
	multiqc --force --outdir data/qc $(sort $(dir $+))

.PHONY: insert-sizes
## Compute insert size distribution and plot their density
insert-sizes: ${insert-sizes}

data/isize/%-insert-size.tsv: data/mapped/%.bam
	@$(mkdir)
	${bsub} $(call mem,16000) \
		"picard CollectInsertSizeMetrics \
		INPUT='$<' \
		OUTPUT='$@' \
		HISTOGRAM_FILE='${@:.tsv=.pdf}' \
		ASSUME_SORTED=false"

.PHONY: read-lengths
## Compute length distributions and plot their density
read-lengths: data/qc/read-lengths/read-length-density.pdf

${read-lengths}: ${trimmed-libraries}

data/qc/read-lengths/%.txt: data/trimmed/%.fastq.gz bin/line-lengths
	@$(mkdir)
	${bsub} "gunzip -c '$<' | sed -n 2~4p | bin/line-lengths > '$@'"

data/qc/read-lengths/read-length-density.pdf: ${read-lengths}
	${bsub} $(call mem,4000) \
		"./scripts/plot-read-length-distribution '$@' $+"

#
# Visualise mapping
#

data/mapped/%.bam.bai: data/mapped/%.bam
	${bsub} "samtools index '$<'"

.PHONY: read-coverage
## Compute per-base read coverage over the bee/viral genomes
read-coverage: ${read-coverage}

data/coverage/%.genomecov: data/mapped/%.bam
	@$(mkdir)
	${bsub} $(call mem,4000) \
		"bedtools genomecov -d -ibam '$<' > '$@'"

.PHONY: bigwig
## Generate BigWig files of the bee/viral read coverage
bigwig: ${bigwig-files}

data/coverage/%.bedgraph: data/mapped/%.bam
	@$(mkdir)
	${bsub} $(call mem,1000) "./scripts/sorted-bedgraph '$<' '$@'"

data/coverage/%.bw: data/coverage/%.bedgraph ${apis-viral-index}
	${bsub} $(call mem,1000) \
		"bedGraphToBigWig '$<' '${apis-viral-index:Genome=chrNameLength.txt}' '$@'"

#
# Feature counts
#

.PHONY: feature-counts
## Count read coverage of genes in bee and viruses
feature-counts: ${feature-counts}

data/quant/%.tsv: data/mapped/%.bam ${apis-viral-annotation}
	@$(mkdir)
	${bsub} $(call mem,1000) \
		"featureCounts -p -t gene -g gene_id -M -O -a '$(lastword $+)' -o '$@' '$<'"

data/quant/royal-jelly-counts.tsv: ${feature-counts}
	./scripts/merge-counts -o '$@' $+

.DEFAULT_GOAL := show-help
# See <https://gist.github.com/klmr/575726c7e05d8780505a> for explanation.
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)";echo;sed -ne"/^## /{h;s/.*//;:d" -e"H;n;s/^## //;td" -e"s/:.*//;G;s/\\n## /---/;s/\\n/ /g;p;}" ${MAKEFILE_LIST}|LC_ALL='C' sort -f|awk -F --- -v n=$$(tput cols) -v i=19 -v a="$$(tput setaf 6)" -v z="$$(tput sgr0)" '{printf"%s%*s%s ",a,-i,$$1,z;m=split($$2,w," ");l=n-i;for(j=1;j<=m;j++){l-=length(w[j])+1;if(l<= 0){l=n-i-length(w[j])-1;printf"\n%*s ",-i," ";}printf"%s ",w[j];}printf"\n";}'|more $(shell test $(shell uname) == Darwin && echo '-Xr')
