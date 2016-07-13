raw-files := $(shell ls raw/*.fastq.gz)
trimmed-files := $(addsuffix .fastq.gz,$(addprefix data/trimmed/,$(basename $(basename $(notdir ${raw-files})))))
lengths := $(addsuffix .lengths.tsv,${trimmed-files:.fastq.gz=})
species := Apis_mellifera
species_v := ${species}.GCA_000002195.1.31
index := data/index/${species_v}/Genome
reference := data/reference/${species_v}.dna.toplevel.fa
annotation := data/reference/${species_v}.gtf
rrna-annotation := ${annotation:.gtf=.rrna.gtf}
mapped-reads := $(addsuffix .bam,$(addprefix data/mapped/,$(basename $(basename $(notdir ${raw-files})))))
rrna-contamination := $(addsuffix .tsv,$(addprefix data/rrna-contamination/,$(basename $(basename $(notdir ${mapped-reads})))))

bsub := scripts/bsub -K -q research-rh7

.PHONY: reference
## Download the genome reference sequence
reference: ${reference}

${reference}:
	@mkdir -p "$(dir $@)"
	wget 'ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/apis_mellifera/dna/${species_v}.dna.toplevel.fa.gz' -O $@.gz
	gunzip $@

.PHONY: annotation
## Download the genome annotation
annotation: ${annotation}

${annotation}:
	@mkdir -p "$(dir $@)"
	wget 'ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/gtf/apis_mellifera/${species_v}.gtf.gz' -O $@.gz
	gunzip $@

${rrna-annotation}: ${annotation}
	fgrep 'gene_biotype "rRNA"' '$<' > '$@'

.PHONY: index
## Create the STAR index of the genome
index: ${index}

${index}: ${reference}
	@mkdir -p "$(dir $@)"
	${bsub} -n 12 -M16000 -R'select[mem>16000]' -R'rusage[mem=16000]' \
		"STAR --runThreadN 12 --runMode genomeGenerate \
		--genomeDir '$(dir $@)' --genomeFastaFiles '$<'"

.PHONY: trimmed-reads
## Perform read trimming
trimmed-reads: ${trimmed-files}

data/trimmed/%.fastq.gz: raw/%.fastq.gz
	@mkdir -p data/qc
	@mkdir -p "$(dir $@)"
	${bsub} "trim_galore --small_rna --output_dir data/trimmed --gzip  \
		--fastqc_args '-o data/qc' '$<'"
	mv '$(basename $(basename $@))_trimmed.fq.gz' '$@'
	rm 'data/qc/$(notdir $(basename $(basename $@)))_trimmed.fq_fastqc.html'
	mkdir -p data/qc/$(notdir $(basename $(basename $@)))
	unzip -p \
		'data/qc/$(notdir $(basename $(basename $@)))_trimmed.fq_fastqc.zip' \
		'$(notdir $(basename $(basename $@)))_trimmed.fq_fastqc/fastqc_data.txt' \
		| sed /^Filename/s/_trimmed// \
		> data/qc/$(notdir $(basename $(basename $@)))/fastqc_data.txt
	rm 'data/qc/$(notdir $(basename $(basename $@)))_trimmed.fq_fastqc.zip'

.PHONY: length-distribution
## Calculate read length distribution
length-distribution: ${lengths}

data/%.lengths.tsv: data/trimmed/%.fastq.gz
	gunzip -c '$<' | ./scripts/size-distribution > '$@'

.PHONY: mapped-reads
## Perform read mapping
mapped-reads: ${mapped-reads}

data/mapped/%.bam: data/trimmed/%.fastq.gz ${index} ${annotation}
	@mkdir -p "$(dir $@)"
	${bsub} -n 6 -M 24000 -R'select[mem>24000]' -R'rusage[mem=24000]' \
		"STAR --runThreadN 6 --genomeDir '$(dir ${index})' \
		--runMode alignReads --alignEndsType Local \
		--sjdbGTFfile '${annotation}' \
		--readFilesIn '$<' --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM SortedByCoordinate --outFileNamePrefix '$(basename $@)'"
	mv '$(basename $@)Aligned.sortedByCoord.out.bam' '$@'

.PHONY: rrna-contamination
## Determine rRNA contamination by assessing what fraction of mapped reads maps
## to rRNA
rrna-contamination: ${rrna-contamination}

data/rrna-contamination/%.tsv: data/mapped/%.bam ${rrna-annotation}
	@mkdir -p "$(dir $@)"
	${bsub} "featureCounts -t exon -g gene_id -M -a ${rrna-annotation} -o '$@' '$<'"

.PHONY: qc-report
## Generate aggregate report from individual tool/QC outputs
qc-report: data/qc/multiqc_report.html

data/qc/multiqc_report.html: ${trimmed-files} ${mapped-reads} ${rrna-contamination}
	multiqc --force --outdir data/qc \
		data/trimmed data/qc data/mapped data/rrna-contamination

.DEFAULT_GOAL := show-help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
	-e "H; \
		n; \
		s/^## //; \
		t doc" \
	-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) == Darwin && echo '--no-init --raw-control-chars')
