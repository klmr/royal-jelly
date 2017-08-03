include dirs.make
include help.make
include bsub.make

$(call dirs,data data/mirna data/mirna/reference data/mirna/reference/index data/mirna/mapped data/mirna/bee-reads data/mirna/summary)

hairpin = data/mirna/reference/hairpin.fasta
hairpin-index = $(dir ${hairpin})index
raw-mapped = $(shell find data/mapped/long -name \*.bam)
reads = $(subst /mapped/long/,/mirna/bee-reads/,${raw-mapped:.bam=.fastq.gz})
mapped = $(subst /bee-reads/,/mapped/,${reads:.fastq.gz=.Aligned.sortedByCoord.out.bam})
summary = $(subst /mapped/,/summary/,${mapped:.Aligned.sortedByCoord.out.bam=.txt})

## Download hairpin miRNA reference.
hairpin: ${hairpin}
.PHONY: hairpin

## Build the STAR index for the hairpin miRNA reference.
hairpin-index: ${hairpin-index}/Genome
.PHONY: hairpin-index

## Generate the summary of found miRNAs.
summary: ${summary}
.PHONY: summary

${hairpin}: | data/mirna/reference
	wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz -O $@.gz
	gunzip $@.gz

${hairpin-index}/Genome: ${hairpin} | ${hairpin-index}
	${bsub} -n3 $(call memreq,16000) \
		"STAR --runMode genomeGenerate \
		--genomeFastaFiles $< \
		--genomeDir ${hairpin-index}"

data/mirna/bee-reads/%.fastq.gz: data/mapped/long/%.bam | data/mirna/bee-reads
	REGIONS=$$(cut -f1 ${apis-annotation} | sort -u | grep -v '^#'); \
	samtools view -hbF4 $< $$REGIONS | bedtools bamtofastq -i - -fq /dev/stdout | gzip -c > $@

data/mirna/mapped/%.Aligned.sortedByCoord.out.bam: data/mirna/bee-reads/%.fastq.gz ${hairpin-index} | data/mirna/mapped
	${bsub} -n8 $(call memreq,16000) \
		"STAR --runThreadN 6 --genomeDir $(lastword $^) \
		--outFilterMismatchNoverLmax 0.15 --outFilterMismatchNmax 1 \
		--alignIntronMax 1 \
		--scoreDelOpen -10000 --scoreInsOpen -10000 \
		--outFilterMultimapNmax 100 \
		--readFilesIn $< --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix ${@:Aligned.sortedByCoord.out.bam=}"

data/mirna/summary/%.txt: data/mirna/mapped/%.Aligned.sortedByCoord.out.bam | data/mirna/summary
	samtools view $< | cut -f 3 | sort | uniq -c | sort -rnk1 > $@

# vim: ft=make
