library(dplyr)
library(tidyr)
library(ggplot2)

annotation = tbl_df(read.table('data/reference/Apis_mellifera.GCA_000002195.1.31.gtf', sep ='\t', comment = '#'))
annotation = filter(annotation, V3 == 'gene')
annotation = extract(annotation, V9, c('gene', 'type'), 'gene_id ([^;]+);.*gene_biotype ([^;]+)') %>% select(gene, type)

counts = tbl_df(read.delim('data/all-quant.tsv'))
colnames(counts) = sub('_R._001.bam', '', sub('data.mapped.', '', colnames(counts)))
counts = inner_join(annotation, counts, by = c(gene = 'Geneid'))
count_summary$all = rowSums(select(count_summary, -type))

count_summary_long = gather(count_summary %>% top_n(5, all), library, count, -type, -all)
p = ggplot(count_summary_long) + aes(type, count) + geom_bar(stat = 'identity') + scale_y_log10() + facet_wrap(~library)
preview(plot(p + theme(axis.text.x = element_text(angle = 90, hjust = 1))))
