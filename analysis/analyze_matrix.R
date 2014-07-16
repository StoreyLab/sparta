library(dplyr)
library(reshape2)
library(ggplot2)
library(limma)

count.matrix = as.matrix(read.table("count_matrix.txt", header=TRUE))

count.matrix.filtered = count.matrix[rowSums(count.matrix) >= 20, ]

count.matrix.norm = voom(count.matrix.filtered, normalize.method = "scale")$E

counts = tbl_df(melt(count.matrix.norm))
colnames(counts) = c("gene", "sample", "count")

counts = counts %>% mutate(allele=gsub("(..)_.*", "\\1", sample)) %>%
                    mutate(condition=gsub(".*_(.).*", "\\1", sample)) %>%
                    mutate(replicate=gsub(".*_.(1|2).*", "\\1", sample)) %>%
                    mutate(mating=gsub(".*[12](a|alpha).*", "\\1", sample)) %>%
                    mutate(time=gsub(".*a(\\d+)$", "\\1", sample))

counts = counts %>% filter(!grepl("Unassigned", sample))

summarized = counts %>% group_by(gene) %>% summarize(count=sum(count)) %>% arrange(desc(count))

counts$time = as.numeric(counts$time)

plotgene = function(genename) {
    ingene = counts %>% filter(gene == "YDR077W")
    ggplot(ingene, aes(x=time, y=count, color=allele, lty=replicate)) + geom_line() + facet_grid(mating ~ condition)    
}

plotgene("YMR303C")

#counts %>% group_by(gene) %>% summarize(count=sum(count)) %>% arrange(desc(count))

heatmap(cor(count.matrix, method="spearman"))
