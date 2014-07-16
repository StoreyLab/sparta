error.rates = NULL

for (folder in list.files("../all_sample_analysis/reports/", full.names=TRUE)) {
    infile = paste0(folder, "/mismatch_prob_info.txt")
    d = read.table(infile, col.names = c("quality", "errorrate", "number"))
    sample = gsub(".*reports\\/\\/(.*)\\/.*", "\\1", infile)
    d$sample = sample
    error.rates = rbind(error.rates, d)
}

error.rates = error.rates %>% mutate(condition=gsub("^(.).*", "\\1", sample)) %>%
    mutate(replicate=gsub(".(1|2).*", "\\1", sample)) %>%
    mutate(mating=gsub(".*[12](a|alpha).*", "\\1", sample)) %>%
    mutate(time=gsub(".*a(\\d+)$", "\\1", sample))

error.rates = error.rates %>% mutate(expected=10^(-(quality-33)/10))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(error.rates, aes(expected, errorrate, col=condition)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_abline(col="black") + scale_color_manual(values=cbPalette)

