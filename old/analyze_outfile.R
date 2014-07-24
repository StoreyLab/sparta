library(ggplot2)

dat = read.csv("verbose_outfile", sep="\t")

ggplot(dat, aes(err.BY, err.RM)) + geom_hex()

dat.noNA = droplevels(dat[complete.cases(dat), ])
ggplot(dat.noNA, aes(err.BY, err.RM)) + geom_hex() + facet_wrap(~ category)

ggplot(dat, aes(prob.BY)) + geom_histogram(binwidth=.005)
