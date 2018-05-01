#!/usr/bin/env Rscript

library("optparse")

suppressMessages({
  library("sleuth")
})
option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, action = "store",
              help="input  filename", metavar="character")
                   )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sample_id <- c('fbfemale_r1', 'fbfemale_r2', 'fbmale_r1', 'fbmale_r2', 'oefemale_r1', 'oefemale_r2', 'oemale_r1', 'oemale_r2')
kal_dirs <- file.path(opt$input_dir, sample_id, 'unsplit')

s2c <- read.table("SleuthDesign.txt", header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

oe <- dplyr::filter(s2c, grepl("oe", sample))
female <- dplyr::filter(s2c, grepl("female", sample))
fb <- dplyr::filter(s2c, grepl("fb", sample))
male <- dplyr::filter(s2c, grepl("(oe|fb)male", sample))

so <- sleuth_prep(oe, extra_bootstrap_summary = TRUE, num_cores=1)
so <- sleuth_fit(so, ~tissue, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

so2<- sleuth_prep(female, extra_bootstrap_summary = TRUE, num_cores=1)
so2<- sleuth_fit(so2, ~tissue, 'full')
so2<- sleuth_fit(so2, ~1, 'reduced')
so2<- sleuth_lrt(so2, 'reduced', 'full')

so3<- sleuth_prep(fb, extra_bootstrap_summary = TRUE, num_cores=1)
so3<- sleuth_fit(so3, ~tissue, 'full')
so3<- sleuth_fit(so3, ~1, 'reduced')
so3<- sleuth_lrt(so3, 'reduced', 'full')

so4<- sleuth_prep(male, extra_bootstrap_summary = TRUE, num_cores=1)
so4<- sleuth_fit(so4, ~tissue, 'full')
so4<- sleuth_fit(so4, ~1, 'reduced')
so4<- sleuth_lrt(so4, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

sleuth_table2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant2 <- dplyr::filter(sleuth_table2, qval <= 0.05)

sleuth_table3 <- sleuth_results(so3, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant3 <- dplyr::filter(sleuth_table3, qval <= 0.05)

sleuth_table4 <- sleuth_results(so4, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant4 <- dplyr::filter(sleuth_table4, qval <= 0.05)

print(length(sleuth_significant$qval))
print(length(sleuth_significant2$qval))

write.table(sleuth_significant,
            file.path(opt$input_dir, 'combined/oe_sleuth.tsv'),
            row.names=FALSE,
            quote=FALSE, sep="\t", )
write.table(sleuth_significant2,
            file.path(opt$input_dir, 'combined/female_sleuth.tsv'),
            row.names=FALSE,
            quote=FALSE, sep="\t", )
write.table(sleuth_significant3,
            file.path(opt$input_dir, 'combined/fb_sleuth.tsv'),
            row.names=FALSE,
            quote=FALSE, sep="\t", )
write.table(sleuth_significant4,
            file.path(opt$input_dir, 'combined/male_sleuth.tsv'),
            row.names=FALSE,
            quote=FALSE, sep="\t", )

write.table(so$obs_norm,
            file.path(opt$input_dir, 'combined/sleuth_oe_obs_norm.tsv'),
            row.names=FALSE,
            quote=FALSE, sep="\t")
write.table(so2$obs_norm,
            file.path(opt$input_dir, 'combined/sleuth_female_obs_norm.tsv'),
            row.names=FALSE,
            quote=FALSE, sep="\t")
write.table(so3$obs_norm,
            file.path(opt$input_dir, 'combined/sleuth_fb_obs_norm.tsv'),
            row.names=FALSE,
            quote=FALSE, sep="\t")
write.table(so4$obs_norm,
            file.path(opt$input_dir, 'combined/sleuth_male_obs_norm.tsv'),
            row.names=FALSE,
            quote=FALSE, sep="\t")
