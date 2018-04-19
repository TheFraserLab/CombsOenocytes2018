suppressMessages({
  library("sleuth")
})

sample_id <- c('fbfemale_r1', 'fbfemale_r2', 'fbmale_r1', 'fbmale_r2', 'oefemale_r1', 'oefemale_r2', 'oemale_r1', 'oemale_r2')
kal_dirs <- file.path('analysis', 'sim', sample_id, 'unsplit')

s2c <- read.table("SleuthDesign.txt", header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

oe <- dplyr::filter(s2c, grepl("oe", sample))
female <- dplyr::filter(s2c, grepl("female", sample))

so <- sleuth_prep(oe, extra_bootstrap_summary = TRUE, num_cores=1)
so <- sleuth_fit(so, ~tissue, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

so2<- sleuth_prep(female, extra_bootstrap_summary = TRUE, num_cores=1)
so2<- sleuth_fit(so2, ~tissue, 'full')
so2<- sleuth_fit(so2, ~1, 'reduced')
so2<- sleuth_lrt(so2, 'reduced', 'full')


sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
sleuth_table2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant2 <- dplyr::filter(sleuth_table2, qval <= 0.05)

print(length(sleuth_significant$qval))
print(length(sleuth_significant2$qval))

write.table(sleuth_significant, 'analysis/sim/combined/oenocyte_sleuth.tsv', 
            quote=FALSE, sep="	", )
write.table(sleuth_significant2, 'analysis/sim/combined/female_sleuth.tsv', 
            quote=FALSE, sep="	", )
