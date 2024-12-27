# quantify VPS45 at rs2027349 site
# quantify multiple regions
# Siwei 20 Dec 2019

# init
library(stats)

#
VPS45_isogenic_counts_samples <- VPS45_isogenic_counts_all_45M[, c(7:14)]

# VPS45_read <- data.frame(c(119, 106, 122, 113, 117, 107, 130, 50.7))
VPS45_read <- data.frame(base_1 = 1:ncol(VPS45_isogenic_counts_samples))
VPS45_read$base_1 <- unlist(VPS45_isogenic_counts_samples[1, ])
VPS45_read$base_7 <- unlist(VPS45_isogenic_counts_samples[2, ])
VPS45_read$base_21 <- unlist(VPS45_isogenic_counts_samples[3, ])
# colnames(VPS45_read)[1] <- "total_reads"
# VPS45_read$normalised_reads <- VPS45_read$raw_reads/VPS45_read$total_reads*115.3
rownames(VPS45_read) <- colnames(VPS45_isogenic_counts_samples)

VPS45_read$normalised_reads <- VPS45_read$total_reads
t.test(x = VPS45_read$base_1[1:3],
       y = VPS45_read$base_1[6:7],
       alternative = "l",
       var.equal = T)

t.test(x = VPS45_read$normalised_reads[1:2],
       y = VPS45_read$normalised_reads[6:7],
       alternative = "l",
       var.equal = T)

t.test(x = VPS45_read$normalised_reads[4:5],
       y = VPS45_read$normalised_reads[6:7],
       alternative = "l",
       var.equal = T)
t.test(x = VPS45_read$normalised_reads[1:3],
       y = VPS45_read$normalised_reads[4:5],
       alternative = "l",
       var.equal = T)


###

plot_VPS45_isogenic(chr = "chr1",
                    start = 150066791,
                    end = 150068025, ylimit = 600)
