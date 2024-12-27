# Chuxuan Li 05/01/2023
# count number of peaks shared by 2 or 3 time points for 018-029 combined data
#to plot Venn diagram

# init ####
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(stringr)
library(S4Vectors)

load("./018-029_macs2_called_new_peaks_sep_by_timextype.RData")

# define aux. functions ####
# from signac peaks.R for combining peaks 
combinePeaks <- function(peakList){
  gr.combined = Reduce(f = c, x = peakList) # does not reduce length
  print(unique(gr.combined$ident))
  gr = reduce(x = gr.combined, with.revmap = TRUE)
  dset.vec = vector(mode = "character", length = length(x = gr))
  ident.vec = gr.combined$ident
  revmap = gr$revmap # here gr is already reduced length
  for (i in seq_len(length.out = length(x = gr))) {
    datasets = ident.vec[revmap[[i]]]
    dset.vec[[i]] = paste(unique(x = datasets), collapse = ",")
  }
  gr$peak_called_in = dset.vec # dset.vec length has to be the same with gr, so gr length is already reduced at this step
  #gr$revmap = NULL
  return(gr)
}

# # clean peaks ####
# for (i in 1:length(peaks_uncombined)) {
#   peaks_uncombined[[i]] <- peaks_uncombined[[i]][str_detect(peaks_uncombined[[i]]@seqnames, "^chr"), ]
# }

# combine all peaks ####
gr <- combinePeaks(peaks_uncombined)
gr <- gr[str_detect(gr@seqnames, "^chr"), ]
gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
gr <- subsetByOverlaps(x = gr, ranges = blacklist_hg38_unified, invert = TRUE)
save(gr, file = "all_peaks_combined_manually_from_peaks_uncombined_object.RData")
unique(gr$peak_called_in)

# count total peaks in each subset ####
txt_counts <- matrix(nrow = 3, ncol = 4, dimnames = list(c("0hr", "1hr", "6hr"), 
                                                         c("GABA", "npglut", "nmglut", "others")))
for (i in 1:length(peaks_uncombined)) {
  tt <- as.character(unique(peaks_uncombined[[i]]$ident))
  print(tt)
  if (i == 1) {
    sorted_tt <- c(tt)
  } else {
    sorted_tt <- c(sorted_tt, tt)
  }
}
sorted_tt <- order(sorted_tt)
peaks_uncombined <- peaks_uncombined[sorted_tt]
for (i in 1:length(peaks_uncombined)) {
  tt <- as.character(unique(peaks_uncombined[[i]]$ident))
  print(tt)
  if (i == 1) {
    sorted_tt <- c(tt)
  } else {
    sorted_tt <- c(sorted_tt, tt)
  }
  txt_counts[i] <- sum(str_detect(gr$peak_called_in, tt))
}
write.table(txt_counts, 
            file = "num_peaks_in_each_typextime_sep_called_peaksets_after_filtering_chr_overlaps.csv",
            quote = F, sep = ",")

# combine peaks loses some peaks for some reason
testr <- combinePeaks(list(peaks_uncombined$GABA_0hr, peaks_uncombined$GABA_1hr))
testr <- testr[str_detect(testr@seqnames, "^chr"), ]
testr <- keepStandardChromosomes(testr, pruning.mode = "coarse")
testr <- subsetByOverlaps(x = testr, ranges = blacklist_hg38_unified, invert = TRUE)
sum(str_detect(testr$peak_called_in, "GABA_0hr")) #163415 -> combining cause some peaks to be lost
sum(str_detect(testr$peak_called_in, "GABA_1hr")) #199655
test <- peaks_uncombined$GABA_0hr[str_detect(peaks_uncombined$GABA_0hr@seqnames, "^chr"), ]
test <- keepStandardChromosomes(test, pruning.mode = "coarse")
test <- subsetByOverlaps(x = test, ranges = blacklist_hg38_unified, invert = TRUE) #164433
reduce(test) #164433
test <- peaks_uncombined$GABA_1hr[str_detect(peaks_uncombined$GABA_1hr@seqnames, "^chr"), ]
test <- keepStandardChromosomes(test, pruning.mode = "coarse")
test <- subsetByOverlaps(x = test, ranges = blacklist_hg38_unified, invert = TRUE) #203916
reduce(test) #203916

test <- Reduce(f = c, x = list(peaks_uncombined$GABA_0hr, peaks_uncombined$GABA_1hr)) #399137
test$ident[revmap[[1]]]
reduce(test, with.revmap = TRUE) # 231823
test <- test[str_detect(test@seqnames, "^chr"), ]
test <- keepStandardChromosomes(test, pruning.mode = "coarse")
test <- subsetByOverlaps(x = test, ranges = blacklist_hg38_unified, invert = TRUE) #368349
203916+164433 #368349 -> the Biogenerics/base Reduce does not change number of peaks 

# count intersections ####
load("./018-029_macs2_called_new_peaks_combined.RData")
peaks <- peaks[str_detect(peaks@seqnames, "^chr"), ]
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

#1. peaks in all cell types
times <- c("0hr", "1hr", "6hr")
for (i in 1:3) {
  timespec_names <- sorted_tt[str_detect(sorted_tt, times[i])]
  print(timespec_names)
  print(sum(str_detect(peaks$peak_called_in, timespec_names[1]) & 
              str_detect(peaks$peak_called_in, timespec_names[2]) & 
              str_detect(peaks$peak_called_in, timespec_names[3]) & 
              str_detect(peaks$peak_called_in, timespec_names[4])))
}

# 2. peaks in only 1 cell type
for (i in 1:3) {
  timespec_names <- sorted_tt[str_detect(sorted_tt, times[i])]
  for (j in 1:length(timespec_names)) {
    print(timespec_names[j])
    print(sum(peaks$peak_called_in == timespec_names[j]))
  }
}

# 3. peaks shared by 2 cell types
for (i in 1:3) {
  cat("\n\n", times[i])
  timespec_names <- sorted_tt[str_detect(sorted_tt, times[i])]
  combs <- combn(timespec_names, 2, simplify = F)
  for (j in 1:length(combs)) {
    print(combs[[j]])
    excluded.idents <- timespec_names[!timespec_names %in% combs[[j]]]
    cat("does not contain: ", excluded.idents)
    print(sum(str_detect(peaks$peak_called_in, combs[[j]][1]) & 
                str_detect(peaks$peak_called_in, combs[[j]][2]) & 
                str_detect(peaks$peak_called_in, excluded.idents[1], negate = T) & 
                str_detect(peaks$peak_called_in, excluded.idents[2], negate = T)))
  }
}

# 4. peaks shared by 3 cell types
library(combinat)
for (i in 1:3) {
  cat("\n\n", times[i])
  timespec_names <- sorted_tt[str_detect(sorted_tt, times[i])]
  combs <- combn(timespec_names, 3, simplify = F)
  for (j in 1:length(combs)) {
    print(combs[[j]])
    excluded.ident <- timespec_names[!timespec_names %in% combs[[j]]]
    cat("does not contain: ", excluded.ident)
    print(sum(str_detect(peaks$peak_called_in, combs[[j]][1]) & 
                str_detect(peaks$peak_called_in, combs[[j]][2]) & 
                str_detect(peaks$peak_called_in, combs[[j]][3]) & 
                str_detect(peaks$peak_called_in, excluded.ident, negate = T)))
  }
}
