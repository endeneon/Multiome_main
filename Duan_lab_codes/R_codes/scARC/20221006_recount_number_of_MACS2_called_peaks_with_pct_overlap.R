# Chuxuan Li 10/06/2022
# count number of peaks in the same way as 11Apr2022_count_MACS2_recalled_peaks
#by combining peaks separately in each cell type, then count the number of peaks
#that came from merging any combinations of 0,1,6hr

# init ####
library(Signac)
library(stringr)
library(GenomicRanges)
library(S4Vectors)

load("./macs2_called_new_peak_by_celltype_separated.RData")
names(peaks_uncombined) <- str_replace(names(peaks_uncombined), "unknown", "Others")

# clean peaks ####
for (i in 1:length(peaks_uncombined)) {
  peaks_uncombined[[i]] <- peaks_uncombined[[i]][str_detect(peaks_uncombined[[i]]@seqnames, "^chr"), ]
  peaks_uncombined[[i]] <- keepStandardChromosomes(peaks_uncombined[[i]], pruning.mode = "coarse")
  peaks_uncombined[[i]] <- subsetByOverlaps(x = peaks_uncombined[[i]], ranges = blacklist_hg38_unified, invert = TRUE)
}

# check the distribution of peak length ####
for (i in 1:length(peaks_uncombined)) {
  png(paste0("./MACS2_called_peak_histograms/MACS2_peak_width_distribution_", 
             names(peaks_uncombined)[i], ".png"), 
      width = 600, height = 480, pointsize = 20)
  hist(peaks_uncombined[[i]]@ranges@width, breaks = 500, 
       main = paste0("Distribution of peak width - ", 
                     str_replace_all(names(peaks_uncombined)[i], "_", " ")),
       xlab = "Peak width")
  dev.off()
}

for (i in 1:length(peaks_uncombined)) {
  cat(names(peaks_uncombined)[i], median(peaks_uncombined[[i]]@ranges@width), "\n")
}

for (i in 1:length(peaks_uncombined)){
  print(paste(i, unique(peaks_uncombined[[i]]$ident), sep = ","))
  names(peaks_uncombined)[i] <- unique(peaks_uncombined[[i]]$ident)
}

# test runs ####

r1 = GRanges(seqnames = c("chr1", "chr1", "chr1"),
             ranges = IRanges(start = c(1, 5, 10), width = 4))
r2 = GRanges(seqnames = c("chr1", "chr1"),
             ranges = IRanges(start = c(4, 12), width = 4))
o = findOverlaps(r1, r2, maxgap = -1L, # -1L allows two ranges start = end, 0L allows start = end + 1
                 minoverlap = 0)
# Problem: which number should we take as the number of peaks shared by two sets?
# if we follow the logic of reduce() with revmap = True, then basically we generated
#a union peak set for query and subject, and we counted which union peak was unique
#and which came from a combination of both sets. In this case, we can still argue
#that the region covered by the union peak was "open" at both time points, and
#with the limit of >25% median peak width, the overlapping peaks should not be
#as many as in the test.
# solution: merge peaks like reduce() but with >25% median peak width overlap,
#either figure out how to manipulate reduce() or write own function.
countOverlaps(r1, r2, minoverlap = 0L)
queryHits(o)
subjectHits(o)
reduce(Reduce(f = c, x = list(r1, r2)), min.gapwidth = -2, with.revmap = T)
# does not support negative min gapwidth

ps1 <- peaks_uncombined$GABA_0hr
ps2 <- peaks_uncombined$GABA_1hr
mbo_res <- mergeByOverlaps(ps1, ps2, minoverlap = as.integer(465*0.25))
mbo_res <- mbo_res[,c(1, 9)]
r1 = mbo_res[1, 1]
r2 = mbo_res[1, 2]
as.vector(r1@seqnames)
r1@ranges@start
mbo_res <- mergeByOverlaps(peaks_uncombined$GABA_0hr[1:5,], 
                           peaks_uncombined$GABA_1hr[1:5,], 
                           minoverlap = as.integer(465 * 0.25))

mergePeaksByMinOverlap(peaks_uncombined$GABA_0hr[1:5,], 
                       peaks_uncombined$GABA_1hr[1:5,], 465)

# main function ####
mergePeaksByMinOverlap <- function(ps1, ps2, med_peak_width) {
  mbo_res <- mergeByOverlaps(ps1, ps2, minoverlap = as.integer(med_peak_width * 0.25))
  mbo_res <- mbo_res[, c(1, 9)]
  merged_seqnames <- c()
  merged_starts <- c()
  merged_widths <- c()
  i = 1
  while (i <= nrow(mbo_res)) {
    cat("\ni: ", i)
    r1_curr <- mbo_res[i, 1]
    r2_curr <- mbo_res[i, 2]
    r1_first <- r1_curr
    r2_first <- r2_curr
    cat("\nr1_first:", r1_first@ranges@start, r1_first@ranges@width)
    cat("\nr2_first:", r2_first@ranges@start, r2_first@ranges@width)
    samer1 <- 0
    samer2 <- 0 
    second_to_last_same_as_last <- F
    if (i < nrow(mbo_res)) {
      r1_next <- mbo_res[i + 1, 1]
      r2_next <- mbo_res[i + 1, 2]
      samer1 <- countOverlaps(granges(r1_curr), granges(r1_next), type = "equal")
      samer2 <- countOverlaps(granges(r2_curr), granges(r2_next), type = "equal")
      
      j <- i
      while (samer1 | samer2 & j < (nrow(mbo_res) - 1)) {
        j <- j + 1
        cat("\nj: ", j)
        r1_curr <- mbo_res[j, 1]
        r2_curr <- mbo_res[j, 2]
        r1_next <- mbo_res[j + 1, 1]
        r2_next <- mbo_res[j + 1, 2]
        samer1 <- countOverlaps(granges(r1_curr), granges(r1_next), type = "equal")
        samer2 <- countOverlaps(granges(r2_curr), granges(r2_next), type = "equal")
      }
    }
    if (samer1 | samer2 & j == (nrow(mbo_res) - 1)) {
      second_to_last_same_as_last <- T
    }
    if (samer1 | samer2) {
      r1_final <- r1_next
      r2_final <- r2_next
    } else {
      r1_final <- r1_curr
      r2_final <- r2_curr
    }
    cat("\nr1_final:", r1_final@ranges@start, r1_final@ranges@width)
    cat("\nr2_final:", r2_final@ranges@start, r2_final@ranges@width)
    firstr <- chooseFirstRange(r1_first, r2_first)
    finalr <- chooseLastRange(r1_final, r2_final)
    merged_row_stats <- mergeTwoRanges(firstr, finalr)
    merged_seqnames <- c(merged_seqnames, merged_row_stats$seqname)
    merged_starts <- c(merged_starts, merged_row_stats$start)
    merged_widths <- c(merged_widths, merged_row_stats$width)
    
    if (j > i) {
      i <- j + 1
    } else {
      i <- i + 1
    }
    if (i == (nrow(mbo_res) - 1) & second_to_last_same_as_last) {
      break
    }
  }
  return(GRanges(seqnames = merged_seqnames, 
                 ranges = IRanges(start = merged_starts, width = merged_widths)))
}

# auxiliary function ####
chooseFirstRange <- function(r1, r2) {
  # Given two GRanges objects, select the one that has the smaller start value
  # Input:
  #   r1, r2: GenomicRanges objects with seqnames, start, and width at least
  # Returns:
  #   the range that has the smaller "start" value 
  
  # check if two ranges are from the same chromosome
  seq1 <- as.vector(r1@seqnames)
  seq2 <- as.vector(r2@seqnames)
  stopifnot(exprs = seq1 == seq2)
  # get start, width from the two ranges, pick the widest range formable by them
  start1 <- r1@ranges@start
  start2 <- r2@ranges@start
  if (start1 < start2) {
    return(r1) # r1 should be used to produce the start in downstream if it starts earlier
  } else {
    return(r2) # r2 starts earlier, or r1 r2 has the same starts
  }
}
chooseLastRange <- function(r1, r2) {
  # Given two GRanges objects, select the one that has the greater end value
  # Input:
  #   r1, r2: GenomicRanges objects with seqnames, start, and width at least
  # Returns:
  #   the range that has the greater "end" value 
  
  # check if two ranges are from the same chromosome
  seq1 <- as.vector(r1@seqnames)
  seq2 <- as.vector(r2@seqnames)
  stopifnot(exprs = seq1 == seq2)
  end1 <- r1@ranges@start + r1@ranges@width - 1
  end2 <- r2@ranges@start + r2@ranges@width - 1
  if (end1 > end2) {
    return(r1) # r1 should be used to produce the end in downstream if it ends later
  } else {
    return(r2) # r2 ends later, or r1 r2 has the same ends
  }
}
mergeTwoRanges <- function(r1, r2) {
  # Given two GRanges objects, merge by the longest width possibly formed by
  #these two ranges
  # Input:
  #   r1, r2: GenomicRanges objects with seqnames, start, and width at least
  # Returns:
  #   values corresponding to the seqname, start, and width of the resulting range
  
  # check if two ranges are from the same chromosome
  seq1 <- as.vector(r1@seqnames)
  seq2 <- as.vector(r2@seqnames)
  stopifnot(exprs = seq1 == seq2)
  # get start, width from the two ranges, pick the widest range formable by them
  start1 <- r1@ranges@start
  start2 <- r2@ranges@start
  start.out <- min(start1, start2)
  end1 <- start1 + r1@ranges@width - 1
  end2 <- start2 + r2@ranges@width - 1
  end.out <- max(end1, end2)
  width.out <- end.out - start.out + 1
  return(data.frame(seqname = seq1, 
                    start = start.out, 
                    width = width.out))
}


# run it####
GABA_overlapping_10 <- mergePeaksByMinOverlap(peaks_uncombined$GABA_0hr, peaks_uncombined$GABA_1hr,
                       465)
GABA_overlapping_60 <- mergePeaksByMinOverlap(peaks_uncombined$GABA_0hr, peaks_uncombined$GABA_6hr,
                                              465)
GABA_overlapping_16 <- mergePeaksByMinOverlap(peaks_uncombined$GABA_1hr, peaks_uncombined$GABA_6hr,
                                              465)
GABA_overlapping_016 <- mergePeaksByMinOverlap(peaks_uncombined$GABA_0hr, GABA_overlapping_16,
                                               465)
GABA_overlapping_610 <- mergePeaksByMinOverlap(peaks_uncombined$GABA_6hr, GABA_overlapping_10,
                                               465)
GABA_overlapping_160 <- mergePeaksByMinOverlap(peaks_uncombined$GABA_1hr, GABA_overlapping_60,
                                               465)
save(GABA_overlapping_016, GABA_overlapping_610, GABA_overlapping_160, file = "GABA_overlapped_ranges_3way.RData")
nmglut_overlapping_10 <- mergePeaksByMinOverlap(peaks_uncombined$NEFM_neg_glut_0hr, 
                                                peaks_uncombined$NEFM_neg_glut_1hr,
                                              446)
nmglut_overlapping_60 <- mergePeaksByMinOverlap(peaks_uncombined$NEFM_neg_glut_0hr, 
                                                peaks_uncombined$NEFM_neg_glut_6hr,
                                              446)
nmglut_overlapping_16 <- mergePeaksByMinOverlap(peaks_uncombined$NEFM_neg_glut_1hr, 
                                                peaks_uncombined$NEFM_neg_glut_6hr,
                                              446)
nmglut_overlapping_016 <- mergePeaksByMinOverlap(peaks_uncombined$NEFM_neg_glut_0hr, 
                                                 nmglut_overlapping_16,
                                               446)

npglut_overlapping_10 <- mergePeaksByMinOverlap(peaks_uncombined$NEFM_pos_glut_0hr, 
                                                peaks_uncombined$NEFM_pos_glut_1hr,
                                                468)
npglut_overlapping_60 <- mergePeaksByMinOverlap(peaks_uncombined$NEFM_pos_glut_0hr, 
                                                peaks_uncombined$NEFM_pos_glut_6hr,
                                                468)
npglut_overlapping_16 <- mergePeaksByMinOverlap(peaks_uncombined$NEFM_pos_glut_1hr, 
                                                peaks_uncombined$NEFM_pos_glut_6hr,
                                                468)
npglut_overlapping_016 <- mergePeaksByMinOverlap(peaks_uncombined$NEFM_pos_glut_0hr, 
                                                 npglut_overlapping_16,
                                                 468)
save(GABA_overlapping_10, GABA_overlapping_16, GABA_overlapping_60,
     file = "GABA_overlapped_ranges_merged.RData")
save(nmglut_overlapping_10, nmglut_overlapping_16, nmglut_overlapping_60,
     file = "nmglut_overlapped_ranges_merged.RData")
save(npglut_overlapping_10, #npglut_overlapping_16, npglut_overlapping_60,
     file = "npglut_overlapped_ranges_merged_10.RData")
save(npglut_overlapping_60, file = "npglut_overlapped_ranges_merged_60.RData")
save(npglut_overlapping_16, file = "npglut_overlapped_ranges_merged_16.RData")
