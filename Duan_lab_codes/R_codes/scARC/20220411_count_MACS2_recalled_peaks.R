# Chuxuan Li 04/11/2022
# count number of peaks in each intersection between cell types for 18-line 
#10x-aggregated data 

# init ####
library(Signac)
library(stringr)
library(S4Vectors)

load("./macs2_called_new_peak_by_celltype_separated.RData")

# count number of peaks unique/intersecting for each cell type by time point
for (i in 1:length(peaks_uncombined)){
  print(paste(i, unique(peaks_uncombined[[i]]$ident), sep = ","))
  names(peaks_uncombined)[i] <- unique(peaks_uncombined[[i]]$ident)
}

# define aux. functions ####
# from signac peaks.R for combining peaks 
combinePeaks <- function(peakList){
  gr.combined = Reduce(f = c, x = peakList)
  print(unique(gr.combined$ident))
  gr = reduce(x = gr.combined, with.revmap = TRUE)
  dset.vec = vector(mode = "character", length = length(x = gr))
  ident.vec = gr.combined$ident
  revmap = gr$revmap
  for (i in seq_len(length.out = length(x = gr))) {
    datasets = ident.vec[revmap[[i]]]
    dset.vec[[i]] = paste(unique(x = datasets), collapse = ",")
  }
  gr$peak_called_in = dset.vec
  gr$revmap = NULL
  return(gr)
}

# # function to generate all 15 possible combinations of idents
# generateAllCombs <- function(peakList){
#   ones <- c()
#   for (l in peakList){
#     ones <- c(ones, unique(l$ident))
#   }
#   print(ones)
#   fours <- 
#   for (i in 1:length(ones)){
#     
#   }
#   four <- paste(ones, sep = ",", collapse = ",")
#   print(four)
#   combs <- c(ones, four)
#   for (i in 1:length(ones)){
#     j = 1
#     while ((i + j) <= length(ones)){
#       two <- paste(ones[i], ones[i + j], sep = ",")
#       print(two)
#       combs <- c(combs, two)
#       j = j + 1
#     }
#     three <- str_remove(four, paste0(ones[i], ","))
#     if (i == 4){
#       three <- str_remove(four, paste0(",", ones[i]))
#     }
#     print(three)
#     combs <- c(combs, three)
#   }
#   return(combs)
# }

generateAllPermut <- function(time){
  idents <- paste(c("GABA", "NEFM_neg_glut", "NEFM_pos_glut", "unknown"),
                  time, sep = "_")
  print(idents)
  combs <- vector("list", 4L)
  combs[[1]] <- idents
  twos <- expand.grid(idents, idents)
  twos <- twos[twos$Var1 != twos$Var2, ]
  combs[[2]] <- do.call(what = paste, args = c(twos, list(sep = ",")))
  threes <- expand.grid(idents, idents, idents)
  threes <- threes[threes$Var1 != threes$Var2 & threes$Var1 != threes$Var3 & threes$Var2 != threes$Var3, ]
  combs[[3]] <- do.call(what = paste, args = c(threes, list(sep = ",")))
  fours <- expand.grid(idents, idents, idents, idents)
  fours <- fours[fours$Var1 != fours$Var2 & 
                   fours$Var1 != fours$Var3 & 
                   fours$Var1 != fours$Var4 &
                   fours$Var2 != fours$Var3 &
                   fours$Var2 != fours$Var4 &
                   fours$Var3 != fours$Var4, ]
  combs[[4]] <- do.call(what = paste, args = c(fours, list(sep = ",")))
  return(combs)
}

# 0hr ####
peaks_0hr <- list(peaks_uncombined[[1]], peaks_uncombined[[3]], #GABA, nmglut
                  peaks_uncombined[[10]], peaks_uncombined[[12]]) #npglut, unk
for(i in peaks_0hr){
  print(unique(i$ident))
}
idents <- c("GABA_0hr", "NEFM_neg_glut_0hr", "NEFM_pos_glut_0hr", "unknown_0hr")
gr <- combinePeaks(peaks_0hr)
unique(gr$peak_called_in)

# count peaks in only 1 cell type
for (i in idents){
  print(i)
  print(sum(gr$peak_called_in == i))
}

combs <- generateAllPermut("0hr")

# count peaks intersecting 4 cell types first
gr[gr$peak_called_in %in% combs[[4]], ] 
#GRanges object with 104734 ranges and 1 metadata column

# count peaks intersecting 2 cell types next
gr[gr$peak_called_in %in% c("NEFM_neg_glut_0hr,GABA_0hr",
                            "GABA_0hr,NEFM_neg_glut_0hr")] #6840
gr[gr$peak_called_in %in% c("NEFM_pos_glut_0hr,GABA_0hr",
                            "GABA_0hr,NEFM_pos_glut_0hr")] #7809
gr[gr$peak_called_in %in% c("unknown_0hr,GABA_0hr",
                            "GABA_0hr,unknown_0hr")] #2073
gr[gr$peak_called_in %in% c("NEFM_neg_glut_0hr,NEFM_pos_glut_0hr",
                            "NEFM_pos_glut_0hr,NEFM_neg_glut_0hr")] #14034
gr[gr$peak_called_in %in% c("NEFM_pos_glut_0hr,unknown_0hr",
                            "unknown_0hr,NEFM_pos_glut_0hr")] #9011
gr[gr$peak_called_in %in% c("unknown_0hr,NEFM_neg_glut_0hr",
                            "NEFM_neg_glut_0hr,unknown_0hr")] #1063

# count peaks intersecting 3 cell types next
combs[[3]]
gr[gr$peak_called_in %in% c("NEFM_pos_glut_0hr,NEFM_neg_glut_0hr,GABA_0hr",
                            "NEFM_neg_glut_0hr,NEFM_pos_glut_0hr,GABA_0hr",
                            "GABA_0hr,NEFM_neg_glut_0hr,NEFM_pos_glut_0hr",
                            "GABA_0hr,NEFM_pos_glut_0hr,NEFM_neg_glut_0hr",
                            "NEFM_pos_glut_0hr,GABA_0hr,NEFM_neg_glut_0hr",
                            "NEFM_neg_glut_0hr,GABA_0hr,NEFM_pos_glut_0hr")] #11583
gr[gr$peak_called_in %in% c("unknown_0hr,NEFM_neg_glut_0hr,GABA_0hr",
                            "NEFM_neg_glut_0hr,unknown_0hr,GABA_0hr",
                            "GABA_0hr,NEFM_neg_glut_0hr,unknown_0hr",
                            "GABA_0hr,unknown_0hr,NEFM_neg_glut_0hr",
                            "unknown_0hr,GABA_0hr,NEFM_neg_glut_0hr",
                            "NEFM_neg_glut_0hr,GABA_0hr,unknown_0hr")] #2392
gr[gr$peak_called_in %in% c("unknown_0hr,NEFM_pos_glut_0hr,GABA_0hr",
                            "NEFM_pos_glut_0hr,unknown_0hr,GABA_0hr",
                            "GABA_0hr,NEFM_pos_glut_0hr,unknown_0hr",
                            "GABA_0hr,unknown_0hr,NEFM_pos_glut_0hr",
                            "unknown_0hr,GABA_0hr,NEFM_pos_glut_0hr",
                            "NEFM_pos_glut_0hr,GABA_0hr,unknown_0hr")] #4836
gr[gr$peak_called_in %in% c("unknown_0hr,NEFM_pos_glut_0hr,NEFM_neg_glut_0hr",
                            "NEFM_pos_glut_0hr,unknown_0hr,NEFM_neg_glut_0hr",
                            "NEFM_neg_glut_0hr,NEFM_pos_glut_0hr,unknown_0hr",
                            "NEFM_neg_glut_0hr,unknown_0hr,NEFM_pos_glut_0hr",
                            "unknown_0hr,NEFM_neg_glut_0hr,NEFM_pos_glut_0hr",
                            "NEFM_pos_glut_0hr,NEFM_neg_glut_0hr,unknown_0hr")] #4836
# combs <- generateAllCombs(peaks_0hr)
# combs
# ones <- combs[1:4]
# for (i in ones){
#   print(i)
#   print(sum(gr$peak_called_in == i))
# }

# 1hr ####
peaks_1hr <- list(peaks_uncombined[[5]], peaks_uncombined[[4]], 
                  peaks_uncombined[[6]], peaks_uncombined[[11]])
for(i in peaks_1hr){
  print(unique(i$ident))
}
idents <- c("GABA_1hr", "NEFM_neg_glut_1hr", "NEFM_pos_glut_1hr", "unknown_1hr")
gr <- combinePeaks(peaks_1hr)
unique(gr$peak_called_in)

# count peaks in only 1 cell type
for (i in idents){
  print(i)
  print(sum(gr$peak_called_in == i))
}

combs <- generateAllPermut("1hr")

# count peaks intersecting 4 cell types first
gr[gr$peak_called_in %in% combs[[4]], ] 
#GRanges object with 104734 ranges and 1 metadata column

# count peaks intersecting 2 cell types next
gr[gr$peak_called_in %in% c("NEFM_neg_glut_1hr,GABA_1hr",
                            "GABA_1hr,NEFM_neg_glut_1hr")] #8245
gr[gr$peak_called_in %in% c("NEFM_pos_glut_1hr,GABA_1hr",
                            "GABA_1hr,NEFM_pos_glut_1hr")] #7456
gr[gr$peak_called_in %in% c("unknown_1hr,GABA_1hr",
                            "GABA_1hr,unknown_1hr")] #3477
gr[gr$peak_called_in %in% c("NEFM_neg_glut_1hr,NEFM_pos_glut_1hr",
                            "NEFM_pos_glut_1hr,NEFM_neg_glut_1hr")] #16212
gr[gr$peak_called_in %in% c("NEFM_pos_glut_1hr,unknown_1hr",
                            "unknown_1hr,NEFM_pos_glut_1hr")] #11400
gr[gr$peak_called_in %in% c("unknown_1hr,NEFM_neg_glut_1hr",
                            "NEFM_neg_glut_1hr,unknown_1hr")] #1695

# count peaks intersecting 3 cell types next
combs[[3]]
gr[gr$peak_called_in %in% c("NEFM_pos_glut_1hr,NEFM_neg_glut_1hr,GABA_1hr",
                            "NEFM_neg_glut_1hr,NEFM_pos_glut_1hr,GABA_1hr",
                            "GABA_1hr,NEFM_neg_glut_1hr,NEFM_pos_glut_1hr",
                            "GABA_1hr,NEFM_pos_glut_1hr,NEFM_neg_glut_1hr",
                            "NEFM_pos_glut_1hr,GABA_1hr,NEFM_neg_glut_1hr",
                            "NEFM_neg_glut_1hr,GABA_1hr,NEFM_pos_glut_1hr")] #9883
gr[gr$peak_called_in %in% c("unknown_1hr,NEFM_neg_glut_1hr,GABA_1hr",
                            "NEFM_neg_glut_1hr,unknown_1hr,GABA_1hr",
                            "GABA_1hr,NEFM_neg_glut_1hr,unknown_1hr",
                            "GABA_1hr,unknown_1hr,NEFM_neg_glut_1hr",
                            "unknown_1hr,GABA_1hr,NEFM_neg_glut_1hr",
                            "NEFM_neg_glut_1hr,GABA_1hr,unknown_1hr")] #5026
gr[gr$peak_called_in %in% c("unknown_1hr,NEFM_pos_glut_1hr,GABA_1hr",
                            "NEFM_pos_glut_1hr,unknown_1hr,GABA_1hr",
                            "GABA_1hr,NEFM_pos_glut_1hr,unknown_1hr",
                            "GABA_1hr,unknown_1hr,NEFM_pos_glut_1hr",
                            "unknown_1hr,GABA_1hr,NEFM_pos_glut_1hr",
                            "NEFM_pos_glut_1hr,GABA_1hr,unknown_1hr")] #6536
gr[gr$peak_called_in %in% c("unknown_1hr,NEFM_pos_glut_1hr,NEFM_neg_glut_1hr",
                            "NEFM_pos_glut_1hr,unknown_1hr,NEFM_neg_glut_1hr",
                            "NEFM_neg_glut_1hr,NEFM_pos_glut_1hr,unknown_1hr",
                            "NEFM_neg_glut_1hr,unknown_1hr,NEFM_pos_glut_1hr",
                            "unknown_1hr,NEFM_neg_glut_1hr,NEFM_pos_glut_1hr",
                            "NEFM_pos_glut_1hr,NEFM_neg_glut_1hr,unknown_1hr")] #22960
# 6hr ####
peaks_6hr <- list(peaks_uncombined[[2]], peaks_uncombined[[7]], 
                  peaks_uncombined[[8]], peaks_uncombined[[9]])
for(i in peaks_6hr){
  print(unique(i$ident))
}
idents <- c("GABA_6hr", "NEFM_neg_glut_6hr", "NEFM_pos_glut_6hr", "unknown_6hr")
gr <- combinePeaks(peaks_6hr)
unique(gr$peak_called_in)

# count peaks in only 1 cell type
for (i in idents){
  print(i)
  print(sum(gr$peak_called_in == i))
}

combs <- generateAllPermut("6hr")

# count peaks intersecting 4 cell types first
gr[gr$peak_called_in %in% combs[[4]], ] 
#GRanges object with 104734 ranges and 1 metadata column

# count peaks intersecting 2 cell types next
gr[gr$peak_called_in %in% c("NEFM_neg_glut_6hr,GABA_6hr",
                            "GABA_6hr,NEFM_neg_glut_6hr")] #9133
gr[gr$peak_called_in %in% c("NEFM_pos_glut_6hr,GABA_6hr",
                            "GABA_6hr,NEFM_pos_glut_6hr")] #8149
gr[gr$peak_called_in %in% c("unknown_6hr,GABA_6hr",
                            "GABA_6hr,unknown_6hr")] #3075
gr[gr$peak_called_in %in% c("NEFM_neg_glut_6hr,NEFM_pos_glut_6hr",
                            "NEFM_pos_glut_6hr,NEFM_neg_glut_6hr")] #16204
gr[gr$peak_called_in %in% c("NEFM_pos_glut_6hr,unknown_6hr",
                            "unknown_6hr,NEFM_pos_glut_6hr")] #9317
gr[gr$peak_called_in %in% c("unknown_6hr,NEFM_neg_glut_6hr",
                            "NEFM_neg_glut_6hr,unknown_6hr")] #1596

# count peaks intersecting 3 cell types next
combs[[3]]
gr[gr$peak_called_in %in% c("NEFM_pos_glut_6hr,NEFM_neg_glut_6hr,GABA_6hr",
                            "NEFM_neg_glut_6hr,NEFM_pos_glut_6hr,GABA_6hr",
                            "GABA_6hr,NEFM_neg_glut_6hr,NEFM_pos_glut_6hr",
                            "GABA_6hr,NEFM_pos_glut_6hr,NEFM_neg_glut_6hr",
                            "NEFM_pos_glut_6hr,GABA_6hr,NEFM_neg_glut_6hr",
                            "NEFM_neg_glut_6hr,GABA_6hr,NEFM_pos_glut_6hr")] #11936
gr[gr$peak_called_in %in% c("unknown_6hr,NEFM_neg_glut_6hr,GABA_6hr",
                            "NEFM_neg_glut_6hr,unknown_6hr,GABA_6hr",
                            "GABA_6hr,NEFM_neg_glut_6hr,unknown_6hr",
                            "GABA_6hr,unknown_6hr,NEFM_neg_glut_6hr",
                            "unknown_6hr,GABA_6hr,NEFM_neg_glut_6hr",
                            "NEFM_neg_glut_6hr,GABA_6hr,unknown_6hr")] #3896
gr[gr$peak_called_in %in% c("unknown_6hr,NEFM_pos_glut_6hr,GABA_6hr",
                            "NEFM_pos_glut_6hr,unknown_6hr,GABA_6hr",
                            "GABA_6hr,NEFM_pos_glut_6hr,unknown_6hr",
                            "GABA_6hr,unknown_6hr,NEFM_pos_glut_6hr",
                            "unknown_6hr,GABA_6hr,NEFM_pos_glut_6hr",
                            "NEFM_pos_glut_6hr,GABA_6hr,unknown_6hr")] #5717
gr[gr$peak_called_in %in% c("unknown_6hr,NEFM_pos_glut_6hr,NEFM_neg_glut_6hr",
                            "NEFM_pos_glut_6hr,unknown_6hr,NEFM_neg_glut_6hr",
                            "NEFM_neg_glut_6hr,NEFM_pos_glut_6hr,unknown_6hr",
                            "NEFM_neg_glut_6hr,unknown_6hr,NEFM_pos_glut_6hr",
                            "unknown_6hr,NEFM_neg_glut_6hr,NEFM_pos_glut_6hr",
                            "NEFM_pos_glut_6hr,NEFM_neg_glut_6hr,unknown_6hr")] #25373


# Count number of RNAseq cells in each of the above categories to check ####
integrated_labeled$broad.type.time.ident <- integrated_labeled$timextype.ident
integrated_labeled$broad.type.time.ident[integrated_labeled$timextype.ident %in%
                                           c("SST_pos_GABA_0hr", "SEMA3E_pos_GABA_0hr",
                                             "GABA_0hr")] <- "GABA_0hr"

integrated_labeled$broad.type.time.ident[integrated_labeled$timextype.ident %in%
                                           c("SST_pos_GABA_1hr", "SEMA3E_pos_GABA_1hr",
                                             "GABA_1hr")] <- "GABA_1hr"
integrated_labeled$broad.type.time.ident[integrated_labeled$timextype.ident %in%
                                           c("SST_pos_GABA_6hr", "SEMA3E_pos_GABA_6hr",
                                             "GABA_6hr")] <- "GABA_6hr"
unique(integrated_labeled$broad.type.time.ident)
for(i in unique(integrated_labeled$broad.type.time.ident)){
  print(i)
  print(sum(integrated_labeled$broad.type.time.ident == i))
}
