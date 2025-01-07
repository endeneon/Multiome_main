sink("eval_covar_inputs.txt")
for (time in c("0hr","1hr","6hr")) {
	  for (cell_type in c("GABA","nmglut","npglut")) {
		      for (num_covar in 0:3) {
			            cat(paste0(time," ",cell_type," ",num_covar,"\n"))
    }
  }
}
sink()
