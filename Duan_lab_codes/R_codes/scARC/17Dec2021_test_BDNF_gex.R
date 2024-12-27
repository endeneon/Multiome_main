SCT_0v1 <- subset(aggr_filtered, subset = time.ident %in% c("0hr", "1hr")) 
vec <- SCT_0v1@assays$SCT@data[SCT_0v1@assays$SCT@data@Dimnames[[1]] == "BDNF"]
df <- data.frame(val = vec,
                 time = SCT_0v1$time.ident,
                 stringsAsFactors = F)
t.test(x = df$val[df$time == "0hr"], 
       y = df$val[df$time == "1hr"], 
       alternative = "t", 
       paired = F, 
       var.equal = F)
ggplot(df,
       aes(x = time,
           y = val,
           colour = time)) +
  geom_boxplot(width = 0.25) +
  geom_violin() +

  theme_classic()

mean(df$val[df$time == "0hr"]) # 0.06600509
mean(df$val[df$time == "1hr"]) # 0.2418952


SCT_1v6 <- subset(aggr_filtered, subset = time.ident %in% c("6hr", "1hr")) 
vec <- SCT_1v6@assays$SCT@data[SCT_1v6@assays$SCT@data@Dimnames[[1]] == "BDNF"]
df <- data.frame(val = vec,
                 time = SCT_1v6$time.ident,
                 stringsAsFactors = F)
t.test(x = df$val[df$time == "1hr"], 
       y = df$val[df$time == "6hr"], 
       alternative = "t", 
       paired = F, 
       var.equal = F)
ggplot(df,
       aes(x = time,
           y = val,
           colour = time)) +
  geom_boxplot(width = 0.25) +
  geom_violin() +
  
  theme_classic()

mean(df$val[df$time == "1hr"]) # 0.2418952
mean(df$val[df$time == "6hr"]) # 0.2478398


SCT_0v6 <- subset(aggr_filtered, subset = time.ident %in% c("6hr", "0hr")) 
vec <- SCT_0v6@assays$SCT@data[SCT_0v6@assays$SCT@data@Dimnames[[1]] == "BDNF"]
df <- data.frame(val = vec,
                 time = SCT_0v6$time.ident,
                 stringsAsFactors = F)
t.test(x = df$val[df$time == "0hr"], 
       y = df$val[df$time == "6hr"], 
       alternative = "t", 
       paired = F, 
       var.equal = F)
ggplot(df,
       aes(x = time,
           y = val,
           colour = time)) +
  geom_boxplot(width = 0.25) +
  geom_violin() +
  
  theme_classic()

mean(df$val[df$time == "0hr"]) # 0.06600509
mean(df$val[df$time == "6hr"]) # 0.2478398
