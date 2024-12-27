# Chuxuan Li 02/07/2022
# plot bargraphs of HOMER annotation data for number of peaks in each of the 
#types of genetic regions

library(ggplot2)
library(readr)
library(stringr)

GABA1 <- 
  read.table("logs/GABA1.txt", 
             header=TRUE, 
             quote="", 
             comment.char="")
ggplot(data = GABA1, aes(x = Annotation,
                         y = Npeaks)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ggtitle("GABA 1v0hr")

GABA6 <- 
  read.table("logs/GABA6.txt", 
             header=TRUE, 
             quote="", 
             comment.char="")
ggplot(data = GABA6, aes(x = Annotation,
                         y = Npeaks)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ggtitle("GABA 6v0hr")

Np1 <- 
  read.table("logs/Np1.txt", 
             header=TRUE, 
             quote="", 
             comment.char="")
ggplot(data = Np1, aes(x = Annotation,
                         y = Npeaks)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ggtitle("NEFM+ glut 1v0hr")

Np6 <- 
  read.table("logs/Np6.txt", 
             header=TRUE, 
             quote="", 
             comment.char="")
ggplot(data = Np6, aes(x = Annotation,
                       y = Npeaks)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ggtitle("NEFM+ glut 6v0hr")


Nm1 <- 
  read.table("logs/Nm1.txt", 
             header=TRUE, 
             quote="", 
             comment.char="")
ggplot(data = Nm1, aes(x = Annotation,
                       y = Npeaks)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ggtitle("NEFM+ glut 1v0hr")


Nm6 <- 
  read.table("logs/Nm6.txt", 
             header=TRUE, 
             quote="", 
             comment.char="")
ggplot(data = Nm6, aes(x = Annotation,
                       y = Npeaks)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ggtitle("NEFM- glut 6v0hr")
