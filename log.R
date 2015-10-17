source("./NvsBp_simpler.R")
ls()
alpha
system("ls SOC*.dat")
as.character(system("ls SOC*.dat"))
files<-as.character(system("ls SOC*.dat"))
files
files<-unlist(as.character(system("ls SOC*.dat")))
files
files<-unlist(as.character(system("ls SOC*.dat")));
files
files<-unlist(as.character(system("ls SOC*.dat",intern=TRUE)));
files
savehistory("./log.R")
