files<-unlist(as.character(system("ls SOC*.dat",intern=TRUE)));

file.create("./alpha_nonneutral.dat");
file.create("./betha_nonneutral.dat");

for(f in 1:length(files)){
	filename<-files[f];
	dat<-read.csv(filename,sep=" ",header=FALSE)
	dat<-dat[,-101]
	nsites<-length(dat)
	
	for(i in 1:nsites){
		bp<-unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(dat[,i]), ",")),"[{]")),"[}]"))[1]
		dp<-unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(dat[,i]), ",")),"[{]")),"[}]"))[2]
		ndp<-unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(dat[,i]), ",")),"[{]")),"[}]"))[3]
		mp<-unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(dat[,i]), ",")),"[{]")),"[}]"))[4]
		n<-unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(dat[,i]), ",")),"[{]")),"[}]"))[5]
		alpha<-as.numeric(bp)/as.numeric(n)
		betha<-as.numeric(dp)/as.numeric(n)
		cat(paste(filename,alpha,sep=" "),sep="\n",file="./alpha_nonneutral.dat",append=TRUE);
		cat(paste(filename,betha,sep=" "),sep="\n",file="./betha_nonneutral.dat",append=TRUE);
	}
}

alpha<-read.csv("./alpha_nonneutral.dat",sep=" ",header=FALSE)
alpha_nums<-((as.numeric(alpha$V2))[!is.nan((as.numeric(alpha$V2)))])[!is.infinite(((as.numeric(alpha$V2))[!is.nan((as.numeric(alpha$V2)))]))]
betha<-read.csv("./betha_nonneutral.dat",sep=" ",header=FALSE)
betha_nums<-((as.numeric(betha$V2))[!is.nan((as.numeric(betha$V2)))])[!is.infinite(((as.numeric(betha$V2))[!is.nan((as.numeric(betha$V2)))]))]

png("./alpha_nonneutral.png");
hist(alpha_nums,nclass=length(alpha_nums),xlim=c(0,0.05),xlab="Alpha values",main="Non-Neutral")
dev.off();

png("./betha_nonneutral.png");
hist(betha_nums,nclass=length(betha_nums),xlim=c(0,0.05),xlab="Betha values",main="Non-Neutral")
dev.off();
