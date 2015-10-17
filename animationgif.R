library(animation)
dat1<-read.csv("./output_species_001_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat2<-read.csv("./output_species_002_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat3<-read.csv("./output_species_003_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat4<-read.csv("./output_species_004_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat5<-read.csv("./output_species_005_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat6<-read.csv("./output_species_006_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat7<-read.csv("./output_species_007_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat8<-read.csv("./output_species_008_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat9<-read.csv("./output_species_009_seed_117_real_1_changes_0.dat",sep=" ",header=FALSE)
dat1<-dat1[,-101];dat2<-dat2[,-101];dat3<-dat3[,-101];dat4<-dat4[,-101];dat5<-dat5[,-101];
dat6<-dat6[,-101];dat7<-dat7[,-101];dat8<-dat8[,-101];dat9<-dat9[,-101];

time<-length(dat1[,1]);
cols<-1:10;rows<-1:10;
grid<-as.matrix(expand.grid(rows,cols));

basecol<-c("black","red","blue","green","brown","cyan","pink","orange","gray");
colors<-c(rep(basecol[1],100),rep(basecol[2],100),rep(basecol[3],100),rep(basecol[4],100),rep(basecol[5],100),rep(basecol[6],100),rep(basecol[7],100),rep(basecol[8],100),rep(basecol[9],100));
xypoints<-rbind(grid, cbind(grid[,1]+0.1,grid[,2]), cbind(grid[,1]+0.2,grid[,2]), cbind(grid[,1]+0.3,grid[,2]), cbind(grid[,1]+0.4,grid[,2]), cbind(grid[,1]+0.5,grid[,2]), cbind(grid[,1]+0.6,grid[,2]), cbind(grid[,1]+0.7,grid[,2]), cbind(grid[,1]+0.8,grid[,2]));

timeinterval<-0.001#time interval for each frame (in seconds)
nframes<-time;
myoptions<-ani.options(interval=timeinterval,outdir =   getwd());

saveGIF({for(t in 1:nframes){
	        sizes<-as.matrix(cbind(dat1[t,],dat2[t,],dat3[t,],dat4[t,],dat5[t,],dat6[t,],dat7[t,],dat8[t,],dat9[t,]))/1000;
	        #The "secret" is to plot the data for each timewindow (windowsize points each time)
	        plot(xypoints,cex=sizes,col=colors,pch=19,main=paste("Frame: ",t,sep=""));
	#       plot(xypoints,cex=sizes,col=colors,pch=19);
	        ani.pause();#a pause before draw the next frame. the time of the pause is defined by the parameter "interval" in function "ani.options"
	}
},movie.name = "animation.gif");
