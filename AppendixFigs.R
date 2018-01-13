

##First run VCFAnalyzer2.R

NewCoverageMatrix=CoverageMatrix[which(rowMeans(CoverageMatrix)>100),]
NewCoverageMatrix=NewCoverageMatrix[-140,]
NewCoverageMatrix=NewCoverageMatrix[-116,]

NewRefMatrix=RefMatrix[which(rowMeans(CoverageMatrix)>100),]
NewRefMatrix=NewRefMatrix[-140,]
NewRefMatrix=NewRefMatrix[-116,]

NewAltMatrix=AltMatrix[which(rowMeans(CoverageMatrix)>100),]
NewAltMatrix=NewAltMatrix[-140,]
NewAltMatrix=NewAltMatrix[-116,]

NewStrains=Strains[which(rowMeans(CoverageMatrix)>100)]
NewStrains=NewStrains[-140]
NewStrains=NewStrains[-116]

NewRefSeqMatrix=RefSeqMatrix[which(rowMeans(CoverageMatrix)>100),]
NewRefSeqMatrix=NewRefSeqMatrix[-140,]
NewRefSeqMatrix=NewRefSeqMatrix[-116,]
NewRefSeqMatrix[which(NewRefSeqMatrix=="")]="."

NewAltSeqMatrix=AltSeqMatrix[which(rowMeans(CoverageMatrix)>100),]
NewAltSeqMatrix=NewAltSeqMatrix[-140,]
NewAltSeqMatrix=NewAltSeqMatrix[-116,]
NewAltSeqMatrix[which(NewAltSeqMatrix=="")]="."

NewStrains=Strains[which(rowMeans(CoverageMatrix)>100)]
NewStrains=NewStrains[-140]
NewStrains=NewStrains[-116]

PercentLdMNPV=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\PercentLdMNPV.txt")
PercentLdMNPV=PercentLdMNPV[,1]
BLASTModel=lm(NucleotideDiversity~PercentLdMNPV)

####To count tri_or_more_allelic sites
Temp=(unlist(lapply(apply(NewAltSeqMatrix[, SegregatingSites], 2, unique), length)))
Temp=which(Temp>2)
apply(NewAltSeqMatrix[,SegregatingSites[Temp]], 2, table)


############################Histogram of mean coverage for each of our samples###############################################
hist(apply(NewCoverageMatrix, 1, mean), xlab="Mean genome coverage", ylab="Number of samples", main="Histogram of mean coverage", breaks=seq(200, 1500, 100))



############################Histogram of SNP locations, and mean coverage for each locus in the genome ###################### 
TextSize=2
AxisSize=1.5

par(mfcol=c(2,1), mar=c(0,0,0,0), oma=c(6,6,2,2))

hist(SegregatingSites, main="", ylab="", xlab="", breaks=seq(0,161048, 1000), axes=FALSE, ylim=c(0,20), xlim=c(0,161000))
axis(2, at=seq(0,20,5), las=1, cex.axis=AxisSize)

#mtext("Histogram of SNP locations", 3, cex=TextSize)
mtext("Segregating sites", 2, cex=TextSize, line=4)
 
 
plot(log10((apply(NewCoverageMatrix, 2, mean))+1), type="l", axes=FALSE, xlab="Genome location", ylab="Mean genome Coverage", ylim=c(0,3.5), xlim=c(0,161000))
axis(1, pretty(c(0, 161048)), cex.axis=AxisSize)

y <- c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5")
aty <- seq(0,3.5,.5)
labels <- sapply(aty,function(i)
				 as.expression(bquote(10^ .(i)))
				 )
axis(2,at=aty,labels=labels, las=1, cex.axis=AxisSize)

#text(SegregatingSites, 3.65, ".")
text(SegregatingSites, 3.65, "|")

mtext("Mean coverage", side=2, line=4, cex=TextSize)
mtext("Location in genome", side=1, line=3.5, cex=TextSize)




##############################Plot of linkage disequilibrium by distance######################################################

GenotypeMatrix= (NewRefMatrix/NewCoverageMatrix)>0.5
#GenotypeMatrix=NewAltSeqMatrix

LinkageMatrix=matrix(,choose(length(SegregatingSites),2), 6)

k=1
for (i in 1:(length(SegregatingSites)-1))
{
	for (j in (i+1):length(SegregatingSites))
	{
		LinkageMatrix[k,1] = (cor(GenotypeMatrix[,SegregatingSites[i]], GenotypeMatrix[,SegregatingSites[j]], use="pairwise.complete.obs"))^2
		LinkageMatrix[k,2] = min((SegregatingSites[j]-SegregatingSites[i]), abs(SegregatingSites[j]+SegregatingSites[i]-161046))
		
		LinkageMatrix[k,3] = SegregatingSites[i]
		LinkageMatrix[k,4] = SegregatingSites[j]

		p_i=mean(GenotypeMatrix[,SegregatingSites[i]]==0, na.rm=TRUE)
		p_j=mean(GenotypeMatrix[,SegregatingSites[j]]==0, na.rm=TRUE)
		
		Temp00 = mean(GenotypeMatrix[,SegregatingSites[i]]==0 & GenotypeMatrix[,SegregatingSites[j]]==0, na.rm=TRUE)
		Temp11 = mean(GenotypeMatrix[,SegregatingSites[i]]==1 & GenotypeMatrix[,SegregatingSites[j]]==1, na.rm=TRUE)
		Temp10 = mean(GenotypeMatrix[,SegregatingSites[i]]==1 & GenotypeMatrix[,SegregatingSites[j]]==0, na.rm=TRUE)
		Temp01 = mean(GenotypeMatrix[,SegregatingSites[i]]==0 & GenotypeMatrix[,SegregatingSites[j]]==1, na.rm=TRUE)
				
		LinkageMatrix[k,5]= Temp00*Temp11-Temp01*Temp10		 #D
		
		if (LinkageMatrix[k,5]>=0)
		{
			LinkageMatrix[k,6]= min(p_i*(1-p_j), (1-p_i)*p_j)
		}
		if (LinkageMatrix[k,5]<0)
		{
			LinkageMatrix[k,6]= max(-p_i*p_j, -(1-p_i)*(1-p_j))
		}
		k=k+1
	}
}

par(mfcol=c(2,1), mar=c(0.5,0.5,0.5,0.5), oma=c(6,6,2,2))

plot(LinkageMatrix[,2], LinkageMatrix[,5]/LinkageMatrix[,6], pch=".", main="", xlab="", ylab="", las=1, axes=FALSE)
#plot(LinkageMatrix[,2], LinkageMatrix[,5], pch=".", main="", xlab="", ylab="", las=1, axes=FALSE)
box()
axis(2,at=pretty(c(0,1)), las=1) 

mtext("D'", 2, line=3, cex=2, las=1)

Spline=smooth.spline(LinkageMatrix[,2], LinkageMatrix[,5]/LinkageMatrix[,6])
lines(Spline, col=2, lwd=2)



plot(LinkageMatrix[,2], LinkageMatrix[,1], pch=".", main="", xlab="", ylab="", las=1)
mtext("Pairwise minimum distance (base pairs)", 1, line=3, cex=2)
mtext(expression(R^2), 2, line=3, cex=2, las=1)

Spline=smooth.spline(LinkageMatrix[,2], LinkageMatrix[,1])
lines(Spline, col=2, lwd=2)

################################Counting InDels################################################################

IndelFinder=function(x)
{
	Temp=grep("[-+]", x)
	if (length(Temp)>6)
	{
		Indel=TRUE
	}
	else
	{
		Indel=FALSE
	}
	return(Indel)
}

IndelLocations=apply(NewAltSeqMatrix[,SegregatingSites], 2, IndelFinder)

sum(IndelLocations)

#NewAltSeqMatrix[,SegregatingSites[IndelLocations]]

##############################Calculating population-wide nucleotide diversity######################################################



x_squared=numeric(length(SegregatingSites))

for (j in 1:length(SegregatingSites))
{
	Alleles=unique(NewAltSeqMatrix[,SegregatingSites[j]])

	x_ij_temp=0

	for (k in 1:length(Alleles))
	{
		x_ij_temp= x_ij_temp + (mean(NewAltSeqMatrix[,SegregatingSites[j]]==Alleles[k], na.rm=TRUE))^2
	}
	
	x_squared[j]=x_ij_temp
}

PopWideHaplotypePi = 1-mean(x_squared)


##############################Pairwise haplotype similarity######################################################

PairwiseSimilarity=matrix(,nrow(NewAltSeqMatrix), nrow(NewAltSeqMatrix))

LocationInfo=read.csv("C:\\Users\\Dave\\Dropbox\\Chapter4\\Figures\\StrainID_CollectionSite_Date.csv")
LocationInfo=LocationInfo[order(LocationInfo$Collection.Site),]

Model=lm(LocationInfo$Nucleotide.diversity ~ LocationInfo$Collection.Site + LocationInfo$Collection.Year + LocationInfo$Collection.Week)

Temp1=sub("[0-9][0-9]_", "", NewStrains)
Temp2=sub("[0-9]_", "", Temp1)

OrderedNewAltSeqMatrix=matrix(,nrow(NewAltSeqMatrix), ncol(NewAltSeqMatrix))		#Will be rearranged

NewNewStrains=numeric(length(NewStrains))
for (i in 1:length(NewStrains))
{
	OrderedNewAltSeqMatrix[i,]=NewAltSeqMatrix[which(Temp2==LocationInfo[i,1]),]
}

for (i in 1:nrow(NewAltSeqMatrix))
{
	for (j in 1:nrow(NewAltSeqMatrix))
	{
		PairwiseSimilarity[i,j]=sum(OrderedNewAltSeqMatrix[i,SegregatingSites]==OrderedNewAltSeqMatrix[j,SegregatingSites])
	}
}

rownames(PairwiseSimilarity)=NewStrains
colnames(PairwiseSimilarity)=NewStrains

par(oma=c(2,8,2,4))
image(PairwiseSimilarity, col= grey(seq(1,0, length=256)), axes=FALSE, ylim=c(1+0.5/(length(NewStrains)-1),-0.5/(length(NewStrains)-1)))
box()

LineWidth=2

rect(-0.5/(length(NewStrains)-1), -0.5/(length(NewStrains)-1), 1.5/(length(NewStrains)-1), 1.5/(length(NewStrains)-1), border="red", lwd=LineWidth)	#2001 samples
rect(1.5/(length(NewStrains)-1), 1.5/(length(NewStrains)-1), 9.5/(length(NewStrains)-1), 9.5/(length(NewStrains)-1), border="red", lwd=LineWidth)	#2002 samples
rect(9.5/(length(NewStrains)-1), 9.5/(length(NewStrains)-1), 15.5/(length(NewStrains)-1), 15.5/(length(NewStrains)-1), border="red", lwd=LineWidth)	#2003 samples

rect(15.5/(length(NewStrains)-1), 15.5/(length(NewStrains)-1), 28.5/(length(NewStrains)-1), 28.5/(length(NewStrains)-1), border="red", lwd=LineWidth)		#2000 samples, AL
rect(28.5/(length(NewStrains)-1), 28.5/(length(NewStrains)-1), 45.5/(length(NewStrains)-1), 45.5/(length(NewStrains)-1), border="red", lwd=LineWidth)		#2000 samples, GL
rect(45.5/(length(NewStrains)-1), 45.5/(length(NewStrains)-1), 51.5/(length(NewStrains)-1), 51.5/(length(NewStrains)-1), border="red", lwd=LineWidth)		#2000 samples, JA
rect(51.5/(length(NewStrains)-1), 51.5/(length(NewStrains)-1), 74.5/(length(NewStrains)-1), 74.5/(length(NewStrains)-1), border="red", lwd=LineWidth)		#2000 samples, KZ
rect(74.5/(length(NewStrains)-1), 74.5/(length(NewStrains)-1), 97.5/(length(NewStrains)-1), 97.5/(length(NewStrains)-1), border="red", lwd=LineWidth)		#2000 samples, MA
rect(97.5/(length(NewStrains)-1), 97.5/(length(NewStrains)-1), 102.5/(length(NewStrains)-1), 102.5/(length(NewStrains)-1), border="red", lwd=LineWidth)		#2000 samples, MU
rect(102.5/(length(NewStrains)-1), 102.5/(length(NewStrains)-1), 142.5/(length(NewStrains)-1), 142.5/(length(NewStrains)-1), border="red", lwd=LineWidth)	#2000 samples, YS

Labels=c("Year 2001", "Year 2002", "Year 2003", "Year 2000, Site AL", "Year 2000, Site GL", "Year 2000, Site JA", "Year 2000, Site KZ", "Year 2000, Site MA", "Year 2000, Site MU", "Year 2000, Site YS")
LabelLocations = c(0.5, 5.5, 13, 22, 37, 48, 63, 86, 100, 122.5)/(length(NewStrains)-1) 
axis(3, at=LabelLocations, labels=Labels, las=2, cex.axis=.7)
axis(4, at=LabelLocations, labels=Labels, las=2, cex.axis=.7)

legend("left", fill=grey(seq(0,1, length=10)), legend=round(c(seq(max(PairwiseSimilarity),min(PairwiseSimilarity), length=10))), xpd=TRUE, inset=-.225, box.col=0, cex=.75)



##############################Calculate Fst between consensus sequences######################################################




p=matrix(0,length(SegregatingSites), 10)	##10 subpopulations
c=matrix(0,length(SegregatingSites), 10)
p_bar=numeric(length(SegregatingSites))		
F_ST=numeric(length(SegregatingSites))

SubpopulationBreaks=c(1,3,11,17,30,47,53,76,99,104,144)

for (i in 1:length(SegregatingSites))
{
	for (j in 1:(length(SubpopulationBreaks)-1))
	{
		p[i,j]= mean(OrderedNewAltSeqMatrix[SubpopulationBreaks[j]:(SubpopulationBreaks[j+1]-1),SegregatingSites[i]]==".", na.rm=TRUE)	
		c[i,j]= length(SubpopulationBreaks[j]:(SubpopulationBreaks[j+1]-1))/sum(!is.na(OrderedNewAltSeqMatrix[,SegregatingSites[i]])) 
	}
	p_bar[i]= mean(OrderedNewAltSeqMatrix[,SegregatingSites[i]]==".", na.rm=TRUE)	
	
	F_ST[i]=  p_bar[i]*(1-p_bar[i]) - sum(c[i,]* p[i,]*(1-p[i,]))
}

hist(F_ST, xlab=expression(F[ST]), ylab="Number of segregating sites")




	
	
	
