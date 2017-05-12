

##Run "VCFAnalyzer2.R" or open saved RDataFile "NucleotideDiversity.RData"##

Pi_calculator=function(Cadaver)
{	
	LociDiversity=2*Cadaver*(1-Cadaver)
	LociDiversity=LociDiversity[!is.na(LociDiversity)]

	NucleotideDiversity=sum(LociDiversity)/length(LociDiversity)
	
	return(NucleotideDiversity)
}

SequencingError=function(AlleleFreq)
{
	
	SampleSize=sample(ReadsPerLocus, ncol(TempData), replace=TRUE)
	
	ErrorRate=0.004			#rgamma(1, shape=5, scale=.001)			#0.005
	Temp=rbinom(length(SampleSize), SampleSize, p=AlleleFreq-(AlleleFreq-.5)*2*ErrorRate)
	return(Temp/SampleSize)
}


set.seed(2)

ReadsPerLocusStor=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CoverageMatrix.txt", sep=",")
#ReadsPerLocusStor[which(ReadsPerLocusStor==0)]=1
#ReadsPerLocusStor=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CoverageMatrix_old.txt", sep=",")


PopDiversity=numeric()


#LowerBound=-0.005
#UpperBound=1.005

LowerBound=0
UpperBound=1

HistMax=20
IntervalWidth=.01
ObservedData=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\Pi_Values.txt", sep=",", header=FALSE)
ObservedData=ObservedData[,1]
#LHoodData=hist(ObservedData, breaks=seq(LowerBound,UpperBound,IntervalWidth))$counts

LHoodStor=numeric()


par(mfrow=c(4,5), mar=c(2,2,2,2), oma=c(6,6, 2, 2))

X_Lim=c(0, .4)
Y_Lim=c(0, 40)

#LHoodData_All=hist(ObservedData, freq=FALSE, breaks=seq(LowerBound,UpperBound,IntervalWidth), xlim=X_Lim, ylim=Y_Lim, xlab="", ylab="Frequency", main="")	
LHoodData_All=hist(ObservedData, breaks=seq(LowerBound,UpperBound,IntervalWidth), plot=FALSE)	
LHoodData=LHoodData_All$counts

Model=0
Model=1
Model=2
r=c(.2)

Seeds=c(1:100)
Seeds=Seeds[-c(3,8,9,15,19,20,22,26,27,30,31,32,35,42,48,49,53,54,55,58,61,64,78,79,82,86,89,90,94,98,100)]

StorBins=array(, c(length(Seeds),length(seq(LowerBound,UpperBound,IntervalWidth))-1,3))


#PlotNames=c("(A -- data)", "(B -- M1)", "(C -- M2)", "(D -- M3)")
PlotNames=c("", "(A)", "(B)", "(C)")
Adjust=.9
TextSize=2
AxisSize=2

#text (.2, max(LHoodData_All$density)*Adjust, PlotNames[1], cex=TextSize)

#mtext("Density", side=2, line=3, cex=AxisSize)

NumStrains=143

for (k in 17)
{
	for (kk in 2)
	{	
		Model=kk

		PersistingStrains=numeric()
		jj=1

		SimulationStrainData=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CodesFromGrid\\OutputFiles\\HostEvolutionMigration_Model", kk, "_Seed", Seeds[k],".txt", sep=""), skip=300, nrows=NumStrains, sep=" ", colClasses="numeric")
		SimulationStrainData=SimulationStrainData[,!is.na(SimulationStrainData[1,])]

		Range=apply(SimulationStrainData, 2, range)

		PersistingStrains[jj]=sum(Range[2,]>0, na.rm=TRUE)
		jj=jj+1

		MajorStrain=apply(SimulationStrainData, 1, max)		

		NumSegregatingSites=712

		AlleleFreq=runif(NumSegregatingSites, min=7/NumStrains, max=(NumStrains-7)/NumStrains)		##Old was 0.05 and 0.95

		Strains=ncol(SimulationStrainData)

		Genotype=matrix(,Strains, NumSegregatingSites)

		for (i in 1:NumSegregatingSites)
		{
			TempOrder=sample(Strains)

			TempBreakPoint=max(which(cumsum(colSums(SimulationStrainData)[sample(Strains)]/nrow(SimulationStrainData))<AlleleFreq[i]))

			if (TempBreakPoint<0)
			{
				TempBreakPoint=1
			}
			if (TempBreakPoint>Strains-1)
			{
				TempBreakPoint=Strains-1
			}

			Genotype[TempOrder[1:TempBreakPoint], i]=1
			Genotype[TempOrder[(TempBreakPoint+1):Strains], i]=0

		}

		SimulatedData=matrix(0,nrow(SimulationStrainData),NumSegregatingSites)

		for (i in 1:nrow(SimulatedData))
		{
			for (j in 1:Strains)
			{
				SimulatedData[i,]=SimulatedData[i,]+SimulationStrainData[i,j]*Genotype[j,]
			}
		}

		SimulatedData[SimulatedData>1]=1
		SimulatedData[SimulatedData<0]=0


		TempData=SimulatedData

		TempData[which(TempData<.5)]=abs(TempData[which(TempData<.5)]-1)


		for (m in 1:nrow(TempData))
		{
			ReadsPerLocus=ReadsPerLocusStor[,sample(ncol(ReadsPerLocusStor),1)]
			TempData[m,]=SequencingError(TempData[m,])
		}

		TempData[which(TempData<.5)]=abs(TempData[which(TempData<.5)]-1)

		PopDiversity=apply(TempData, 1, Pi_calculator)

	}

	FocalHosts=sample(143, 20)
	for (i in 1:length(FocalHosts))
	{
		hist(TempData[FocalHosts[i],], main=PopDiversity[FocalHosts[i]], xlim=c(0.5,1), ylim=c(0,NumSegregatingSites), breaks=seq(0,1,.01))
	}

}






LabelSize=3
LabelLocation_x=1
LabelLocation_y=1.25


Cex_Pi_Size=2
TextSize=2.5
TextLine_y=3.5
TextLine_x=4.5

Temp1=cbind(PopDiversity, c(1:NumStrains))
Temp1=Temp1[order(Temp1[,1]),]

Temp2=cbind(NucleotideDiversity, c(1:NumStrains))
Temp2=Temp2[order(Temp2[,1]),]


par(mfcol=c(2,5))

i=1
hist(TempData[Temp1[i,2],], xlim=c(0.5,1), ylim=c(0,NumSegregatingSites), breaks=seq(-0, 1, .01), main="")
text(0.75, 400, bquote(pi== .(round(Temp1[i,1], digits=2))), cex=Cex_Pi_Size)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(A)", cex=LabelSize)



i=7
FocalGenome=GoodCoverageGenomes[i]
TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main="")
print(Pi_calculator(MajAlleleFreq))
text(0.75, 400, bquote(pi== .(round(NucleotideDiversity[i], digits=2))), cex=Cex_Pi_Size)


mtext("Number of Loci", side=2, line=TextLine_y, cex=TextSize, at=850)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(F)", cex=LabelSize)


##########################################################################


i=47	##412
hist(TempData[Temp1[i,2],], xlim=c(0.5,1), ylim=c(0,NumSegregatingSites), breaks=seq(-0, 1, .01), main="")
text(0.75, 400, bquote(pi== .(round(Temp1[i,1], digits=2))), cex=Cex_Pi_Size)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(B)", cex=LabelSize)

i=10
FocalGenome=GoodCoverageGenomes[i]
TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main="")
print(Pi_calculator(MajAlleleFreq))
text(0.75, 400, bquote(pi== .(round(NucleotideDiversity[i], digits=2))), cex=Cex_Pi_Size)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(G)", cex=LabelSize)




##########################################################################


i=65	##412
hist(TempData[Temp1[i,2],], xlim=c(0.5,1), ylim=c(0,NumSegregatingSites), breaks=seq(-0, 1, .01), main="")
text(0.75, 400, bquote(pi== .(round(Temp1[i,1], digits=2))), cex=Cex_Pi_Size)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(C)", cex=LabelSize)

i=64
FocalGenome=GoodCoverageGenomes[i]
TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main="")
print(Pi_calculator(MajAlleleFreq))
text(0.75, 400, bquote(pi== .(round(NucleotideDiversity[i], digits=2))), cex=Cex_Pi_Size)


mtext("Major allele frequency ", side=1, line=TextLine_x, cex=TextSize)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(H)", cex=LabelSize)

##########################################################################


i=87	##172
hist(TempData[Temp1[i,2],], xlim=c(0.5,1), ylim=c(0,NumSegregatingSites), breaks=seq(-0, 1, .01), main="")
text(0.75, 400, bquote(pi== .(round(Temp1[i,1], digits=2))), cex=Cex_Pi_Size)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(D)", cex=LabelSize)


i=43
FocalGenome=GoodCoverageGenomes[i]
TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main="")
print(Pi_calculator(MajAlleleFreq))
text(0.75, 400, bquote(pi== .(round(NucleotideDiversity[i], digits=2))), cex=Cex_Pi_Size)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(I)", cex=LabelSize)




##########################################################################





i=132
hist(TempData[Temp1[i,2],], xlim=c(0.5,1), ylim=c(0,NumSegregatingSites), breaks=seq(-0, 1, .01), main="")
text(0.75, 400, bquote(pi== .(round(Temp1[i,1], digits=2))), cex=Cex_Pi_Size)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(E)", cex=LabelSize)


i=129
FocalGenome=GoodCoverageGenomes[i]
TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main="")
print(Pi_calculator(MajAlleleFreq))
text(0.75, 400, bquote(pi== .(round(NucleotideDiversity[i], digits=2))), cex=Cex_Pi_Size)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
text(LabelLocation_x,LabelLocation_y, "(J)", cex=LabelSize)


