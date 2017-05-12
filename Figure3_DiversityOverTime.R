

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

PopDiversity=numeric()
NucleotideDiversityWithinHosts=numeric()
NucleotideDiversityTotal=numeric()
S_Stor=numeric()
I_Stor=numeric()
nu_Stor=numeric()
Final_S_Stor=numeric()

Seeds=c(7)

NumStrains=143
SegregatingSites=712
Strains=50

Years=120

NumCadavers=1000

LowerBound=0
UpperBound=1
IntervalWidth=0.01
X_Lim=c(0,0.35)
Y_Lim=c(0, 600)

par(mfrow=c(3,2), mar=c(4,4,2,2))

for (k in 1:length(Seeds))
{

	Genotype=matrix(,Strains, SegregatingSites)
	
	AlleleFreq=runif(SegregatingSites, min=0.05, max=0.95)		##Old was 0.05 and 0.95

	for (i in 1:SegregatingSites)
	{
		TempOrder=sample(Strains)

		TempBreakPoint=max(which(cumsum(rep(1, Strains)/Strains)<AlleleFreq[i]))

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

	for (kk in 1:Years)
	{
	
		PopulationDynamics=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CodesFromGrid\\OutputFiles\\Figure3Code_Seed", Seeds[k], ".txt",sep=""), skip=12 + (kk-1)*(NumCadavers+3), nrows=1, sep=",", colClasses="numeric")
		PopulationDynamics=unlist(PopulationDynamics)

		Final_S=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CodesFromGrid\\OutputFiles\\Figure3Code_Seed", Seeds[k], ".txt",sep=""), skip=14 + (kk-1)*(NumCadavers+3), nrows=1, sep=",", colClasses="character")
		Final_S=as.character(Final_S)
		Final_S=strsplit(Final_S, "=")
		Final_S=unlist(Final_S)
		Final_S=as.numeric(Final_S[2])
		
		S_Stor[kk]=PopulationDynamics[1]
		I_Stor[kk]=PopulationDynamics[2]
		nu_Stor[kk]=PopulationDynamics[3]
		Final_S_Stor[kk]=Final_S

		SimulationStrainData=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CodesFromGrid\\OutputFiles\\Figure3Code_Seed", Seeds[k], ".txt", sep=""), skip=15 + (kk-1)*(NumCadavers+3), nrows=NumCadavers, sep=" ", colClasses="numeric")
		SimulationStrainData=SimulationStrainData[,1:50]
		
		SimulatedData=matrix(0,nrow(SimulationStrainData),SegregatingSites)

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

		PopDiversity=apply(TempData, 1, Pi_calculator)
		
		NucleotideDiversityWithinHosts[kk]=mean(PopDiversity)
		NucleotideDiversityTotal[kk]=Pi_calculator(colMeans(SimulatedData))

		HistStor=hist(PopDiversity, freq=TRUE, breaks=seq(LowerBound,UpperBound,IntervalWidth), xlim=X_Lim, xlab="", ylab="", main=kk, ylim=Y_Lim)	

		text(0.15, 400, round(NucleotideDiversityWithinHosts[kk], digits=2))
		text(0.15, 300, round(NucleotideDiversityTotal[kk], digits=2))


	}
}	
	
par(mfrow=c(3,2), mar=c(4,4,2,2))

plot(log10(S_Stor), type="l")
plot(log10(I_Stor), type="l")
plot(NucleotideDiversityWithinHosts, type="l")
plot(NucleotideDiversityTotal, type="l")

plot(log10(nu_Stor), type="l")




par(mfcol=c(2,1), oma=c(3,3,3,3), mar=c(1,4,1,2), las=1)	

AxisSize=2
LineWidth=2
LabelSize=3

FocalYears=c(85:115)
plot(FocalYears, log10(S_Stor[FocalYears]), type="l", main="", ylab="", ylim=c(3,6), xlab="", axes=FALSE, lwd=LineWidth)
#lines(FocalYears, log10(I_Stor[FocalYears]), col=2)
box()

mtext("Population size", side=2, cex=AxisSize, las=0, line=3)

#aty <- axTicks(2)
aty=c(3:6)
labels <- sapply(aty,function(i)
				 as.expression(bquote(10^ .(i)))
				 )
axis(2,at=aty,labels=labels)
#mtext("Fraction infected", side=4, cex=AxisSize, las=0, line=3)

text(86,6, "(A)", cex=LabelSize, pos=1)

#par(new=TRUE)
#plot(FocalYears, (S_Stor[FocalYears]-Final_S_Stor[FocalYears])/S_Stor[FocalYears], type="l", ylim=c(0,.6), col=3, axes=FALSE, main="", xlab="", ylab="", lwd=LineWidth)

#aty=seq(0,.6, .2)
#axis(4, at=aty)
#text(119,.3, "Fraction infected", cex=AxisSize, srt=270, xpd=NA)

plot(FocalYears, NucleotideDiversityWithinHosts[FocalYears], type="l", main="", ylab="", xlab="", ylim=c(0,0.12), lwd=LineWidth, col="red")
mtext("Mean nucleotide diversity", side=2, cex=AxisSize, las=0, line=3)


par(new=TRUE)

plot(FocalYears, log10(I_Stor[FocalYears]), type="l", main="", ylab="", xlab="", axes=FALSE, ylim=c(3,7), lwd=LineWidth, col="blue")

text(119,5, "Infectious cadavers", cex=AxisSize, srt=270, xpd=NA)
aty=c(3:7)
labels <- sapply(aty,function(i)
				 as.expression(bquote(10^ .(i)))
				 )
axis(4,at=aty,labels=labels)

mtext("Year", side=1, cex=AxisSize, las=0, line=3)

text(86,7, "(B)", cex=LabelSize, pos=1)


