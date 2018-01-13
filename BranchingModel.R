

set.seed(1)

LowerBound=0
UpperBound=0.5

HistMax=20
IntervalWidth=.01
ObservedData=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\Pi_Values.txt", sep=",", header=FALSE)
ObservedData=ObservedData[,1]

LHoodStor=numeric()


X_Lim=c(0, .4)
Y_Lim=c(0, 40)

LHoodData_All=hist(ObservedData, breaks=seq(LowerBound,UpperBound,IntervalWidth), plot=FALSE)	
LHoodData=LHoodData_All$counts

MutationRates=c(1e-1, 1e-2, 1e-3, 1e-7)

StorBins=array(, c(1,length(seq(LowerBound,UpperBound,IntervalWidth))-1,4))


PlotNames=c("", "(A)", "(B)", "(C)", "(D)")
Adjust=.9
TextSize=2
AxisSize=2


PercentileCutoffs=function(Vector)
{
	Temp=quantile(Vector, p=c(.025, .5, .975))
	return(Temp)
}


NucleotideDiversity=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\Pi_Values.txt")
NucleotideDiversity=NucleotideDiversity[,1]
DataStor=hist(NucleotideDiversity, breaks=seq(LowerBound,UpperBound,IntervalWidth), plot=FALSE)

X_axis=seq(LowerBound,UpperBound,IntervalWidth)+.5*IntervalWidth
X_axis=X_axis[-length(X_axis)]

StorRoundsOfReplication=c(5,10,15)

Labels=numeric(12)
Labels[1:4] = c(expression(paste(Mu, "=", 10^{-1}, ", B=5", sep="")), expression(paste(Mu, "=", 10^{-2}, ", B=5", sep="")), expression(paste(Mu, "=", 10^{-3}, ", B=5", sep="")), expression(paste(Mu, "=", 10^{-7}, ", B=5", sep="")))
Labels[5:8] = c(expression(paste(Mu, "=", 10^{-1}, ", B=10", sep="")), expression(paste(Mu, "=", 10^{-2}, ", B=10", sep="")), expression(paste(Mu, "=", 10^{-3}, ", B=10", sep="")), expression(paste(Mu, "=", 10^{-7}, ", B=10", sep="")))
Labels[9:12] = c(expression(paste(Mu, "=", 10^{-1}, ", B=15", sep="")), expression(paste(Mu, "=", 10^{-2}, ", B=15", sep="")), expression(paste(Mu, "=", 10^{-3}, ", B=15", sep="")), expression(paste(Mu, "=", 10^{-7}, ", B=15", sep="")))

par(mfcol=c(4,3), mar=c(1,1,1,1), oma=c(6,6, 2, 2))
LabelSize=2
LabelLocation_x=1
LabelLocation_y=1.35

Y_Lim=140
X_Lim=0.485
FillColor="grey"
BorderColor="black"
LTY=2
Axis_Size=2
TextSize=2
Line_y=4
Line_x=4

Y_AxisBreaks=c(0,30,60,90,120)


Simulations=143
NucleotideDiversityStor=numeric(Simulations)


for (ii in 1:3)
{
	for (kk in 0:3)
	{
		for (i in 1:Simulations)
		{
			Loci=712

			StartingSize=1
			RoundsOfReplication=StorRoundsOfReplication[ii]

			MutationRate=MutationRates[kk+1]	#MutationRate=1e-7

			PopSize=StartingSize
			Haplotype=matrix(,StartingSize*2^RoundsOfReplication, Loci)

			Haplotype[1,]=0

			for (j in 1:RoundsOfReplication)
			{	
				Haplotype[(PopSize+1):(PopSize*2),]=Haplotype[1:PopSize,]

				for (k in 1:(PopSize*2))
				{
					Haplotype[k,]=Haplotype[k,]+rbinom(Loci, 1, MutationRate)
				}
				PopSize=PopSize*2
			}

			RevertedHaplotype=Haplotype %% 2

			AlleleFreq=colMeans(RevertedHaplotype)

			NucleotideDiversity = 1 - sum(AlleleFreq^2 + (1-AlleleFreq)^2)/Loci

			print(NucleotideDiversity)
			NucleotideDiversityStor[i]=NucleotideDiversity

		}
	
		HistStor=hist(NucleotideDiversityStor, density=FALSE, breaks=seq(0,.5,.01), xlim=c(0,X_Lim), ylim=c(0,Y_Lim), cex.axis=Axis_Size, axes=FALSE, main="", xlab="", ylab="")
		StorBins[1,,kk+1]=HistStor$counts


	#	points(X_axis, StorBins[,,kk+1]  , xlim=c(0,X_Lim), ylim=c(0,Y_Lim), cex.axis=Axis_Size, axes=FALSE, type="p", col=2, cex=1.5)
		box()
		
		if (ii==1)
		{
			axis(2, Y_AxisBreaks, cex.axis=Axis_Size, las=1)
		}
		if (kk==3)
		{
			axis(1, c(0,.1,.2,.3,.4), cex.axis=Axis_Size)
		}
		
		points(X_axis, DataStor$counts, pch=20, cex=2)
		#mtext("Count", 2, cex=TextSize, line=Line_y)

		par(new=TRUE)
		plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
		text(LabelLocation_x,LabelLocation_y, Labels[(kk+1)+4*(ii-1)], cex=LabelSize)


		if (kk==3 && ii==1)
		{
			mtext("Count", 2, cex=TextSize, line=Line_y, at=c(2.5))
		}

		if (kk==3 && ii==2)
		{
			mtext("Nucleotide diversity", 1, cex=TextSize, line=Line_x)
		}
	}

}




