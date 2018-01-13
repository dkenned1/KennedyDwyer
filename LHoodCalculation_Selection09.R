



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


Model1_LHoodStor=numeric()
Model2_LHoodStor=numeric()
Model3_LHoodStor=numeric()
Model4_LHoodStor=numeric()

#for (zzz in 1)
for (zzz in 1:100)
{
print(zzz)

set.seed(zzz)

ReadsPerLocusStor=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CoverageMatrix.txt", sep=",")
#ReadsPerLocusStor[which(ReadsPerLocusStor==0)]=1
#ReadsPerLocusStor=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CoverageMatrix_old.txt", sep=",")


PopDiversity=numeric()


#LowerBound=-0.005
#UpperBound=1.005

LowerBound=0
UpperBound=0.5

HistMax=20
IntervalWidth=.01
ObservedData=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\Pi_Values.txt", sep=",", header=FALSE)
ObservedData=ObservedData[,1]
#LHoodData=hist(ObservedData, breaks=seq(LowerBound,UpperBound,IntervalWidth))$counts

LHoodStor=numeric()


par(mfrow=c(1,3), mar=c(2,2,2,2), oma=c(6,6, 2, 2))

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
Seeds=Seeds[-c(2,3,4,5,6,10,14,17,18,20,21,24,25,27,30,31,34,38,43,47,59,63,68,72,80,81,87,88,91,92,96,98)]		#If s0.9

StorBins=array(, c(length(Seeds),length(seq(LowerBound,UpperBound,IntervalWidth))-1,4))


#PlotNames=c("(A -- data)", "(B -- M1)", "(C -- M2)", "(D -- M3)")
PlotNames=c("", "(A)", "(B)", "(C)")
Adjust=.9
TextSize=2
AxisSize=2

#text (.2, max(LHoodData_All$density)*Adjust, PlotNames[1], cex=TextSize)

#mtext("Density", side=2, line=3, cex=AxisSize)

#NumStrains=143
NumStrains=1000

for (k in 1:length(Seeds))
{
	for (kk in 3)
	{	
		Model=kk

		PersistingStrains=numeric()
		jj=1

		SimulationStrainData=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CodesFromGrid\\OutputFiles\\HostEvolutionMigration_Selection_Model", kk, "_s0.9_Seed", Seeds[k],".txt", sep=""), skip=300, nrows=NumStrains, sep=" ", colClasses="numeric")
		SimulationStrainData=SimulationStrainData[,!is.na(SimulationStrainData[1,])]

		Range=apply(SimulationStrainData, 2, range)

		PersistingStrains[jj]=sum(Range[2,]>0, na.rm=TRUE)
		jj=jj+1

		MajorStrain=apply(SimulationStrainData, 1, max)		

		SegregatingSites=712

		AlleleFreq=runif(SegregatingSites, min=7/NumStrains, max=(NumStrains-7)/NumStrains)		##Old was 0.05 and 0.95

		Strains=ncol(SimulationStrainData)

		Genotype=matrix(,Strains, SegregatingSites)

		for (i in 1:SegregatingSites)
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


		for (m in 1:nrow(TempData))
		{
			ReadsPerLocus=ReadsPerLocusStor[,sample(ncol(ReadsPerLocusStor),1)]
			TempData[m,]=SequencingError(TempData[m,])
		}

		TempData[which(TempData<.5)]=abs(TempData[which(TempData<.5)]-1)

		PopDiversity=apply(TempData, 1, Pi_calculator)

#		HistStor=hist(PopDiversity, freq=TRUE, breaks=seq(LowerBound,UpperBound,IntervalWidth), xlim=X_Lim, ylim=Y_Lim, xlab="", main="")	
		HistStor=hist(PopDiversity, freq=TRUE, breaks=seq(LowerBound,UpperBound,IntervalWidth), xlim=X_Lim, xlab="", ylab="", main="")	

#		print(PopDiversity)

		Sims=sum(HistStor$counts)
		Bins=length(HistStor$counts)

		StorBins[k,,kk+1]=HistStor$counts

#		Probs=HistStor$density * IntervalWidth * ((Sims-1)/Sims) + 1/(Sims)*IntervalWidth
		Probs=HistStor$count/(Sims+length(seq(LowerBound,UpperBound,IntervalWidth))-1) + 1/(Sims+length(seq(LowerBound,UpperBound,IntervalWidth))-1)


		LHood=dmultinom(LHoodData, prob=Probs, log=TRUE)
		LHoodStor[length(LHoodStor)+1]=LHood

		text (.2, max(HistStor$density)*Adjust, PlotNames[kk+2], cex=TextSize)

		if (kk==0)
		{
			mtext("Density", side=2, line=3, cex=AxisSize)
		}
		if (kk==1)
		{
			mtext("Nucleotide diversity", side=1, line=3, cex=AxisSize)
		}
	}
	print(k)
}


TempMaxLHood=max(LHoodStor)
MeanLHood= log(mean(exp(LHoodStor-TempMaxLHood)))+TempMaxLHood


print(MeanLHood)

PercentileCutoffs=function(Vector)
{
	Temp=quantile(Vector, p=c(.025, .5, .975))
	return(Temp)
}

#Bounds1=apply(StorBins[,,1], 2, PercentileCutoffs)
#Bounds2=apply(StorBins[,,2], 2, PercentileCutoffs)
#Bounds3=apply(StorBins[,,3], 2, PercentileCutoffs)
Bounds4=apply(StorBins[,,4], 2, PercentileCutoffs)

NucleotideDiversity=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\Pi_Values.txt")
NucleotideDiversity=NucleotideDiversity[,1]
DataStor=hist(NucleotideDiversity, breaks=seq(LowerBound,UpperBound,IntervalWidth), plot=FALSE)

X_axis=seq(LowerBound,UpperBound,IntervalWidth)+.5*IntervalWidth
X_axis=X_axis[-length(X_axis)]



par(mfrow=c(3,1), mar=c(1,1,1,1), oma=c(6,6, 2, 2))
LabelSize=3
LabelLocation_x=1
LabelLocation_y=1.25

Y_Lim=95
X_Lim=0.4
FillColor="grey"
BorderColor="black"
LTY=2
Axis_Size=2
TextSize=2
Line_y=4
Line_x=4

Y_AxisBreaks=c(0,30,60,90)

if (0)
{
	plot(X_axis, Bounds1[1,]  , xlim=c(0,X_Lim), ylim=c(0,Y_Lim), type="l", col=0, cex.axis=Axis_Size, axes=FALSE)
	box()
	axis(2, Y_AxisBreaks, cex.axis=Axis_Size, las=1)
	polygon(c(X_axis, rev(X_axis)), c(Bounds1[3,], rev(Bounds1[1,])), col = FillColor, border = NA)
	lines(X_axis, Bounds1[1,], type="l", col=BorderColor, lty=LTY)
	lines(X_axis, Bounds1[3,], type="l", col=BorderColor, lty=LTY)
	points(X_axis, DataStor$counts, pch=20, cex=2)
	#mtext("Count", 2, cex=TextSize, line=Line_y)

	par(new=TRUE)
	plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
	text(LabelLocation_x,LabelLocation_y, "(A)", cex=LabelSize)
}

if (0)
{
	#plot(X_axis, Bounds2[1,]  , xlim=c(0,X_Lim), ylim=c(0,45), type="l", col=0, cex.axis=Axis_Size, axes=FALSE)
	plot(X_axis, Bounds2[1,]  , xlim=c(0,X_Lim), ylim=c(0,Y_Lim), type="l", col=0, cex.axis=Axis_Size, axes=FALSE)
	box()
	axis(2, Y_AxisBreaks, cex.axis=Axis_Size, las=1)
	polygon(c(X_axis, rev(X_axis)), c(Bounds2[3,], rev(Bounds2[1,])), col = FillColor, border = NA)
	lines(X_axis, Bounds2[1,], type="l", col=BorderColor, lty=LTY)
	lines(X_axis, Bounds2[3,], type="l", col=BorderColor, lty=LTY)
	points(X_axis, DataStor$counts, pch=20, cex=2)
	#axis(1, c(0,.1,.2,.3,.4), cex.axis=Axis_Size)
	#mtext("Nucleotide diversity", 1, cex=TextSize, line=Line_x)
	mtext("Count", 2, cex=TextSize, line=Line_y)

	par(new=TRUE)
	plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
	text(LabelLocation_x,LabelLocation_y, "(B)", cex=LabelSize)
}
if (0)
{
	#plot(X_axis, Bounds3[1,]  , xlim=c(0,X_Lim), ylim=c(0,90), type="l", col=0, cex.axis=Axis_Size, axes=FALSE)
	plot(X_axis, Bounds3[1,]  , xlim=c(0,X_Lim), ylim=c(0,Y_Lim), type="l", col=0, cex.axis=Axis_Size, axes=FALSE)
	box()
	axis(2, Y_AxisBreaks, cex.axis=Axis_Size, las=1)
	axis(1, c(0,.1,.2,.3,.4), cex.axis=Axis_Size)
	polygon(c(X_axis, rev(X_axis)), c(Bounds3[3,], rev(Bounds3[1,])), col = FillColor, border = NA)
	lines(X_axis, Bounds3[1,], type="l", col=BorderColor, lty=LTY)
	lines(X_axis, Bounds3[3,], type="l", col=BorderColor, lty=LTY)
	points(X_axis, DataStor$counts, pch=20, cex=2)
	mtext("Nucleotide diversity", 1, cex=TextSize, line=Line_x)

	par(new=TRUE)
	plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")
	text(LabelLocation_x,LabelLocation_y, "(C)", cex=LabelSize)
}

plot(X_axis, Bounds4[1,]  , xlim=c(0,X_Lim), ylim=c(0,Y_Lim), type="l", col=0, cex.axis=Axis_Size, axes=FALSE)
box()
axis(2, Y_AxisBreaks, cex.axis=Axis_Size, las=1)
polygon(c(X_axis, rev(X_axis)), c(Bounds4[3,], rev(Bounds4[1,])), col = FillColor, border = NA)
lines(X_axis, Bounds4[1,], type="l", col=BorderColor, lty=LTY)
lines(X_axis, Bounds4[3,], type="l", col=BorderColor, lty=LTY)
points(X_axis, DataStor$counts, pch=20, cex=2)
axis(1, c(0,.1,.2,.3,.4), cex.axis=Axis_Size)
mtext("Nucleotide diversity", 1, cex=TextSize, line=Line_x)
mtext("Count", 2, cex=TextSize, line=Line_y)

par(new=TRUE)
plot(1,1, col=0, axes=FALSE, main="", xlab="", ylab="")


#MeanLHoodModel1=log(mean(exp(LHoodStor[(0:68 * 3)+1])))
#MeanLHoodModel2=log(mean(exp(LHoodStor[1:66])))
#MeanLHoodModel3=log(mean(exp(LHoodStor[(0:68 * 3)+3])))
MeanLHoodModel4=log(mean(exp(LHoodStor[1:68])))


#print(MeanLHoodModel1)
#print(MeanLHoodModel2)
#print(MeanLHoodModel3)
print(MeanLHoodModel4)

#Model1_LHoodStor[zzz]=MeanLHoodModel1
#Model2_LHoodStor[zzz]=MeanLHoodModel2
#Model3_LHoodStor[zzz]=MeanLHoodModel3
Model4_LHoodStor[zzz]=MeanLHoodModel4

}
