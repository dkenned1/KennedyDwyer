

ParamSet=""
#ParamSet="_Short"

par(mfcol=c(2,4))

GenomeSize=161048
MinimumCoverage=100
CutoffForUnanimous=0.975
MinimumNumberOfCalls=7		# Using 10 gives a nicer picture.  As the number goes up the curve shifts to the right.

Strains=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\SequenceNames.txt")
Strains=Strains[,1]
SegregatingSites=numeric()
IntermediateFrequencySites=numeric()

PerSequenceAlternative=numeric()
PerSequenceIntermediate=numeric()

CoverageMatrix=matrix(0,length(Strains), GenomeSize)
RefMatrix=matrix(0,length(Strains), GenomeSize)
AltMatrix=matrix(0,length(Strains), GenomeSize)


for (i in 1:length(Strains))
{
	VCF=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\Bowtie2Output\\", Strains[i], ParamSet, ".vcf", sep=""), sep="\t", header=TRUE)

	for (j in 1:nrow(VCF))
	{
		Index=VCF[j,2]
		Output=as.numeric(unlist(strsplit(as.character(VCF[j,5]), ":")))
		CoverageMatrix[i,Index]=Output[2]
		RefMatrix[i,Index]=Output[3]
		AltMatrix[i,Index]=Output[4]
	}	
	print(i)
}



GoodCoverageGenomes=which(rowMeans(CoverageMatrix)>100)
GoodCoverageGenomes=GoodCoverageGenomes[-140]
GoodCoverageGenomes=GoodCoverageGenomes[-116]

SegregatingSites=numeric()
for (i in 1:length(GoodCoverageGenomes))
{
	SegregatingSites=c(SegregatingSites, which(AltMatrix[GoodCoverageGenomes[i],]/CoverageMatrix[GoodCoverageGenomes[i],]>CutoffForUnanimous & CoverageMatrix[GoodCoverageGenomes[i],]>MinimumCoverage))
}
ReferenceSites=numeric()
for (i in 1:length(Strains))
{
	ReferenceSites=c(ReferenceSites, which(RefMatrix[GoodCoverageGenomes[i],]/CoverageMatrix[GoodCoverageGenomes[i],]>CutoffForUnanimous & CoverageMatrix[GoodCoverageGenomes[i],]>MinimumCoverage))
}

AltCalls=names(which(table(SegregatingSites)>MinimumNumberOfCalls))
ReferenceCalls=names(which(table(ReferenceSites)>MinimumNumberOfCalls))

SegregatingSites=intersect(AltCalls, ReferenceCalls)

SegregatingSites=as.numeric(SegregatingSites)

### Added here so that it could be added to plots
NucleotideDiversity=numeric()
for (i in 1:length(GoodCoverageGenomes))
{
	AlleleFreq=RefMatrix[GoodCoverageGenomes[i],SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]
	NucleotideDiversity[i]=mean(2*AlleleFreq*(1-AlleleFreq), na.rm=TRUE)
}
###

par(mfrow=c(4,5))

for (i in 1:20)
{
#	TempSegregatingSites=which(CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]>MinimumCoverage)
#	hist(RefMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]], xlim=c(0,1), breaks=seq(-.005, 1.005, .01), main=Strains[GoodCoverageGenomes[i]])

	FocalGenome=GoodCoverageGenomes[i]
	TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

#	hist(RefAlleleFreq, xlim=c(0,1), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))

}



for (i in 21:40)
{
#	TempSegregatingSites=which(CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]>MinimumCoverage)
#	hist(RefMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]], xlim=c(0,1), breaks=seq(-.005, 1.005, .01), main=Strains[GoodCoverageGenomes[i]])

	FocalGenome=GoodCoverageGenomes[i]
	TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

#	hist(RefAlleleFreq, xlim=c(0,1), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
}


for (i in 41:60)
{
#	TempSegregatingSites=which(CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]>MinimumCoverage)
#	hist(RefMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]], xlim=c(0,1), breaks=seq(-.005, 1.005, .01), main=Strains[GoodCoverageGenomes[i]])

	FocalGenome=GoodCoverageGenomes[i]
	TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

#	hist(RefAlleleFreq, xlim=c(0,1), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
}


for (i in 61:80)
{
#	TempSegregatingSites=which(CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]>MinimumCoverage)
#	hist(RefMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]], xlim=c(0,1), breaks=seq(-.005, 1.005, .01), main=Strains[GoodCoverageGenomes[i]])

	FocalGenome=GoodCoverageGenomes[i]
	TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

#	hist(RefAlleleFreq, xlim=c(0,1), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
}


for (i in 81:100)
{
#	TempSegregatingSites=which(CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]>MinimumCoverage)
#	hist(RefMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]], xlim=c(0,1), breaks=seq(-.005, 1.005, .01), main=Strains[GoodCoverageGenomes[i]])

	FocalGenome=GoodCoverageGenomes[i]
	TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

#	hist(RefAlleleFreq, xlim=c(0,1), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
}

for (i in 101:120)
{
#	TempSegregatingSites=which(CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]>MinimumCoverage)
#	hist(RefMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]], xlim=c(0,1), breaks=seq(-.005, 1.005, .01), main=Strains[GoodCoverageGenomes[i]])

	FocalGenome=GoodCoverageGenomes[i]
	TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

#	hist(RefAlleleFreq, xlim=c(0,1), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
}

for (i in 121:140)
{
#	TempSegregatingSites=which(CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]>MinimumCoverage)
#	hist(RefMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]], xlim=c(0,1), breaks=seq(-.005, 1.005, .01), main=Strains[GoodCoverageGenomes[i]])

	FocalGenome=GoodCoverageGenomes[i]
	TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

#	hist(RefAlleleFreq, xlim=c(0,1), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
}


for (i in 141:143)
{
#	TempSegregatingSites=which(CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]>MinimumCoverage)
#	hist(RefMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites[TempSegregatingSites]], xlim=c(0,1), breaks=seq(-.005, 1.005, .01), main=Strains[GoodCoverageGenomes[i]])

	FocalGenome=GoodCoverageGenomes[i]
	TempCoverage=RefMatrix[FocalGenome,SegregatingSites]+AltMatrix[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrix[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrix[FocalGenome,SegregatingSites], AltMatrix[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

#	hist(RefAlleleFreq, xlim=c(0,1), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main=paste(Strains[FocalGenome], ", pi=", round(NucleotideDiversity[i], digits=2), sep=""))
}



###############Doesn't work because these genomes were removed just before plotting histograms
#Duplicates1=c(116, 75)
#Duplicates2=c(140, 77)

#plot(RefMatrix[GoodCoverageGenomes[Duplicates1[1]],SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[Duplicates1[1]],SegregatingSites], RefMatrix[GoodCoverageGenomes[Duplicates1[2]],SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[Duplicates1[2]],SegregatingSites])
#plot(RefMatrix[GoodCoverageGenomes[Duplicates2[1]],SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[Duplicates2[1]],SegregatingSites], RefMatrix[GoodCoverageGenomes[Duplicates2[2]],SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[Duplicates2[2]],SegregatingSites])
#plot(RefMatrix[GoodCoverageGenomes[Duplicates1[1]],SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[Duplicates1[1]],SegregatingSites], RefMatrix[GoodCoverageGenomes[Duplicates2[2]],SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[Duplicates2[2]],SegregatingSites])


NucleotideDiversity=numeric()
for (i in 1:length(GoodCoverageGenomes))
{
	AlleleFreq=RefMatrix[GoodCoverageGenomes[i],SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[i],SegregatingSites]
	NucleotideDiversity[i]=mean(2*AlleleFreq*(1-AlleleFreq), na.rm=TRUE)
}
hist(NucleotideDiversity, breaks=seq(0, .405, .01), xlim=c(0,.3))


if (0)
{
	NucleotideDiversity=numeric()
	for (i in 1:length(GoodCoverageGenomes))
	{
		AlleleFreq=RefMatrix[GoodCoverageGenomes[i],-SegregatingSites]/CoverageMatrix[GoodCoverageGenomes[i],-SegregatingSites]
		NucleotideDiversity[i]=mean(2*AlleleFreq*(1-AlleleFreq), na.rm=TRUE)
	}
	hist(NucleotideDiversity, breaks=seq(0, .4, .001), xlim=c(0,.3))
}

##To generate a coverage matrix
if (0)
{
	Temp=t(CoverageMatrix[GoodCoverageGenomes,SegregatingSites])
	
	write.table(file="C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\CoverageMatrix.txt", Temp, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=",")
}




if(0)
{
	TempStrains=c("GeneratedData1", "TempMultiStrain")
	CoverageMatrixTemp=matrix(0,length(TempStrains), GenomeSize)
	RefMatrixTemp=matrix(0,length(TempStrains), GenomeSize)
	AltMatrixTemp=matrix(0,length(TempStrains), GenomeSize)

	for (i in 1:length(TempStrains))
	{
		VCFTemp=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\Bowtie2Output\\", TempStrains[i], ParamSet, ".vcf", sep=""), sep="\t", header=TRUE)

		for (j in 1:nrow(VCFTemp))
		{
			IndexTemp=VCFTemp[j,2]
			OutputTemp=as.numeric(unlist(strsplit(as.character(VCFTemp[j,5]), ":")))
			CoverageMatrixTemp[i,IndexTemp]=OutputTemp[2]
			RefMatrixTemp[i,IndexTemp]=OutputTemp[3]
			AltMatrixTemp[i,IndexTemp]=OutputTemp[4]
		}	
		print(i)
	}
	FocalGenome=1
	TempCoverage=RefMatrixTemp[FocalGenome,SegregatingSites]+AltMatrixTemp[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrixTemp[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrixTemp[FocalGenome,SegregatingSites], AltMatrixTemp[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

	par(mar=c(6,6,2,2))
	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main="", xlab="Allele frequency", ylab="Count", cex.lab=2)

	FocalGenome=2
	TempCoverage=RefMatrixTemp[FocalGenome,SegregatingSites]+AltMatrixTemp[FocalGenome,SegregatingSites]
	RefAlleleFreq=RefMatrixTemp[FocalGenome,SegregatingSites]/TempCoverage

	MajAlleleFreq=apply(cbind(RefMatrixTemp[FocalGenome,SegregatingSites], AltMatrixTemp[FocalGenome,SegregatingSites]), 1, max)/TempCoverage

	hist(MajAlleleFreq, xlim=c(0.5,1), ylim=c(0,length(SegregatingSites)), breaks=seq(-0, 1, .01), main="", xlab="Allele frequency", ylab="Count", cex.lab=2)


}