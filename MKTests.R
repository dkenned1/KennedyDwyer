

Stor_Output=matrix(, length(Samples), 6)

Samples=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\Samples143.txt")
Samples=Samples[,1]
Samples=as.character(Samples)

SynonymousSites=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\snpEff_latest_core\\snpEff\\SynonymousSites.txt")
SynonymousSites=SynonymousSites[,1]

NonSynonymousSites=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\snpEff_latest_core\\snpEff\\NonSynonymousSites.txt")
NonSynonymousSites=NonSynonymousSites[,1]


for (j in 1:length(Samples))
{
	TempData=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\Bowtie2Output\\", Samples[j], ".vcf", sep=""), header=TRUE)


	SynOutput=matrix(,length(SynonymousSites), 3)
	for (i in 1:length(SynonymousSites))
	{
		SynOutput[i,]=as.numeric(TempData[(TempData$Position==SynonymousSites[i]),7:9])
	}

	NonSynOutput=matrix(,length(NonSynonymousSites), 3)
	for (i in 1:length(NonSynonymousSites))
	{
		NonSynOutput[i,]=as.numeric(TempData[(TempData$Position==NonSynonymousSites[i]),7:9])
	}


	FixedSyn=sum(SynOutput[,c(1,3)], na.rm=TRUE)
	SegregatingSyn=sum(SynOutput[,2], na.rm=TRUE)

	FixedNonSyn=sum(NonSynOutput[,c(1,3)], na.rm=TRUE)
	SegregatingNonSyn=sum(NonSynOutput[,2], na.rm=TRUE)

	MKTable=matrix(,2,2)
	MKTable[1,1]=FixedSyn
	MKTable[1,2]=SegregatingSyn
	MKTable[2,1]=FixedNonSyn
	MKTable[2,2]=SegregatingNonSyn


	Results=fisher.test(MKTable)

	print(Samples[j])
	print(MKTable)
	print(Results$p.value)


	Stor_Output[j,] = c(Samples[j], Results$p.value, MKTable[1,1], MKTable[1,2], MKTable[2,1], MKTable[2,2])

}


# Define the function
ggd.qqplot = function(pvector, main=NULL, ...) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
		 xlab=expression(Expected~~-log[10](italic(p))),
		 ylab=expression(Observed~~-log[10](italic(p))),
		 xlim=c(0,max(e)), ylim=c(0,max(o)))
    abline(0,1, col="red")
}



ggd.qqplot(Stor_Output[,2], "QQ-plot, regression")


