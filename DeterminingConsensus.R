





par(mfrow=c(4,5), mar=c(2,2,2,2))


#####Find poorly sequenced genomes

PoorCoverage=100
GenomeLength=161048
CumulativeCoverage=numeric(GenomeLength)

FocalStrain=read.table("C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\SamplesUsed4.txt")
FocalStrain=as.character(FocalStrain[,1])

StorGenomes=matrix("", GenomeLength, length(FocalStrain))

for (i in 1:length(FocalStrain))
{
	TempA=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\AlleleCounts\\NumberOfA", FocalStrain[i], ".txt", sep=""))
	TempT=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\AlleleCounts\\NumberOfT", FocalStrain[i], ".txt", sep=""))
	TempC=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\AlleleCounts\\NumberOfC", FocalStrain[i], ".txt", sep=""))
	TempG=read.table(paste("C:\\Users\\Dave\\Dropbox\\Chapter4\\AlleleCounts\\NumberOfG", FocalStrain[i], ".txt", sep=""))

	Data=cbind(0, TempA, TempT, TempC, TempG)
	Coverage=rowSums(Data)

	MajorAlleleCount=apply(Data, 1, which.max)

	CumulativeCoverage=CumulativeCoverage+Coverage
	BaseCallCode=c("N", "A", "T", "C", "G")
	
	StorGenomes[,i]=BaseCallCode[MajorAlleleCount]
}

StorPoorCoverageSites=which((CumulativeCoverage/length(FocalStrain))<PoorCoverage)

write.table(StorGenomes, "C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\ConsensusWholeGenome.txt", sep=",", row.names=FALSE, col.names=FALSE)
write.table(StorGenomes[-StorPoorCoverageSites,], "C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\ConsensusStrongCalls.txt", sep=",", row.names=FALSE, col.names=FALSE)

write.table(StorPoorCoverageSites, "C:\\Users\\Dave\\Dropbox\\Chapter4\\R_Code\\PoorCoverageSites.txt", sep=",", row.names=FALSE, col.names=FALSE)

if (0)
{
	PairwiseConsensus=matrix(, length(FocalStrain), length(FocalStrain))
	for (i in 1:(length(FocalStrain)-1))
	{
		for (j in i:length(FocalStrain))
		{
			PairwiseConsensus[i,j]=sum(StorGenomes[-StorPoorCoverageSites,i]!=StorGenomes[-StorPoorCoverageSites,j])
		}

	}
}






