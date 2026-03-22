#Copyright 2026 Doug Speed.

#    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

###############################################

#' @title GigaPRS Step 2
#
#' @description This function uses genotypes and SNP effect sizes to compute PRS for target individuals, possibly taking into account inferred ancestry
#
#' @details You need to provide genotypes (in PLINK 1 binary format) and either a single-ancestry or multi-ancestry scorefile
#
#' @param bedstem The stem of the genotype files (e.g., if the genotypes are stored in file.bed, file.bim and file.fam, use bedstem="file")
#' @param scores A file containing single-ancestry SNP effect sizes (e.g., the output from running GigaPRS Step 1 with one focal ancestry)
#' @param multiscores A file containing multi-ancestry SNP effect sizes  (e.g., the output from running GigaPRS Step 1 with multiple focal ancestries)
#' @param outstem The desired prefix for the output files
#' @param allowMissing (TRUE or FALSE; default FALSE) Whether to allow missing genotypes for some score predictors
#' @param extractfile (Optional) A file containing a subset of predictors
#' @param keepfile (Optional) A file containing a subset of samples
#' @param phenofile (Optional) A file containing phenotypes for one or more samples
#' @param saveCounts (TRUE or FALSE; default FALSE) Whether to save number of non-missing values for each sample

#' @export

#' @examples
#' #These examples use the gwas results that come with the GigaPRS package
#' ex_gwas_file=system.file("extdata", "ex_gwas.txt.gz", package="MegaPRS")
#' 
#' #Example 1 - Provide only the two required arguments
#' format_sumstats(gwasfile=ex_gwas_file, outstem="ex_out")
#' 
#' #Example 2 - Run the function based on the on-screen instructions from Example 1
#' format_sumstats(gwasfile=ex_gwas_file, outstem="ex_out", NameCol="Predictor", A1Col="A1", A2Col="A2", EffectCol="Beta", SECol="SE", nCol="n", FreqCol="A1Freq")
GigaPRS_Step2=function(bedstem=NULL, scores=NULL, multiscores=NULL,outstem=NULL,allowMissing=FALSE,extractfile=NULL,keepfile=NULL,phenofile=NULL,saveCounts=FALSE)
{
################
#get start time
start_time=Sys.time()
cat(paste0("Start at ",start_time,"\n\n"))


################
#check arguments

if(is.null(bedstem))
{return(paste0("Error, you must use the argument bedstem to provide the stem of the genotype data, stored in PLINK 1 binary format (e.g., if the genotypes are stored in file.bed, file.bim and file.fam, use bedstem=\"file\""))}

bedfile=paste0(bedstem,".bed")
bimfile=paste0(bedstem,".bim")
famfile=paste0(bedstem,".fam")

if(file.exists(bedfile)==FALSE)
{return(paste0("Error, unable to find the bedfile ",bedfile, ", please check the argument bedstem"))}
if(file.exists(bimfile)==FALSE)
{return(paste0("Error, unable to find the bimfile ",bimfile, ", please check the argument bedstem"))}
if(file.exists(bedfile)==FALSE)
{return(paste0("Error, unable to find the famfile ",famfile, ", please check the argument bedstem"))}

if(is.null(scores)&&is.null(multiscores))
{return(paste0("Error, you must use either the argument scores or the argument multiscores to provide SNP effect sizes from GigaPRS Step 1 (use scores if you ran Step 1 with one focal ancestry, and multiscores if you ran Step 1 with multiple focal ancestries)"))}

if(is.null(outstem))
{return(paste0("Error, you must use the argument outstem to specify the stem for the output files"))}

if(allowMissing!=TRUE&&allowMissing!=FALSE)
{return(paste0("Error, allowMissing must be True or False (not ",allowMissing,")"))}

if(!is.null(extractfile)&&file.exists(extractfile)==FALSE)
{return(paste0("Error, unable to find the extractfile ",extractfile))}

if(!is.null(keepfile)&&file.exists(keepfile)==FALSE)
{return(paste0("Error, unable to find the keepfile ",keepfile))}

if(!is.null(phenofile)&&file.exists(phenofile)==FALSE)
{return(paste0("Error, unable to find the phenofile ",phenofile))}

if(saveCounts!=TRUE&&saveCounts!=FALSE)
{return(paste0("Error, saveCounts must be True or False (not ",saveCounts,")"))}


################
#set multiflag and scorefile

multiflag=is.null(scores)
if(multiflag==FALSE){scorefile=scores}
else{scorefile=multiscores}

if(multiflag==FALSE)    #read scorefile and set num_scores
{
cat(paste0("Reading scores from ",scorefile,"\n"))

all=read.table(scorefile,head=TRUE)
if(nrow(all)==0)
{return(paste0("Error, ", scorefile," contains no predictors"))}

if(ncol(all)<5)
{return(paste0("Error, ", scorefile," should have at least five columns (not ",ncol(all),")"))}

num_scores=ncol(all)-4
if(num_scores==1){cat(paste0("There is one score and ", nrow(all)," predictors\n\n"))}
if(num_scores>1){cat(paste0("There are ", num_scores," scores and ", nrow(all)," predictors\n\n"))}

preds=all[,1]
al1=all[,2]
al2=all[,3]
centres=all[,4]
effects=as.matrix(all[,4+1:num_scores],ncol=num_scores)

if(sum(centres<0,na.rm=T)>0||sum(centres>2,na.rm=T)>0)
{return(paste0("Error reading ", scorefile, "; all centres should be between 0 and 2 (or NA)"))}

if(sum(is.na(effects))>0)
{return(paste0("Error reading ", scorefile, "; ", sum(is.na(effects)), " effect sizes are NA"))}
}

if(multiflag==TRUE)    #read multiscorefile and set num_scores
{
cat(paste0("Reading scores from ",scorefile,"\n"))

all=read.table(scorefile,head=TRUE)
if(nrow(all)==0)
{return(paste0("Error, ", scorefile," contains no predictors"))}

if(ncol(all)<5)
{return(paste0("Error, ", scorefile," should have at least five columns (not ",ncol(all),")"))}

if(ncol(all)%%2!=1)
{return(paste0("Error, ", scorefile," should have an odd number of columns (not ",ncol(all),")"))}

num_scores=(ncol(all)-3)/2

preds=all[,1]
al1=all[,2]
al2=all[,3]
centres=as.matrix(all[,2+2*1:num_scores],ncol=num_scores)
effects=as.matrix(all[,3+2*1:num_scores],ncol=num_scores)

if(sum(centres<0,na.rm=T)>0||sum(centres>2,na.rm=T)>0)
{return(paste0("Error reading ", scorefile, "; all centres should be between 0 and 2 (or NA)"))}

miss1=which(is.na(centres))
miss2=which(is.na(effects))

diff=setdiff(miss1,miss2)
if(length(diff)>0)
{return(paste0("Error reading ", scorefile,"; ", length(diff)," predictors have effects but not centres"))}

diff=setdiff(miss2,miss1)
if(length(diff)>0)
{return(paste0("Error reading ", scorefile,"; ", length(diff)," predictors have centres but not effects"))}

miss_count=apply(centres,1,function(x) sum(is.na(x)))
all_present=which(miss_count==0)

if(length(all_present)==0)
{return(paste0("Error reading ", scorefile, "; no predictors have centres and effects for all scores"))}

if(length(all_present)<nrow(all))
{cat(paste0("Warning, only ", length(all_present)," predictors have centres and effects for all scores\n\n"))}

if(length(all_present)<100)
{cat(paste0("Warning, it will be hard to reliably esitmate ancestries\n\n"))}

#convenient to set NA effects to zero
effects[miss2]=0
}

rm(all)


################
#see if reducing scores

if(!is.null(extractfile))
{
cat(paste0("Reading extract predictors from ",extractfile,"\n"))
extract=read.table(extractfile,head=FALSE)

red=which(preds %in% extract[,1])

if(length(red)==0)
{return(paste0("Error, none of the ", nrow(extract), " predictors in ", extractfile," are also in ", scorefile))}

if(length(red)==nrow(extract))
{cat(paste0("All ", nrow(extract)," predictors are in ", scorefile,"\n\n"))}

if(length(red)<nrow(extract))
{cat(paste0("Warning, only ", length(red)," of the ", nrow(extract)," predictors in ", extractfile," are also in ", scorefile,"\n\n"))}

preds=preds[red]
al1=al1[red]
al2=al2[red]
if(multiflag==FALSE){centres=centres[red]}
else{centres=matrix(centres[red,],ncol=num_scores)}
effects=matrix(effects[red,],ncol=num_scores)
}


################
#read famfile and set num_samples

cat(paste0("Reading sample details from ",famfile,"\n"))
fam=read.table(famfile,head=FALSE)
num_samples=nrow(fam)
cat(paste0("In total, there are ",num_samples," samples\n\n"))


################
#which samples are we using

if(is.null(keepfile))
{use_samples=1:num_samples}
else
{
cat(paste0("Reading keep samples from ",keepfile,"\n"))
keep=read.table(keepfile,head=FALSE)
if(ncol(keep)<2)
{return(paste0("Error, ", keepfile," should have at least two columns (not ", ncol(keep),")"))}

use_samples=which(paste0(fam[,1],"___",fam[,2]) %in% paste0(keep[,1],"___",keep[,2]))

if(length(use_samples)==0)
{return(paste0("Error, none of the ", nrow(keep), " samples in ", keepfile," are also in ",famfile))}

if(length(use_samples)==nrow(keep))
{cat(paste0("All ", nrow(keep)," samples are in ", famfile,"\n\n"))}

if(length(use_samples)<nrow(keep))
{cat(paste0("Warning, only ", length(use_samples)," of the ", nrow(keep)," samples in ", keepfile," are also in ",famfile,"\n\n"))}
}


################
#set phenotypes

resps=rep(NA,length(use_samples))
if(!is.null(phenofile))
{
cat(paste0("Reading phenotypes from ",phenofile,"\n"))
phens=read.table(phenofile,head=FALSE)
if(ncol(phens)!=3)
{return(paste0("Error, ", phenofile," should have exactly three columns (not ", ncol(phens),")"))}

overlap=intersect(paste0(fam[use_samples,1],"___",fam[use_samples,2]),paste0(phens[,1],"___",phens[,2]))
found1=match(overlap,paste0(fam[use_samples,1],"___",fam[use_samples,2]))
found2=match(overlap,paste0(phens[,1],"___",phens[,2]))
resps[found1]=as.numeric(phens[found2,3])

if(length(overlap)==0)
{
if(is.null(keepfile)){return(paste0("Error, none of the ", nrow(phens), " samples in ", phenofile," are also in ",famfile))}
else{return(paste0("Error, none of the ", nrow(phens), " samples in ", phenofile," are also in ",keepfile))}
}

if(length(overlap)==length(use_samples))
{cat(paste0("All ", length(overlap)," samples have phenotypes\n\n"))}

if(length(overlap)<length(use_samples))
{cat(paste0("Warning, only ", length(overlap)," of the ", length(use_samples)," samples have phenotypes\n\n"))}
}


################
#read bimfile and set num_snps

cat(paste0("Reading predictor details from ",bimfile,"\n"))
bim=read.table(bimfile,head=FALSE)
num_snps=nrow(bim)
cat(paste0("In total, there are ",num_snps," predictors\n\n"))


################
#which predictors are in the data

use=intersect(bim[,2],preds)
if(length(use)==0)
{return(paste0("Error, none of the ", length(preds)," score predictors are in the data"))}

use_snps=match(use,bim[,2])

if(length(use)<length(preds))
{
if(allowMissing==FALSE)
{return(paste0("Error, only ", length(use)," of the ", length(preds)," score predictors are in the data (to contrinue regardless, use allowMissing=TRUE)"))}

cat(paste0("Warning, only ", length(use)," of the ", length(preds)," score predictors are in the data\n\n"))

red=match(use,preds)
preds=preds[red]
al1=al1[red]
al2=al2[red]
if(multiflag==FALSE){centres=centres[red]}
else{centres=matrix(centres[red,],ncol=num_scores)}
effects=matrix(effects[red,],ncol=num_scores)
}


################
#which predictors have consistent alleles, and which need to be flipped

comb1=paste0(al1,al2)
comb2=paste0(al2,al1)
comb3=paste0(bim[use_snps,5],bim[use_snps,6])

cons=which(comb1==comb3|comb2==comb3)

if(length(cons)==0)
{return(paste0("Error, none of the ", length(use)," score predictors have consistent alleles"))}

if(length(cons)<length(use))
{
if(allowMissing==FALSE)
{return(paste0("Error, only ", length(cons)," of the ", length(use)," score predictors have consistent alleles (to contrinue regardless, use allowMissing=TRUE)"))}

cat(paste0("Warning, only ", length(cons)," of the ", length(use)," score predictors have consistent alleles\n\n"))

preds=preds[cons]
al1=al1[cons]
al2=al2[cons]
if(multiflag==FALSE){centres=centres[cons]}
else{centres=matrix(centres[cons,],ncol=num_scores)}
effects=matrix(effects[cons,],ncol=num_scores)

use=use[cons]
use_snps=use_snps[cons]
}

flip=which(comb1[cons]!=comb3[cons])
if(length(flip)>0)
{
if(multiflag==FALSE){centres[flip]=2-centres[flip]}
else{centres[flip,]=2-centres[flip,]}
effects[flip,]=-effects[flip,]
}

################
#open bedfile and check headers

filecon=file(bedfile, "rb")
header=readBin(filecon, what = "raw", n = 3)
if(header[1]!=as.raw(0x6C)||header[2]!=as.raw(0x1B))
{return(paste0("Error, reading ", bedfile, "; incorrect magic numbers"))}

if(header[3]!=as.raw(0x01))
{return(paste0("Error, reading ", bedfile, "; not SNP-major mode (unsupported)"))}


################
#get some indexes
byte_indexes=ceiling(use_samples/4)
byte_shifts=2*(use_samples-1)%%4


################
#guesses and tallies start at NA or zero, work out read size, and set genotype codes
guesses=matrix(0,length(use_samples),ncol=num_scores)
if(saveCounts==FALSE){tallies=matrix(NA,length(use_samples),ncol=num_scores)}
else{tallies=matrix(0,length(use_samples),ncol=num_scores)}
bytes_per_snp=ceiling(num_samples/4)
codes=c(2,NA,1,0)


################
#STS and STM are used with multiscores
STS=matrix(0,nrow=num_scores,ncol=num_scores)
STM=matrix(0,nrow=length(use_samples),ncol=num_scores)


################
#loop through snps

location=1
num_snps_use=length(use_snps)
wcount=0
for(j in 1:num_snps_use)
{
if(j%%50000==1)
{cat(paste0("Processing SNP ", j, " out of ", num_snps_use,"\n"))}

if(location!=use_snps[j]) #get to location of next snp
{
offset=3+(use_snps[j]-1)*bytes_per_snp
seek(filecon, where = offset, origin = "start")
}

#read raw snp data
raw_bytes=readBin(filecon, what = "raw", n = bytes_per_snp)
location=use_snps[j]+1

#convert to 0, 1, 2 and NA for extract individuals we are using
bits=as.integer(raw_bytes[byte_indexes])
gens=codes[1+bitwAnd(bitwShiftR(bits, byte_shifts), 3)]

#find non-missing
nonmiss=!is.na(gens)
sum_nonmiss=sum(nonmiss)

if(sum_nonmiss>0)   #will use SNP
{
aver=sum(gens[nonmiss])/sum_nonmiss
gens[!nonmiss]=aver

if(multiflag==FALSE)    #update centre (if missing)
{
if(is.na(centres[j])){centres[j]=aver}
}

#update guesses (this works for NA multiscores because they have effect zero)
guesses=guesses+outer(gens, effects[j,])

if(saveCounts==TRUE)    #update tallies
{tallies[nonmiss,]=tallies[nonmiss,]+matrix(as.numeric(effects[j,]!=0),nrow=sum(nonmiss),ncol=num_scores,byrow=TRUE)}

if(multiflag==TRUE) #add on contributions to STS and STM (if all centres non-missing)
{
if(sum(is.na(centres[j,]))==0)
{
STS=STS+outer(centres[j,],centres[j,])
STM=STM+outer(gens,centres[j,])
}
}
}
else    #ignore SNP
{
if(wcount<5){cat(paste0("Warning, SNP ", preds[j], " is trivial (all values the same or missing)\n"))}
wcount=wcount+1
}
}
if(wcount>5){cat(paste0("In total, there were ", wcount, "trivial,SNPs\n"))}
cat(paste("\n"))


################
#shut file and subtract means from guesses

close(filecon)

if(multiflag==FALSE)    #only one set of means and no missing values
{
for(k in 1:num_scores)
{
guesses[,k]=guesses[,k]-sum(centres*effects[,k])
}
}
else   #multiple sets of means, and can be missing values
{
for(k in 1:num_scores)
{
nonmiss=which(!is.na(centres[,k]))
guesses[,k]=guesses[,k]-sum(centres[nonmiss,k]*effects[nonmiss,k])
}
}


################
#save profiles
profile=cbind(fam[use_samples,1:2],resps,NA)
names=c("ID1","ID2","Phenotype","BLANK")
for(k in 1:num_scores)
{
profile=cbind(profile,round(guesses[,k],6),tallies[,k])
names=c(names,c(paste0("Profile",k),paste0("Count",k)))
}
colnames(profile)=names

if(multiflag==FALSE){write.table(profile,paste0(outstem,".profile"),row=F,col=T,quote=F,sep="\t")}
else{write.table(profile,paste0(outstem,".profile.single"),row=F,col=T,quote=F,sep="\t")}


################
#deal with ancestries (if using) then make combined profile
if(multiflag==TRUE)
{
#get raw estimates, then force to be within [0,1] and sum to one, then round when one value is above 0.9
ests_raw=STM%*%solve(STS)
ests_trun=ests_raw

ests_trun[ests_trun<0]=0
ests_trun[ests_trun>1]=1
sums=apply(ests_trun,1,sum)
ests_trun=ests_trun/outer(sums,rep(1,num_scores))

maxes=apply(ests_trun,1,max)
find=which(maxes>0.9)
ests_trun[find,]=round(ests_trun[find,])

#save
ancestry=cbind(fam[use_samples,1:2],round(ests_raw,4),round(ests_trun,4))
colnames(ancestry)=c("ID1","ID2",paste0("Raw",1:num_scores),paste0("Processed",1:num_scores))
write.table(ancestry,paste0(outstem,".ancestry"),row=F,col=T,quote=F,sep="\t")

#get combined profile and save
combined=apply(guesses*ests_trun,1,sum)
profile=cbind(fam[use_samples,1:2],resps,NA,round(combined,6),NA)
colnames(profile)=c("ID1","ID2","Phenotype","BLANK","Combined_Profile","BLANK")
write.table(profile,paste0(outstem,".profile"),row=F,col=T,quote=F,sep="\t")
}


print(head(profile))


################
#get total time
end_time=Sys.time()
cat(paste0("End at ",end_time,"\n"))
cat(paste0("Total run time was ", round(difftime(end_time,start_time,units="hours"),2), " hours\n\n"))
}

