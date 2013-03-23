
# ------------------------------------------------------ #
# ---------- CLEAN OBJECTS AND LOAD PACKAGES ----------- #
# ------------------------------------------------------ #

ls()
rm(list=ls())
library(foreign)


setwd("/Users/Administrateur/Desktop/MultilocusTDT_Parallelization")

set.seed(pi)


# ------------------------------------------------------ #
# ---------- LOAD DATA FILES --------------------------- #
# ------------------------------------------------------ #

for (phenotype in 0:0){
for (village1 in 1:2){
for (village2 in village1:2){

# to set the path and the number of simulation for empirical pvalues

path = "/Users/Administrateur/Desktop/MultilocusTDT_Parallelization"
# nbsimul = 1e4
nbsimul = 1e1

ped=read.table(paste(path,"/MultilocusTDT_Output/Data/CheikhSNPs_2013_02_28.txt",sep=""), sep="\t", h=T)
ped=ped[ped$village==village1 | ped$village==village2,]

# ------------------------------------------------------ #
# ------- CHOICE OF MARKERS TO ANALYZE AND PHENOTYPE --- #
# ------------------------------------------------------ #

ped$phen=ped$pfaidbin

gene="HBB_HBG2_HBA"
L1="M_HbS"
L2="xmn1"
L3="alpha_37del"

ped$idlocus1=ped[,L1]
ped$idlocus2=ped[,L2]
ped$idlocus3=ped[,L3]

l=3
locus_on_X=c(0)

# ------------------------------------------------------ #
# ---------- DATA FRAME OF GENOTYPES ------------------- #
# ------------------------------------------------------ #

father=data.frame(unique(ped$fatherfbat))
names(father)=c("idDN")
father=unique(merge(father,ped[,c("fm","idDN",paste("idlocus",1:l,sep=""))], by="idDN"))
names(father)=c("fatherfbat","fm",paste("falocus",1:l,sep=""))

mother=data.frame(unique(ped$motherfbat))
names(mother)=c("idDN")
mother=unique(merge(mother,ped[,c("fm","idDN",paste("idlocus",1:l,sep=""))], by="idDN"))
names(mother)=c("motherfbat","fm",paste("molocus",1:l,sep=""))

gendata=unique(ped[,c("fm","idDN","fatherfbat","motherfbat","sex","phen",paste("idlocus",1:l,sep=""))])

gendata=merge(gendata, father, by=c("fm","fatherfbat"), all.x=T)
gendata=merge(gendata, mother, by=c("fm","motherfbat"), all.x=T)
rm(father,mother)

gendata[,(dim(gendata)[2]-3*l+1):dim(gendata)[2]][is.na(gendata[,(dim(gendata)[2]-3*l+1):dim(gendata)[2]])==TRUE]="0 0"
gendata=gendata[,c("fm","idDN","fatherfbat","motherfbat","sex","phen",paste("idlocus",1:l,sep=""),paste("falocus",1:l,sep=""),paste("molocus",1:l,sep=""))]
gendata=gendata[order(gendata$fm,gendata$idDN),]

gendata2=unique(gendata[is.na(gendata$phen)==FALSE & gendata$phen==phenotype & 
			gendata$falocus1!="0 0" & gendata$molocus1!="0 0" & gendata$idlocus1!="0 0" & 
		 	gendata$falocus2!="0 0" & gendata$molocus2!="0 0" & gendata$idlocus2!="0 0" & 
		 	gendata$falocus3!="0 0" & gendata$molocus3!="0 0" & gendata$idlocus3!="0 0" , 
			c("falocus1","molocus1","idlocus1","falocus2","molocus2","idlocus2","falocus3","molocus3","idlocus3","fatherfbat","motherfbat","idDN","sex")])

# ------------------------------------------------------ #

write.table(gendata2, file=paste(path,"/MultilocusTDT_Output/Data/TempData/gendata2.txt",sep=""), sep=" ", quote=F, row.names=F, col.names=F)
gendata2=read.table(paste(path,"/MultilocusTDT_Output/Data/TempData/gendata2.txt",sep=""), sep=" ")
names(gendata2)=c("fal1a1","fal1a2","mol1a1","mol1a2","chl1a1","chl1a2",
			"fal2a1","fal2a2","mol2a1","mol2a2","chl2a1","chl2a2",
			"fal3a1","fal3a2","mol3a1","mol3a2","chl3a1","chl3a2",
			"father","mother","child","sex")

######################################################
# --------- CANCELLING OFFSPRING GENOTYPE ---------- #
# --------    IN THE SIMULATED DATASET ------------- #
######################################################

gendata2saved = gendata2
# nbloci=3

# for (l in 1:nbloci){
# 	gendata2saved[,6*(l-1)+5] = NA
# 	gendata2saved[,6*(l-1)+6] = NA
# }

######################################################

# ------------------------------------------------------ #
# ---------- LISTE OF POSSIBLE K-UPLET ----------------- #
# ------------------------------------------------------ #

nbloci=3
nballeles=2

taballeles=matrix(NA,nbloci,nballeles)
rownames(taballeles)=c(paste("locus",1:nbloci,sep=""))
colnames(taballeles)=c(paste("allele",1:nballeles,sep=""))
for (l in 1:nbloci){
for (a in 1:nballeles){
taballeles[l,a]=unique(sort(c(as.matrix(gendata2[,((l-1)*6+1):(6*l)]))))[a]
}}
rm(a,l)

kuplet=NULL
l=0
for (l1 in taballeles[1,]){
for (l2 in taballeles[2,]){
for (l3 in taballeles[3,]){
	if (is.na(l1)==FALSE & is.na(l2)==FALSE & is.na(l3)==FALSE){
		l=l+1
		kuplet[l]=paste(l1,l2,l3, sep="")
	}
}}}
nbkuplet=length(kuplet)
rm(l,l1,l2,l3)

# ----------------------------------------------------------- #
# SIMULATION OF POSSIBLE CHILDREN FOR AMBIGUOUS TRANSMISSIONS #
# ----------------------------------------------------------- #

gendata2$countw=1
gendata2$realchild=1

for (n in 1:nrow(gendata2)){

nbdoubt=0
locusdoubt=0
for (l in 1:nbloci){
	if (gendata2[n,6*(l-1)+1]==gendata2[n,6*(l-1)+3] & gendata2[n,6*(l-1)+1]==gendata2[n,6*(l-1)+5] &
	    gendata2[n,6*(l-1)+2]==gendata2[n,6*(l-1)+4] & gendata2[n,6*(l-1)+2]==gendata2[n,6*(l-1)+6] &
	    gendata2[n,6*(l-1)+1]!=gendata2[n,6*(l-1)+2] &
	    gendata2$realchild[n]==1){
	nbdoubt=nbdoubt+1
	locusdoubt[nbdoubt]=l
	}
}

if (nbdoubt>0) {
gensimchild=NULL
for (p in (nbdoubt-1):0) {
gensimchild = c(gensimchild,rep(c(rep(1,2^p),rep(2,2^p)),2^(nbdoubt-p-1)))
}
gensimchild=matrix(gensimchild,2^nbdoubt,nbdoubt)

for (i in 1:2^nbdoubt){
gendata2=rbind(gendata2,gendata2[n,])
gendata2[nrow(gendata2),6*(locusdoubt-1)+5]=gensimchild[i,]
gendata2[nrow(gendata2),6*(locusdoubt-1)+6]=gensimchild[i,]
}

gendata2$countw[n]=0
gendata2$countw[(nrow(gendata2)-2^nbdoubt+1):nrow(gendata2)]=1/2^nbdoubt
gendata2$realchild[(nrow(gendata2)-2^nbdoubt+1):nrow(gendata2)]=0

}
}
rm(n,l,p,i)

dim(gendata2)
sum(gendata2$countw)

# ------------------------------------------------------ #
# --- TO COMPUTE MATRIX OF SUMULTANEOUS TRANSMISSION --- #
# ------------------------------------------------------ #

transmat=matrix(0,length(kuplet),length(kuplet))
rownames(transmat)=kuplet
colnames(transmat)=kuplet

for (n in 1:nrow(gendata2)){

gendata2$idDN=gendata2$child
gendata2=merge(gendata2, ped[,c("idDN","pfaid")], by=c("idDN"), all.x=T)
gendata2=gendata2[,-1]

kuplet_fa_T=NULL; kuplet_mo_T=NULL; kuplet_fa_NT=NULL; kuplet_mo_NT=NULL

for (l in 1:nbloci){

if ((length(setdiff(locus_on_X,l))==length(locus_on_X)) | (length(setdiff(locus_on_X,l))!=length(locus_on_X) & gendata2$sex[n]==2)){

for (i in taballeles[l,][is.na(taballeles[l,])==FALSE]){
for (j in i:max(taballeles[l,][is.na(taballeles[l,])==FALSE])){
for (u in taballeles[l,][is.na(taballeles[l,])==FALSE]){
for (v in u:max(taballeles[l,][is.na(taballeles[l,])==FALSE])){
	if (gendata2[n,6*(l-1)+1]==i & gendata2[n,6*(l-1)+2]==j & gendata2[n,6*(l-1)+3]==u & gendata2[n,6*(l-1)+4]==v & ((gendata2[n,6*(l-1)+5]==i & gendata2[n,6*(l-1)+6]==u) | (gendata2[n,6*(l-1)+5]==u & gendata2[n,6*(l-1)+6]==i))){
	kuplet_fa_T=paste(kuplet_fa_T,i, sep=""); kuplet_mo_T=paste(kuplet_mo_T,u, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,j, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,v, sep="")}
	else {
		if (gendata2[n,6*(l-1)+1]==i & gendata2[n,6*(l-1)+2]==j & gendata2[n,6*(l-1)+3]==u & gendata2[n,6*(l-1)+4]==v & ((gendata2[n,6*(l-1)+5]==i & gendata2[n,6*(l-1)+6]==v) | (gendata2[n,6*(l-1)+5]==v & gendata2[n,6*(l-1)+6]==i))){
		kuplet_fa_T=paste(kuplet_fa_T,i, sep=""); kuplet_mo_T=paste(kuplet_mo_T,v, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,j, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,u, sep="")}
		else {
			if (gendata2[n,6*(l-1)+1]==i & gendata2[n,6*(l-1)+2]==j & gendata2[n,6*(l-1)+3]==u & gendata2[n,6*(l-1)+4]==v & ((gendata2[n,6*(l-1)+5]==j & gendata2[n,6*(l-1)+6]==u) | (gendata2[n,6*(l-1)+5]==u & gendata2[n,6*(l-1)+6]==j))){
			kuplet_fa_T=paste(kuplet_fa_T,j, sep=""); kuplet_mo_T=paste(kuplet_mo_T,u, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,i, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,v, sep="")}
			else {
				if (gendata2[n,6*(l-1)+1]==i & gendata2[n,6*(l-1)+2]==j & gendata2[n,6*(l-1)+3]==u & gendata2[n,6*(l-1)+4]==v & ((gendata2[n,6*(l-1)+5]==j & gendata2[n,6*(l-1)+6]==v) | (gendata2[n,6*(l-1)+5]==v & gendata2[n,6*(l-1)+6]==j))){
				kuplet_fa_T=paste(kuplet_fa_T,j, sep=""); kuplet_mo_T=paste(kuplet_mo_T,v, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,i, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,u, sep="")}
				}
			}
		}
}}}}
}

if (length(setdiff(locus_on_X,l))!=length(locus_on_X) & gendata2$sex[n]==1){

if (gendata2[n,6*(l-1)+1]==1 & gendata2[n,6*(l-1)+2]==1 & gendata2[n,6*(l-1)+3]==1 & gendata2[n,6*(l-1)+4]==2 & gendata2[n,6*(l-1)+5]==2 & gendata2[n,6*(l-1)+6]==2){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,2, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,1, sep="")}

if (gendata2[n,6*(l-1)+1]==1 & gendata2[n,6*(l-1)+2]==1 & gendata2[n,6*(l-1)+3]==1 & gendata2[n,6*(l-1)+4]==2 & gendata2[n,6*(l-1)+5]==1 & gendata2[n,6*(l-1)+6]==1){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,1, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,2, sep="")}

if (gendata2[n,6*(l-1)+1]==2 & gendata2[n,6*(l-1)+2]==2 & gendata2[n,6*(l-1)+3]==1 & gendata2[n,6*(l-1)+4]==1 & gendata2[n,6*(l-1)+5]==1 & gendata2[n,6*(l-1)+6]==1){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,1, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,1, sep="")}

if (gendata2[n,6*(l-1)+1]==2 & gendata2[n,6*(l-1)+2]==2 & gendata2[n,6*(l-1)+3]==1 & gendata2[n,6*(l-1)+4]==2 & gendata2[n,6*(l-1)+5]==2 & gendata2[n,6*(l-1)+6]==2){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,2, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,1, sep="")}

if (gendata2[n,6*(l-1)+1]==2 & gendata2[n,6*(l-1)+2]==2 & gendata2[n,6*(l-1)+3]==1 & gendata2[n,6*(l-1)+4]==2 & gendata2[n,6*(l-1)+5]==1 & gendata2[n,6*(l-1)+6]==1){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,1, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,2, sep="")}

if (gendata2[n,6*(l-1)+1]==1 & gendata2[n,6*(l-1)+2]==1 & gendata2[n,6*(l-1)+3]==2 & gendata2[n,6*(l-1)+4]==2 & gendata2[n,6*(l-1)+5]==2 & gendata2[n,6*(l-1)+6]==2){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,2, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,2, sep="")}

if (gendata2[n,6*(l-1)+1]==2 & gendata2[n,6*(l-1)+2]==2 & gendata2[n,6*(l-1)+3]==2 & gendata2[n,6*(l-1)+4]==2 & gendata2[n,6*(l-1)+5]==2 & gendata2[n,6*(l-1)+6]==2){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,2, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,2, sep="")}

if (gendata2[n,6*(l-1)+1]==1 & gendata2[n,6*(l-1)+2]==1 & gendata2[n,6*(l-1)+3]==1 & gendata2[n,6*(l-1)+4]==1 & gendata2[n,6*(l-1)+5]==1 & gendata2[n,6*(l-1)+6]==1){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,1, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,1, sep="")}
}
}

if (length(setdiff(kuplet,kuplet_fa_T))!=length(kuplet) & length(setdiff(kuplet,kuplet_fa_NT))!=length(kuplet) & length(setdiff(kuplet,kuplet_mo_T))!=length(kuplet) & length(setdiff(kuplet,kuplet_mo_NT))!=length(kuplet)){
transmat[kuplet_fa_T,kuplet_fa_NT] = transmat[kuplet_fa_T,kuplet_fa_NT] + gendata2$countw[n]*abs(gendata2$pfaid[n])/sum(abs(gendata2$pfaid)*gendata2$countw)*sum(gendata2$countw)
transmat[kuplet_mo_T,kuplet_mo_NT] = transmat[kuplet_mo_T,kuplet_mo_NT] + gendata2$countw[n]*abs(gendata2$pfaid[n])/sum(abs(gendata2$pfaid)*gendata2$countw)*sum(gendata2$countw)}
gendata2=gendata2[,-ncol(gendata2)]

}
rm(i,j,l,n,u,v,kuplet_fa_NT,kuplet_fa_T,kuplet_mo_NT,kuplet_mo_T)

transmat[transmat==0]=0.0000001
sum(transmat)

# ------------------------------------------------------ #
# - TRANSMISSION INTENSITY OF ALLELES AT SINGLE LOCUS -- #
# ------------------------------------------------------ #

alpha=matrix(NA,nbloci,nballeles)
rownames(alpha)=c(paste("locus",1:nbloci,sep=""))
colnames(alpha)=c(paste("allele",1:nballeles,sep=""))

for (l in 1:nbloci){
for (a in 1:nballeles){
	alpha[l,a]=sum(transmat[substr(rownames(transmat),l,l)==paste(a),substr(colnames(transmat),l,l)!=paste(a)])/(sum(transmat[substr(rownames(transmat),l,l)==paste(a),substr(colnames(transmat),l,l)!=paste(a)])+sum(transmat[substr(rownames(transmat),l,l)!=paste(a),substr(colnames(transmat),l,l)==paste(a)]))
}}
rm(a,l)

# ------------------------------------------------------ #
# ----- NUMBER OF TRANSMITTED AND NOT-TRANSMITTED ------ #
# ------------------------------------------------------ #

k=0
nT=0
nNT=0
for (i in 1:(dim(transmat)[2]-1)){
for (j in (i+1):dim(transmat)[1]){
k=k+1
nT[k]=transmat[i,j]
nNT[k]=transmat[j,i]
}
}
rm(i,j,k)


# ----------------------------------------------------------- #
# ----- The Number of possible alternative hypotheses   ----- #
# ----------------------------------------------------------- #
fact=function(m){
 fm=1
 while (m>=2){
  fm=fm*m
  m=m-1}
 return(fm)}
# ----------------------------------------------------------- #
comb=function(n,p){
 while (n>=p){
  return(fact(n)/(fact(p)*fact(n-p)))}}
# ----------------------------------------------------------- #
nbmodel= nbloci
if(nbloci>=2){
for (i in 2:nbloci){
 nbmodel= nbmodel+ 2*comb(nbloci,i)
}}
nbmodel

rm(comb,fact,i)

# ------------------------------------------------------ #
# ------- SINGLE TRANSMISSION PROBABILITIES ------------ #
# ------------------------------------------------------ #

tau=matrix(0,nbmodel,length(nT))

ddl=0
nb_info_fam=NULL
for (m in 1:nbloci){
ddl[m]=length(taballeles[m,])-1

nb_info_fam[m] = (
nrow(gendata2[
gendata2[,paste("fal",m,"a",1, sep="")]!=
gendata2[,paste("fal",m,"a",2, sep="")]
,]) +
nrow(gendata2[
gendata2[,paste("mol",m,"a",1, sep="")]!=
gendata2[,paste("mol",m,"a",2, sep="")]
,]))

k=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	a=as.numeric(substr(rownames(transmat)[i],m,m))
	b=as.numeric(substr(colnames(transmat)[j],m,m))
	tau[m,k]=alpha[m,a]/(alpha[m,a]+alpha[m,b])
}}
}
rm(a,b,i,j,k)

# ------------------------------------------------------ #
# ------- 2-UPLET TRANSMISSION PROBABILITIES ----------- #
# ------------------------------------------------------ #

# MULTIPLICATIVE #

for (lm in 1:(nbloci-1)){
for (ln in (lm+1):nbloci){
m=m+1
ddl[m]=length(taballeles[lm,])-1 + length(taballeles[ln,])-1

nb_info_fam[m] = (
nrow(gendata2[
paste(gendata2[,paste("fal",lm,"a",1, sep="")],gendata2[,paste("fal",ln,"a",1, sep="")],sep="")!=
paste(gendata2[,paste("fal",lm,"a",2, sep="")],gendata2[,paste("fal",ln,"a",2, sep="")],sep="")
,]) +
nrow(gendata2[
paste(gendata2[,paste("mol",lm,"a",1, sep="")],gendata2[,paste("mol",ln,"a",1, sep="")],sep="")!=
paste(gendata2[,paste("mol",lm,"a",2, sep="")],gendata2[,paste("mol",ln,"a",2, sep="")],sep="")
,]))

k=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	a=as.numeric(substr(rownames(transmat)[i],lm,lm))
	b=as.numeric(substr(colnames(transmat)[j],lm,lm))
	c=as.numeric(substr(rownames(transmat)[i],ln,ln))
	d=as.numeric(substr(colnames(transmat)[j],ln,ln))
	tau[m,k]=alpha[lm,a]*alpha[ln,c]/(alpha[lm,a]*alpha[ln,c]+alpha[lm,b]*alpha[ln,d])
}}
}}
rm(a,b,c,d,i,j,k,lm,ln)

# EPISTASIS #

for (lm in 1:(nbloci-1)){
for (ln in (lm+1):nbloci){
m=m+1
ddl[m]=length(taballeles[lm,])*length(taballeles[ln,])-1

nb_info_fam[m] = (
nrow(gendata2[
paste(gendata2[,paste("fal",lm,"a",1, sep="")],gendata2[,paste("fal",ln,"a",1, sep="")],sep="")!=
paste(gendata2[,paste("fal",lm,"a",2, sep="")],gendata2[,paste("fal",ln,"a",2, sep="")],sep="")
,]) +
nrow(gendata2[
paste(gendata2[,paste("mol",lm,"a",1, sep="")],gendata2[,paste("mol",ln,"a",1, sep="")],sep="")!=
paste(gendata2[,paste("mol",lm,"a",2, sep="")],gendata2[,paste("mol",ln,"a",2, sep="")],sep="")
,]))


k=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	a=as.numeric(substr(rownames(transmat)[i],lm,lm))
	b=as.numeric(substr(colnames(transmat)[j],lm,lm))
	c=as.numeric(substr(rownames(transmat)[i],ln,ln))
	d=as.numeric(substr(colnames(transmat)[j],ln,ln))

	x=sum(transmat[substr(rownames(transmat),lm,lm)==paste(a) & substr(rownames(transmat),ln,ln)==paste(c),
			   substr(colnames(transmat),lm,lm)!=paste(a) | substr(colnames(transmat),ln,ln)!=paste(c)])/
	 (sum(transmat[substr(rownames(transmat),lm,lm)==paste(a) & substr(rownames(transmat),ln,ln)==paste(c),
			   substr(colnames(transmat),lm,lm)!=paste(a) | substr(colnames(transmat),ln,ln)!=paste(c)])
	 +sum(transmat[substr(rownames(transmat),lm,lm)!=paste(a) | substr(rownames(transmat),ln,ln)!=paste(c),
			   substr(colnames(transmat),lm,lm)==paste(a) & substr(colnames(transmat),ln,ln)==paste(c)]))

	y=sum(transmat[substr(rownames(transmat),lm,lm)==paste(b) & substr(rownames(transmat),ln,ln)==paste(d),
			   substr(colnames(transmat),lm,lm)!=paste(b) | substr(colnames(transmat),ln,ln)!=paste(d)])/
	 (sum(transmat[substr(rownames(transmat),lm,lm)==paste(b) & substr(rownames(transmat),ln,ln)==paste(d),
			   substr(colnames(transmat),lm,lm)!=paste(b) | substr(colnames(transmat),ln,ln)!=paste(d)])
	 +sum(transmat[substr(rownames(transmat),lm,lm)!=paste(b) | substr(rownames(transmat),ln,ln)!=paste(d),
			   substr(colnames(transmat),lm,lm)==paste(b) & substr(colnames(transmat),ln,ln)==paste(d)]))

	tau[m,k]=x/(x+y)
}}
}}
rm(a,b,c,d,i,j,k,x,y,lm,ln)

# ------------------------------------------------------ #
# ------- L-UPLET TRANSMISSION PROBABILITIES ----------- #
# ------------------------------------------------------ #

# MULTIPLICATIVE #

m=m+1
k=a=b=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	x=y=1
	df=0
	for (l in 1:nbloci){
	a[l]=as.numeric(substr(rownames(transmat)[i],l,l))
	b[l]=as.numeric(substr(colnames(transmat)[j],l,l))
	x=x*alpha[l,a[l]]
	y=y*alpha[l,b[l]]
	df = df + length(taballeles[l,])-1
				}
	tau[m,k]=x/(x+y)
	ddl[m]=df
	nb_info_fam[m]=sum(nT+nNT)
}}
rm(a,b,i,j,k,l,x,y,df)

# EPISTASIS #

m=m+1
k=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	tau[m,k]=(sum(transmat[i,-i])/(sum(transmat[i,-i])+sum(transmat[-i,i])))/((sum(transmat[i,-i])/(sum(transmat[i,-i])+sum(transmat[-i,i])))+(sum(transmat[j,-j])/(sum(transmat[j,-j])+sum(transmat[-j,j]))))
}}
rm(i,j,k)

df=1
for (l in 1:nbloci){df=df*length(taballeles[l,])}
ddl[m]=df-1
nb_info_fam[m]=sum(nT+nNT)
rm(l,df)

# ------------------------------------------------------ #
# ------------------ LOG-LIKELIHOODS ------------------- #
# ------------------------------------------------------ #

# MODEL0: WHITE MODEL

LL0= -log(2)*sum(nT+nNT)
LL0

# MODEL 1 to m

LL=0
ETDT=0
pvalETDT=0

for (l in 1:m){
	LL[l]= sum(nT*log(tau[l,]/(1-tau[l,]))) + sum((nT+nNT)*log(1-tau[l,]))
	ETDT[l]=2*(LL[l]-LL0)
	pvalETDT[l]=1-pchisq(ETDT[l],ddl[l])
}
rm(l)

# ------------------------------------------------------ #
# ---------------- TO DISPLAY RESULTS ------------------ #
# ------------------------------------------------------ #

models=c(L1,L2,L3,

paste("Multplicative",L1,L2,sep="*"),paste("Multplicative",L1,L3,sep="*"),paste("Multplicative",L2,L3,sep="*"),
paste("Epistasis",L1,L2,sep="*"),paste("Epistasis",L1,L3,sep="*"),paste("Epistasis",L2,L3,sep="*"),
paste("Multplicative",L1,L2,L3,sep="*"), paste("Epistasis",L1,L2,L3,sep="*"))

ResultsMultilocus=data.frame(data.frame(models,nb_info_fam=round(nb_info_fam,1),LL1=round(LL,1),LL0=round(LL0,1),ETDT=round(ETDT,2),ddl,pvalETDT))

# ------------------------------------------------------ #

ResultsMultilocus$ord = 1:nrow(ResultsMultilocus)
ResultsMultilocus=ResultsMultilocus[order(ResultsMultilocus[,"pvalETDT"]),]
p = ResultsMultilocus[,"pvalETDT"]
mtests=length(p)
q=0
for (j in 1:mtests){
	q[j]=min(mtests*p[j:mtests]/(j:mtests))
}
ResultsMultilocus$FDR = q
ResultsMultilocus=ResultsMultilocus[order(ResultsMultilocus[,"ord"]),]
ResultsMultilocus=ResultsMultilocus[,-(ncol(ResultsMultilocus)-1)]

rm(p,q,j)

# ------------------------------------------------------ #

write.csv2(ResultsMultilocus, file = paste(path,"/MultilocusTDT_Output/Results/","weighted_res_multilocus_",gene,"_phen",phenotype,"_","village",paste(village1,village2,sep=""),".csv",sep=""), quote=F, row.names=F)
write.table(round(transmat,1), file = paste(path,"/MultilocusTDT_Output/Results/","weighted_transmat_",gene,"_phen",phenotype,"_","village",paste(village1,village2,sep=""),".txt",sep=""), sep="\t", quote=F)

######################################################
######################################################
######################################################
######################################################
######################################################

SimTestStat = matrix(NA,nbsimul,nrow(ResultsMultilocus))

for (sim in 1:nbsimul){

	gendata2simulated = gendata2saved

	for (n in 1:nrow(gendata2simulated)){
	for (l in 1:nbloci){
		faa1a2 = sample(gendata2simulated[n,(6*(l-1)+1):(6*(l-1)+2)],1)
		moa1a2 = sample(gendata2simulated[n,(6*(l-1)+3):(6*(l-1)+4)],1)
		gendata2simulated[n,6*(l-1)+5] = min(faa1a2,moa1a2)
		gendata2simulated[n,6*(l-1)+6] = max(faa1a2,moa1a2)
	}}
	rm(n,l)

######################################################

# ----------------------------------------------------------- #
# SIMULATION OF POSSIBLE CHILDREN FOR AMBIGUOUS TRANSMISSIONS #
# ----------------------------------------------------------- #
gendata2simulated$countw=1
gendata2simulated$realchild=1

for (n in 1:nrow(gendata2simulated)){
nbdoubt=0
locusdoubt=0
for (l in 1:nbloci){
	if (gendata2simulated[n,6*(l-1)+1]==gendata2simulated[n,6*(l-1)+3] & gendata2simulated[n,6*(l-1)+1]==gendata2simulated[n,6*(l-1)+5] &
	    gendata2simulated[n,6*(l-1)+2]==gendata2simulated[n,6*(l-1)+4] & gendata2simulated[n,6*(l-1)+2]==gendata2simulated[n,6*(l-1)+6] &
	    gendata2simulated[n,6*(l-1)+1]!=gendata2simulated[n,6*(l-1)+2] &
	    gendata2simulated$realchild[n]==1){
	nbdoubt=nbdoubt+1
	locusdoubt[nbdoubt]=l
	}
}

if (nbdoubt>0) {
gensimchild=NULL
for (p in (nbdoubt-1):0) {
gensimchild = c(gensimchild,rep(c(rep(1,2^p),rep(2,2^p)),2^(nbdoubt-p-1)))
}
gensimchild=matrix(gensimchild,2^nbdoubt,nbdoubt)

for (i in 1:2^nbdoubt){
gendata2simulated=rbind(gendata2simulated,gendata2simulated[n,])
gendata2simulated[nrow(gendata2simulated),6*(locusdoubt-1)+5]=gensimchild[i,]
gendata2simulated[nrow(gendata2simulated),6*(locusdoubt-1)+6]=gensimchild[i,]
}

gendata2simulated$countw[n]=0
gendata2simulated$countw[(nrow(gendata2simulated)-2^nbdoubt+1):nrow(gendata2simulated)]=1/2^nbdoubt
gendata2simulated$realchild[(nrow(gendata2simulated)-2^nbdoubt+1):nrow(gendata2simulated)]=0
}}
rm(n,l,p,i)

# ------------------------------------------------------ #
# --- TO COMPUTE MATRIX OF SUMULTANEOUS TRANSMISSION --- #
# ------------------------------------------------------ #
transmat=matrix(0,length(kuplet),length(kuplet))
rownames(transmat)=kuplet
colnames(transmat)=kuplet

for (n in 1:nrow(gendata2simulated)){
gendata2simulated$idDN=gendata2simulated$child
gendata2simulated=merge(gendata2simulated, ped[,c("idDN","pfaid")], by=c("idDN"), all.x=T)
gendata2simulated=gendata2simulated[,-1]

kuplet_fa_T=NULL; kuplet_mo_T=NULL; kuplet_fa_NT=NULL; kuplet_mo_NT=NULL

for (l in 1:nbloci){

if ((length(setdiff(locus_on_X,l))==length(locus_on_X)) | (length(setdiff(locus_on_X,l))!=length(locus_on_X) & gendata2simulated$sex[n]==2)){

for (i in taballeles[l,][is.na(taballeles[l,])==FALSE]){
for (j in i:max(taballeles[l,][is.na(taballeles[l,])==FALSE])){
for (u in taballeles[l,][is.na(taballeles[l,])==FALSE]){
for (v in u:max(taballeles[l,][is.na(taballeles[l,])==FALSE])){
	if (gendata2simulated[n,6*(l-1)+1]==i & gendata2simulated[n,6*(l-1)+2]==j & gendata2simulated[n,6*(l-1)+3]==u & gendata2simulated[n,6*(l-1)+4]==v & ((gendata2simulated[n,6*(l-1)+5]==i & gendata2simulated[n,6*(l-1)+6]==u) | (gendata2simulated[n,6*(l-1)+5]==u & gendata2simulated[n,6*(l-1)+6]==i))){
	kuplet_fa_T=paste(kuplet_fa_T,i, sep=""); kuplet_mo_T=paste(kuplet_mo_T,u, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,j, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,v, sep="")}
	else {
		if (gendata2simulated[n,6*(l-1)+1]==i & gendata2simulated[n,6*(l-1)+2]==j & gendata2simulated[n,6*(l-1)+3]==u & gendata2simulated[n,6*(l-1)+4]==v & ((gendata2simulated[n,6*(l-1)+5]==i & gendata2simulated[n,6*(l-1)+6]==v) | (gendata2simulated[n,6*(l-1)+5]==v & gendata2simulated[n,6*(l-1)+6]==i))){
		kuplet_fa_T=paste(kuplet_fa_T,i, sep=""); kuplet_mo_T=paste(kuplet_mo_T,v, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,j, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,u, sep="")}
		else {
			if (gendata2simulated[n,6*(l-1)+1]==i & gendata2simulated[n,6*(l-1)+2]==j & gendata2simulated[n,6*(l-1)+3]==u & gendata2simulated[n,6*(l-1)+4]==v & ((gendata2simulated[n,6*(l-1)+5]==j & gendata2simulated[n,6*(l-1)+6]==u) | (gendata2simulated[n,6*(l-1)+5]==u & gendata2simulated[n,6*(l-1)+6]==j))){
			kuplet_fa_T=paste(kuplet_fa_T,j, sep=""); kuplet_mo_T=paste(kuplet_mo_T,u, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,i, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,v, sep="")}
			else {
				if (gendata2simulated[n,6*(l-1)+1]==i & gendata2simulated[n,6*(l-1)+2]==j & gendata2simulated[n,6*(l-1)+3]==u & gendata2simulated[n,6*(l-1)+4]==v & ((gendata2simulated[n,6*(l-1)+5]==j & gendata2simulated[n,6*(l-1)+6]==v) | (gendata2simulated[n,6*(l-1)+5]==v & gendata2simulated[n,6*(l-1)+6]==j))){
				kuplet_fa_T=paste(kuplet_fa_T,j, sep=""); kuplet_mo_T=paste(kuplet_mo_T,v, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,i, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,u, sep="")}
				}
			}
		}
}}}}
}

if (length(setdiff(locus_on_X,l))!=length(locus_on_X) & gendata2simulated$sex[n]==1){

if (gendata2simulated[n,6*(l-1)+1]==1 & gendata2simulated[n,6*(l-1)+2]==1 & gendata2simulated[n,6*(l-1)+3]==1 & gendata2simulated[n,6*(l-1)+4]==2 & gendata2simulated[n,6*(l-1)+5]==2 & gendata2simulated[n,6*(l-1)+6]==2){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,2, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,1, sep="")}

if (gendata2simulated[n,6*(l-1)+1]==1 & gendata2simulated[n,6*(l-1)+2]==1 & gendata2simulated[n,6*(l-1)+3]==1 & gendata2simulated[n,6*(l-1)+4]==2 & gendata2simulated[n,6*(l-1)+5]==1 & gendata2simulated[n,6*(l-1)+6]==1){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,1, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,2, sep="")}

if (gendata2simulated[n,6*(l-1)+1]==2 & gendata2simulated[n,6*(l-1)+2]==2 & gendata2simulated[n,6*(l-1)+3]==1 & gendata2simulated[n,6*(l-1)+4]==1 & gendata2simulated[n,6*(l-1)+5]==1 & gendata2simulated[n,6*(l-1)+6]==1){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,1, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,1, sep="")}

if (gendata2simulated[n,6*(l-1)+1]==2 & gendata2simulated[n,6*(l-1)+2]==2 & gendata2simulated[n,6*(l-1)+3]==1 & gendata2simulated[n,6*(l-1)+4]==2 & gendata2simulated[n,6*(l-1)+5]==2 & gendata2simulated[n,6*(l-1)+6]==2){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,2, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,1, sep="")}

if (gendata2simulated[n,6*(l-1)+1]==2 & gendata2simulated[n,6*(l-1)+2]==2 & gendata2simulated[n,6*(l-1)+3]==1 & gendata2simulated[n,6*(l-1)+4]==2 & gendata2simulated[n,6*(l-1)+5]==1 & gendata2simulated[n,6*(l-1)+6]==1){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,1, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,2, sep="")}

if (gendata2simulated[n,6*(l-1)+1]==1 & gendata2simulated[n,6*(l-1)+2]==1 & gendata2simulated[n,6*(l-1)+3]==2 & gendata2simulated[n,6*(l-1)+4]==2 & gendata2simulated[n,6*(l-1)+5]==2 & gendata2simulated[n,6*(l-1)+6]==2){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,2, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,2, sep="")}

if (gendata2simulated[n,6*(l-1)+1]==2 & gendata2simulated[n,6*(l-1)+2]==2 & gendata2simulated[n,6*(l-1)+3]==2 & gendata2simulated[n,6*(l-1)+4]==2 & gendata2simulated[n,6*(l-1)+5]==2 & gendata2simulated[n,6*(l-1)+6]==2){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,2, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,2, sep="")}

if (gendata2simulated[n,6*(l-1)+1]==1 & gendata2simulated[n,6*(l-1)+2]==1 & gendata2simulated[n,6*(l-1)+3]==1 & gendata2simulated[n,6*(l-1)+4]==1 & gendata2simulated[n,6*(l-1)+5]==1 & gendata2simulated[n,6*(l-1)+6]==1){
	kuplet_fa_T=paste(kuplet_fa_T,1, sep=""); kuplet_mo_T=paste(kuplet_mo_T,1, sep=""); kuplet_fa_NT=paste(kuplet_fa_NT,1, sep=""); kuplet_mo_NT=paste(kuplet_mo_NT,1, sep="")}
}}

if (length(setdiff(kuplet,kuplet_fa_T))!=length(kuplet) & length(setdiff(kuplet,kuplet_fa_NT))!=length(kuplet) & length(setdiff(kuplet,kuplet_mo_T))!=length(kuplet) & length(setdiff(kuplet,kuplet_mo_NT))!=length(kuplet)){
transmat[kuplet_fa_T,kuplet_fa_NT] = transmat[kuplet_fa_T,kuplet_fa_NT] + gendata2simulated$countw[n]*abs(gendata2simulated$pfaid[n])/sum(abs(gendata2simulated$pfaid)*gendata2simulated$countw)*sum(gendata2simulated$countw)
transmat[kuplet_mo_T,kuplet_mo_NT] = transmat[kuplet_mo_T,kuplet_mo_NT] + gendata2simulated$countw[n]*abs(gendata2simulated$pfaid[n])/sum(abs(gendata2simulated$pfaid)*gendata2simulated$countw)*sum(gendata2simulated$countw)}
gendata2simulated=gendata2simulated[,-ncol(gendata2simulated)]
}
rm(i,j,l,n,u,v,kuplet_fa_NT,kuplet_fa_T,kuplet_mo_NT,kuplet_mo_T)

transmat[transmat==0]=0.0000001

# ------------------------------------------------------ #
# - TRANSMISSION INTENSITY OF ALLELES AT SINGLE LOCUS -- #
# ------------------------------------------------------ #
alpha=matrix(NA,nbloci,nballeles)
rownames(alpha)=c(paste("locus",1:nbloci,sep=""))
colnames(alpha)=c(paste("allele",1:nballeles,sep=""))

for (l in 1:nbloci){
for (a in 1:nballeles){
	alpha[l,a]=sum(transmat[substr(rownames(transmat),l,l)==paste(a),substr(colnames(transmat),l,l)!=paste(a)])/(sum(transmat[substr(rownames(transmat),l,l)==paste(a),substr(colnames(transmat),l,l)!=paste(a)])+sum(transmat[substr(rownames(transmat),l,l)!=paste(a),substr(colnames(transmat),l,l)==paste(a)]))
}}
rm(a,l)

# ------------------------------------------------------ #
# ----- NUMBER OF TRANSMITTED AND NOT-TRANSMITTED ------ #
# ------------------------------------------------------ #
k=0
nT=0
nNT=0
for (i in 1:(dim(transmat)[2]-1)){
for (j in (i+1):dim(transmat)[1]){
k=k+1
nT[k]=transmat[i,j]
nNT[k]=transmat[j,i]
}}
rm(i,j,k)

# ----------------------------------------------------------- #
# ----- The Number of possible alternative hypotheses   ----- #
# ----------------------------------------------------------- #
fact=function(m){
 fm=1
 while (m>=2){
  fm=fm*m
  m=m-1}
 return(fm)}
# ----------------------------------------------------------- #
comb=function(n,p){
 while (n>=p){
  return(fact(n)/(fact(p)*fact(n-p)))}}
# ----------------------------------------------------------- #
nbmodel= nbloci
if(nbloci>=2){
for (i in 2:nbloci){
 nbmodel= nbmodel+ 2*comb(nbloci,i)
}}
nbmodel
rm(comb,fact,i)

# ------------------------------------------------------ #
# ------- SINGLE TRANSMISSION PROBABILITIES ------------ #
# ------------------------------------------------------ #
tau=matrix(0,nbmodel,length(nT))
ddl=0
nb_info_fam=NULL
for (m in 1:nbloci){
ddl[m]=length(taballeles[m,])-1

k=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	a=as.numeric(substr(rownames(transmat)[i],m,m))
	b=as.numeric(substr(colnames(transmat)[j],m,m))
	tau[m,k]=alpha[m,a]/(alpha[m,a]+alpha[m,b])
}}
}
rm(a,b,i,j,k)

# ------------------------------------------------------ #
# ------- 2-UPLET TRANSMISSION PROBABILITIES ----------- #
# ------------------------------------------------------ #
# MULTIPLICATIVE #

for (lm in 1:(nbloci-1)){
for (ln in (lm+1):nbloci){
m=m+1
ddl[m]=length(taballeles[lm,])-1 + length(taballeles[ln,])-1

k=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	a=as.numeric(substr(rownames(transmat)[i],lm,lm))
	b=as.numeric(substr(colnames(transmat)[j],lm,lm))
	c=as.numeric(substr(rownames(transmat)[i],ln,ln))
	d=as.numeric(substr(colnames(transmat)[j],ln,ln))
	tau[m,k]=alpha[lm,a]*alpha[ln,c]/(alpha[lm,a]*alpha[ln,c]+alpha[lm,b]*alpha[ln,d])
}}
}}
rm(a,b,c,d,i,j,k,lm,ln)

# EPISTASIS #

for (lm in 1:(nbloci-1)){
for (ln in (lm+1):nbloci){
m=m+1
ddl[m]=length(taballeles[lm,])*length(taballeles[ln,])-1

k=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	a=as.numeric(substr(rownames(transmat)[i],lm,lm))
	b=as.numeric(substr(colnames(transmat)[j],lm,lm))
	c=as.numeric(substr(rownames(transmat)[i],ln,ln))
	d=as.numeric(substr(colnames(transmat)[j],ln,ln))

	x=sum(transmat[substr(rownames(transmat),lm,lm)==paste(a) & substr(rownames(transmat),ln,ln)==paste(c),
			   substr(colnames(transmat),lm,lm)!=paste(a) | substr(colnames(transmat),ln,ln)!=paste(c)])/
	 (sum(transmat[substr(rownames(transmat),lm,lm)==paste(a) & substr(rownames(transmat),ln,ln)==paste(c),
			   substr(colnames(transmat),lm,lm)!=paste(a) | substr(colnames(transmat),ln,ln)!=paste(c)])
	 +sum(transmat[substr(rownames(transmat),lm,lm)!=paste(a) | substr(rownames(transmat),ln,ln)!=paste(c),
			   substr(colnames(transmat),lm,lm)==paste(a) & substr(colnames(transmat),ln,ln)==paste(c)]))

	y=sum(transmat[substr(rownames(transmat),lm,lm)==paste(b) & substr(rownames(transmat),ln,ln)==paste(d),
			   substr(colnames(transmat),lm,lm)!=paste(b) | substr(colnames(transmat),ln,ln)!=paste(d)])/
	 (sum(transmat[substr(rownames(transmat),lm,lm)==paste(b) & substr(rownames(transmat),ln,ln)==paste(d),
			   substr(colnames(transmat),lm,lm)!=paste(b) | substr(colnames(transmat),ln,ln)!=paste(d)])
	 +sum(transmat[substr(rownames(transmat),lm,lm)!=paste(b) | substr(rownames(transmat),ln,ln)!=paste(d),
			   substr(colnames(transmat),lm,lm)==paste(b) & substr(colnames(transmat),ln,ln)==paste(d)]))
	tau[m,k]=x/(x+y)
}}
}}
rm(a,b,c,d,i,j,k,x,y,lm,ln)

# ------------------------------------------------------ #
# ------- L-UPLET TRANSMISSION PROBABILITIES ----------- #
# ------------------------------------------------------ #
# MULTIPLICATIVE #

m=m+1
k=a=b=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	x=y=1
	df=0
	for (l in 1:nbloci){
	a[l]=as.numeric(substr(rownames(transmat)[i],l,l))
	b[l]=as.numeric(substr(colnames(transmat)[j],l,l))
	x=x*alpha[l,a[l]]
	y=y*alpha[l,b[l]]
	df = df + length(taballeles[l,])-1
				}
	tau[m,k]=x/(x+y)
	ddl[m]=df
	nb_info_fam[m]=sum(nT+nNT)
}}
rm(a,b,i,j,k,l,x,y,df)

# EPISTASIS #

m=m+1
k=0
for (i in 1:(nbkuplet-1)){
for (j in (i+1):nbkuplet){
	k=k+1
	tau[m,k]=(sum(transmat[i,-i])/(sum(transmat[i,-i])+sum(transmat[-i,i])))/((sum(transmat[i,-i])/(sum(transmat[i,-i])+sum(transmat[-i,i])))+(sum(transmat[j,-j])/(sum(transmat[j,-j])+sum(transmat[-j,j]))))
}}
rm(i,j,k)

df=1
for (l in 1:nbloci){df=df*length(taballeles[l,])}
ddl[m]=df-1
rm(l,df)

# ------------------------------------------------------ #
# ------------------ LOG-LIKELIHOODS ------------------- #
# ------------------------------------------------------ #
# MODEL0: WHITE MODEL

LL0= -log(2)*sum(nT+nNT)
LL0

# MODEL 1 to m
LL=0
ETDT=0
pvalETDT=0

for (l in 1:m){
	LL[l]= sum(nT*log(tau[l,]/(1-tau[l,]))) + sum((nT+nNT)*log(1-tau[l,]))
	ETDT[l]=2*(LL[l]-LL0)
	pvalETDT[l]=1-pchisq(ETDT[l],ddl[l])
}
rm(l)

######################################################

	SimTestStat[sim,] = ETDT
}
write.table(SimTestStat, file = paste(path,"/MultilocusTDT_Output/Results/","weighted_SimTestStat_multilocus_",gene,"_phen",phenotype,"_","village",paste(village1,village2,sep=""),".txt",sep=""), quote=F, row.names=F, col.names=F)

# ------------------------------------------------------ #
empPval = NULL
for (i in 1:ncol(SimTestStat)){
	empPval[i] = sum(ifelse(SimTestStat[,i]>ETDT[i],1,0)) / nbsimul 
}
ResultsMultilocus$empPval = empPval

# ------------------------------------------------------ #
ResultsMultilocus$ord = 1:nrow(ResultsMultilocus)
ResultsMultilocus=ResultsMultilocus[order(ResultsMultilocus[,"empPval"]),]
p = ResultsMultilocus[,"empPval"]
mtests=length(p)
q=0
for (j in 1:mtests){
	q[j]=min(mtests*p[j:mtests]/(j:mtests))
}
ResultsMultilocus$FDRemp = q
ResultsMultilocus=ResultsMultilocus[order(ResultsMultilocus[,"ord"]),]
ResultsMultilocus=ResultsMultilocus[,-(ncol(ResultsMultilocus)-1)]

rm(p,q,j)

# ------------------------------------------------------ #
write.csv2(ResultsMultilocus, file = paste(path,"/MultilocusTDT_Output/Results/","weighted_res_multilocus_",gene,"_phen",phenotype,"_","village",paste(village1,village2,sep=""),".csv",sep=""), quote=F, row.names=F)

######################################################
######################################################
######################################################
######################################################
######################################################

}}}




