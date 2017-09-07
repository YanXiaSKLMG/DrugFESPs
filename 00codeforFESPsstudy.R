#R code for methylation study 
#Title:Risperidone-induced DNA methylation alterations in first-episode drug-naïve schizophrenia patients
#      and their relationship with neurophysiological phenotypes
#Code author: Yan Xia
#email:yanxia.sklmg@gmail.com

##
setwd("/zs32/home/yxia/04methyZong/zongidat")
baseDir <- "/zs32/home/yxia/04methyZong"

## load packages
library(minfi)
library(ChAMP)# version 2.0.1
library(WGCNA)
library(DMRcate)
library(coin)

## preprocessing 
###01 load data and probes filtering
myload<-champ.load(directory = getwd(), 
	                 methValue = "B", resultsDir = paste(baseDir,"resultsChamp", sep = "/"), 
	                 filterXY = FALSE, 
	                 QCimages = TRUE, 
	                 filterDetP = TRUE,
                     detPcut = 0.01, 
                     remeDetP = 0, 
                     filterBeads=TRUE, 
                     beadCutoff=0.05, 
                     filterNoCG=FALSE,
                     filterSNPs=TRUE,
                     filterMultiHit=TRUE)

###02 remove type II bias by BMIQ
Norm<-champ.norm(beta = myload$beta, 
	             rgSet = myload$rgSet, 
	             pd = myload$pd, 
	             mset = myload$mset,
                 sampleSheet = "samplesheet.csv", 
                 resultsDir = paste(getwd(), "resultsChamp",sep = "/"), 
                 methValue = "B", 
                 fromIDAT = TRUE, 
                 norm = "BMIQ", 
                 fromFile = FALSE, betaFile,
                 filter = TRUE, 
                 filterXY = FALSE, 
                 QCimages = FALSE, 
                 plotBMIQ = TRUE)

###03 estimate the confounders by SVD
champ.SVD(beta = Norm$beta, 
	      rgSet = myload$rgSet, 
	      detP = myload$detP, 
	      pd = myload$pd,
          loadFile = FALSE, 
          methProfile = FALSE,
          controlProfile = FALSE, 
          studyInfo = FALSE, 
          infoFactor = c(), 
          resultsDir = paste(getwd(), "resultsChamp", sep = "/"))

###04 remove the batch effect by Combat
Cbeta<-champ.runCombat(beta.c = Norm$beta, 
	                   pd = myload$pd, 
	                   logitTrans = TRUE)

###05 cell type composion estimate
Cell<-champ.refbase(beta=Cbeta$beta,
	                arraytype="450K")
names(Cell)

###06 remove the confounders by linear regression
mmage=model.matrix(~-1+myload$pd$Age+myload$pd$Sex+myload$pd$Smoke+myload$pd$Drink)
z=apply(Cell$CorrectedBeta,1,function(x){residuals(lm(x~mmage))})
Lmbeta<-t(z)+rowMeans(Cell$CorrectedBeta)

###07 save the data
save(myload,Norm,Cbeta,Cell,Lmbeta,file="20160713myloadNormCbetaCellLmbeta.RData")

###08 filter probes identified in Nordlund et al.



##Identification of differentially methylated individual CpG sites and genomic regions
##00 normal distribution test


###01 DMP between FESPs and Controls
load(file="/zs32/home/yxia/04methyZong/20161215datatraits/BefO_AftO_CtlO.RData")
group<-factor(c(rep("Bef",38),rep("Ctl",38)))
befctl<-cbind(BefO,CtlO)
BefCtl_tVal=rep(10,nrow(befctl))
BefCtl_tPval=rep(10,nrow(befctl))
BefMean=rep(10,nrow(befctl))
CtlMean=rep(10,nrow(befctl))
BefCtl_permutationP=rep(10,nrow(befctl))
for (i in 1:nrow(befctl)) 
{   ttest=t.test(befctl[i,]~group,data=rbind(group,befctl[1,]),var.equal=TRUE)
	BefCtl_tVal[i]=ttest$statistic; 
	BefCtl_tPval[i]=ttest$p.value;
	BefMean[i]=ttest$estimate[1];
	CtlMean[i]=ttest$estimate[2];
	BefCtl_permutationP[i]=pvalue(oneway_test(befctl[i,]~group,data=data.frame(group,befctl[i,])))
}
result<-data.frame(probeID=rownames(befctl),
	               BefCtl_tVal=BefCtl_tVal,
	               BefCtl_tPval=BefCtl_tPval,
	               BefCtl_permutationP=BefCtl_permutationP,
	               BefMean=BefMean,
	               CtlMean=CtlMean)
save(result,"20170505permuationresults.RData")

###02 DMP between pre- and post-treatment samples
#### paired t test
load("20160713myloadNormCbetaCellLmbeta.RData")
bef<-subset(myload$pd,Sample_Group == "BEF")
aft<-subset(myload$pd,Sample_Group == "AFT")
befbeta<-Lmbeta[,bef$Sample_Name]
aftbeta<-Lmbeta[,aft$Sample_Name]
data<-cbind(befbeta,aftbeta)
dataO<-data[,order(colnames(data))]
 a<-seq(from=1,to=75,by=2)
 b<-seq(from=2,to=76,by=2)
tval<-apply(dataO,1,function(x){t.test(x[a],x[b],paired=TRUE)$statistic})
pval<-apply(dataO,1,function(x){t.test(x[a],x[b],paired=TRUE)$p.value})
qval<-p.adjust(pval,method="BH")

BEFAFT<-data.frame(tval=tval,
	               pval=pval,
	               BHpval=qval)

#### permuation test ??
# BefCtl_permutationP[i]=pvalue(oneway_test(befctl[i,]~group,data=data.frame(group,befctl[i,])))

####DMR for pre- vs post- treatment
infor<-read.csv("Sampleinfor.csv",header=T)
dim(infor)
rownames(infor)<-infor[,1]
patient<-factor(infor$Ind)
Finfor<-subset(infor,Sample_Group!="CTL")
patient<-factor(Finfor$Ind)
type<-factor(Finfor$Sample_Group)
FLmbeta<-Lmbeta[,rownames(Finfor)]
dim(FLmbeta) #[1] 445283     76
design <- model.matrix(~patient + type) 
DMP <- cpg.annotate(FLmbeta,  
	                analysis.type="differential",
                    design=design, 
                    contrasts=FALSE, 
                    cont.matrix=NULL, 
                    coef=39)###
DMR <- dmrcate(DMP, lambda=1000, C=2,pcutoff=0.05)
names(DMR)
save(DMR,file="20160725DMRcateresults_DMR.RData")


## normatlization test
load(file="/zs32/home/yxia/04methyZong/20161215datatraits/BefO_AftO_CtlO.RData")
enableWGCNAThreads()
### unify_conn
Conn_Bef=softConnectivity(t(BefO),corFnc = "bicor",power=6) 
Conn_Aft=softConnectivity(t(AftO),corFnc = "bicor",power=6) 
Conn_Ctl=softConnectivity(t(CtlO),corFnc = "bicor",power=6)
Conn_diff0<-abs(Conn_Bef - Conn_Aft)
Conn_diff1<-abs(Conn_Aft - Conn_Ctl)
Csign=sign((Conn_Ctl -Conn_Bef)*(Conn_Aft - Conn_Bef))
Conn_diff<-Csign*Conn_diff0/(Conn_diff1+0.1)
unify_conn<-(Conn_diff-min(Conn_diff))/(max(Conn_diff)-min(Conn_diff))
### unify_t
#BefAft_tVal[i]=t.test(as.numeric(BefO[i,]),as.numeric(AftO[i,]),paired=TRUE)$statistic;
#AftCtl_tVal[i]=t.test(as.numeric(CtlO[i,]),as.numeric(AftO[i,]))$statistic; 
Tsign=sign((Ctl_Mean -Bef_Mean)*(Aft_Mean - Bef_Mean))
t_diff=Tsign*abs(BefAft_tVal)/(abs(AftCtl_tVal))
unify_t=(t_diff - min(t_diff))/(max(t_diff) - min(t_diff))
### weight
weight=(unify_conn + unify_t)/2
Conn=cbind(rownames(BefO),
	       Conn_case1,
	       Conn_case2,
	       Conn_ctrl,
	       Conn_diff,
	       unify_diff_Conn,
	       tVal,
	       tPval,
	       unify_diff,unify)
dim(Conn)
write.csv(Conn,file="H.ConnectivityDiff_Tvalue.csv")


## replication

##correlation with clinical phenotypes







