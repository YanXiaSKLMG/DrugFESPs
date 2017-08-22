
setwd("/zs32/home/yxia/04methyZong/20161215datatraits/20170413")
load(file="/zs32/home/yxia/04methyZong/20161215datatraits/BefO_AftO_CtlO.RData")
enableWGCNAThreads()
Conn_Bef=softConnectivity(t(BefO),corFnc = "bicor",power=6) 
Conn_Aft=softConnectivity(t(AftO),corFnc = "bicor",power=6) 
Conn_Ctl=softConnectivity(t(CtlO),corFnc = "bicor",power=6)
save(Conn_Bef,Conn_Aft,Conn_Ctl,file="Conn.RData")

Conn_diff0<-abs(Conn_Bef - Conn_Aft)
Conn_diff1<-abs(Conn_Aft - Conn_Ctl)

Csign=sign((Conn_Ctl -Conn_Bef)*(Conn_Aft - Conn_Bef))

Conn_diff<-Csign*Conn_diff0/(Conn_diff1+0.1)

unify_conn<-(Conn_diff-min(Conn_diff))/(max(Conn_diff)-min(Conn_diff))



#BefAft_tVal[i]=t.test(as.numeric(BefO[i,]),as.numeric(AftO[i,]),paired=TRUE)$statistic;
#AftCtl_tVal[i]=t.test(as.numeric(CtlO[i,]),as.numeric(AftO[i,]))$statistic; 

Tsign=sign((Ctl_Mean -Bef_Mean)*(Aft_Mean - Bef_Mean))

t_diff=Tsign*abs(BefAft_tVal)/(abs(AftCtl_tVal))
unify_t=(t_diff - min(t_diff))/(max(t_diff) - min(t_diff))

Conn_Bef=softConnectivity(t(Bef),corFnc = "bicor",power=6) 
Conn_Aft=softConnectivity(t(Aft),corFnc = "bicor",power=6) 
Conn_Ctl=softConnectivity(t(Ctl),corFnc = "bicor",power=6)

Conn_diff0<-abs(Conn_Bef - Conn_Aft)
Conn_diff1<-abs(Conn_Aft - Conn_Ctl)
Conn
Conn_diff<-Csign*Conn_diff0/(Conn_diff1+0.1)
unify_conn<-(Conn_diff-min(Conn_diff))/(max(Conn_diff)-min(Conn_diff))

weight=(unify_conn + unify_t)/2


####permuation test
# /zs32/home/yxia/R-3.3.2/R-3.3.2/bin/R

setwd("/zs32/home/yxia/04methyZong/20161215datatraits/20170413")
setwd("/zs32/home/yxia/04methyZong/20161215datatraits/20170413")
load(file="/zs32/home/yxia/04methyZong/20161215datatraits/BefO_AftO_CtlO.RData")
library(coin)

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

for (i in 1:nrow(befctl)) 
{
	BefCtl_permutationP[i]=pvalue(oneway_test(befctl[i,]~group,data=data.frame(group,befctl[i,])))
}

result<-data.frame(probeID=rownames(befctl),
	               BefCtl_tVal=BefCtl_tVal,
	               BefCtl_tPval=BefCtl_tPval,
	               BefCtl_permutationP=BefCtl_permutationP,
	               BefMean=BefMean,
	               CtlMean=CtlMean)

save(result,"20170505permuationresults2.RData")