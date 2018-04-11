library(plyr)
library(dplyr)
library(tidyr)

a<- read.table("./file.txt",stringsAsFactors = F)
data=data.frame(name=character(), dnds=numeric(), stringsAsFactors = F)

i<-a$V1[1]


file= data.frame(V1=character(), stringsAsFactors = F)
for (i in a$V1){
  fas<-paste0( as.character(i))
  fas1<-paste0( as.character(i))
  output<-readLines(paste0(fas1))
  n=length(output)
  x=substr(output[n], 1,9)
  if (x=="Time used") {
    file[nrow(file)+1,]<-c(i)
  }
  
}

data<-data.frame(gene.ID=character(),taxan=numeric(),t=numeric(),N=numeric(),S=numeric(),dnds=numeric(),dN=numeric(),ds=numeric(),NdN=numeric(),  SdS=numeric(),lnl=numeric(), stringsAsFactors = F)
for (i in file$V1){
  fas<-paste0( as.character(i))
  fas1<-paste0(as.character(i))
  output<-readLines(paste0(fas1))
  n<-length(output)
  for (j in 1:n) { #looking for the line that has the log likelihood value
    subset<-substr(output[j],1,7)
    if (subset==" branch") {
      l<-j
    }
  }
  l<-l+2
  for (j in 1:n) { #looking for the line that has the log likelihood value
    subset<-substr(output[j],1,3)
    if (subset=="lnL") {
      lnl<-j
    }
  }
  
  ltxt<-output[lnl]
  like<-as.numeric(unlist(regmatches(ltxt,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",ltxt))))
  lnl<- -1*like[3]
  
  for (k in l:(l+5)){
    txt<-output[k]
    txt<-substr(txt,7,nchar(txt))
    num<-as.numeric(unlist(regmatches(txt,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",txt))))
    
    data[nrow(data)+1,]<-c(substr(i,1,nchar(i)-4),num,lnl)
    
  }
  
}

taxa<-read.csv("../taxan.csv", stringsAsFactors = F)
d<-data

d<- merge(d, taxa, by="taxan")

zGGA<- d[d$taxa=="GGA",]
zTGU<- d[d$taxa=="TGU",]
zASI<- d[d$taxa=="ASI",]
zGGH<- d[d$taxa=="GGH",]


write.csv(d, "adata.csv")
autoset<-read.csv("adata.csv", stringsAsFactors = F)



xset<-read.csv("../Zdata.csv", stringsAsFactors = F)

autoset<-autoset[autoset$ds<2,]

asum<-autoset %>% group_by(taxan)%>% summarise(N=sum(N),S=sum(S),NdN=sum(NdN),SdS=sum(SdS))
asum$dN<-asum$NdN/asum$N
asum$dS<-asum$SdS/asum$S
asum$dNds<-asum$dN/asum$dS
adnds<-asum$dNds
adn<-asum$dN
ads<-asum$dS
asum<- merge(asum, taxa, by="taxan")

xset<-xset[xset$ds<2,]
xsum<-xset %>% group_by(taxan)%>% summarise(N=sum(N),S=sum(S),NdN=sum(NdN),SdS=sum(SdS))
xsum$dN<-xsum$NdN/xsum$N
xsum$dS<-xsum$SdS/xsum$S
xsum$dNds<-xsum$dN/xsum$dS
xsum$dataset<-c("Z")
asum$dataset<-c("A")

xsum<- merge(xsum, taxa, by="taxan")
ss<-rbind(asum,xsum)

# replace TGU with species of interest

adnds<-asum$dNds[asum$taxa=="ASI"]
adn<-asum$dN[asum$taxa=="ASI"]
ads<-asum$dS[asum$taxa=="ASI"]
zdnds<-xsum$dNds[xsum$taxa=="ASI"]
zdn<-xsum$dN[xsum$taxa=="ASI"]
zds<-xsum$dS[xsum$taxa=="ASI"]
zdnds
adnds
aa<-autoset[autoset$taxa=="ASI",]
xx<-xset[xset$taxa=="ASI",]
data<-rbind(aa,xx)
data$taxa<-NULL

data$gene.ID<-NULL

p1<-0
p2<-0
for (i in 1:10000){
  n<-sample(1:nrow(data),nrow(aa)) 
  ap<-as.matrix.data.frame(data[c(n),])
  zp<-as.matrix.data.frame(data[-c(n),])
  #ap2<-as.matrix.data.frame(para[c(n2),])
  #zp2<-as.matrix.data.frame(para[-c(n2),])
  
  aasum<-colSums(ap)
  zzsum<-colSums(zp)
  
  
  padnds<-as.numeric((aasum[9]/aasum[4])/(aasum[10]/aasum[5]))
  pzdnds<-as.numeric((zzsum[9]/zzsum[4])/(zzsum[10]/zzsum[5]))
  
  #padnds<-asum[4]
  #pzdnds<-zsum[4]
  
  if ((pzdnds-padnds)>abs(zdnds-adnds)){
    p1=p1+1
  }
  if ((pzdnds-padnds)<(-abs(zdnds-adnds))){
    p2=p2+1
  }
}
p1
p2
pvaluednds<-(p1+p2)/10000
#
#
#
#
#
#
#dn
p1<-0
p2<-0
for (i in 1:10000){
  n<-sample(1:nrow(data),nrow(aa)) 
  ap<-as.matrix.data.frame(data[c(n),])
  zp<-as.matrix.data.frame(data[-c(n),])
  #ap2<-as.matrix.data.frame(para[c(n2),])
  #zp2<-as.matrix.data.frame(para[-c(n2),])
  
  aasum<-colSums(ap)
  zzsum<-colSums(zp)
  
  
  padn<-as.numeric((aasum[9]/aasum[4]))
  pzdn<-as.numeric((zzsum[9]/zzsum[4]))
  
  #padnds<-asum[4]
  #pzdnds<-zsum[4]
  
  if ((pzdn-padn)>abs(zdn-adn)){
    p1=p1+1
  }
  if ((pzdn-padn)<(-abs(zdn-adn))){
    p2=p2+1
  }
}
p1
p2
pvaluedn<-(p1+p2)/10000


#ds
p1<-0
p2<-0
for (i in 1:10000){
  n<-sample(1:nrow(data),nrow(aa)) 
  ap<-as.matrix.data.frame(data[c(n),])
  zp<-as.matrix.data.frame(data[-c(n),])
  #ap2<-as.matrix.data.frame(para[c(n2),])
  #zp2<-as.matrix.data.frame(para[-c(n2),])
  
  aasum<-colSums(ap)
  zzsum<-colSums(zp)
  
  
  pads<-as.numeric((aasum[10]/aasum[5]))
  pzds<-as.numeric((zzsum[10]/zzsum[5]))
  
  #padnds<-asum[4]
  #pzdnds<-zsum[4]
  
  if ((pzds-pads)>abs(zds-ads)){
    p1=p1+1
  }
  if ((pzds-pads)<(-abs(zds-ads))){
    p2=p2+1
  }
}
p1
p2
pvalueds<-(p1+p2)/10000



pvaluednds
pvaluedn
pvalueds


#bootstrap
aboot<-data.frame(dnds=as.numeric(),dn=as.numeric(), ds=as.numeric(), stringsAsFactors = F)
zboot<-data.frame(dnds=as.numeric(),dn=as.numeric(), ds=as.numeric(),  stringsAsFactors = F)


aa<-autoset[autoset$taxa=="ASI",]
zz<-xset[xset$taxa=="ASI",]

aa$taxa<-NULL
aa$gene.ID<-NULL

zz$taxa<-NULL
zz$gene.ID<-NULL

for (i in 1:10000){
  n<-sample(1:nrow(aa),replace = T)
  ap<-as.matrix.data.frame(aa[c(n),])
  n<-sample(1:nrow(zz),replace = T)
  zp<-as.matrix.data.frame(zz[c(n),])
  aasum<-colSums(ap)
  zzsum<-colSums(zp)
  
  badnds<-c(dnds=as.numeric((aasum[9]/aasum[4])/(aasum[10]/aasum[5])),dn=as.numeric((aasum[9]/aasum[4])),ds=as.numeric((aasum[10]/aasum[5])))
  bzdnds<-c(dnds=as.numeric((zzsum[9]/zzsum[4])/(zzsum[10]/zzsum[5])),dn=as.numeric((zzsum[9]/zzsum[4])),ds=as.numeric((zzsum[10]/zzsum[5])))
  aboot[nrow(aboot)+1,]<-c(badnds)
  zboot[nrow(zboot)+1,]<-c(bzdnds)
  
}
a<-aboot$dnds[order(aboot$dnds)]
a[250]
a[10000-250]

a<-aboot$dn[order(aboot$dn)]
a[250]
a[10000-250]

a<-aboot$ds[order(aboot$ds)]
a[250]
a[10000-250]


z<-zboot$dnds[order(zboot$dnds)]
z[250]
z[10000-250]

z<-zboot$dn[order(zboot$dn)]
z[250]
z[10000-250]

z<-zboot$ds[order(zboot$ds)]
z[250]
z[10000-250]
