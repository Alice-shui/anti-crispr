setwd("E:/Anti-CRISPR-Plasmid/part2")
pacr<-read.csv("E:/Anti-CRISPR-Plasmid/part2/pacr.csv")[,1:4]
dacr<-read.table("E:/Anti-CRISPR-Plasmid/part2/acr_details",fill = F,sep = "|")
colnames(dacr)<-c("acr","protein","acrfamily","cctype","verified","source")


library(dplyr)
library(tidyr)
library(reshape2)
library(ggsci)
library(ggplot2)
library(ggpubr)

setwd("D:/wechatpublicstation/川大/序列/part2/")
pacr<-read.csv("pacr.csv")
dacr<-read.table("acr_details",fill = F,sep = "|")
colnames(dacr)<-c("acr","protein","acrfamily","cctype","verified","source")
#delete space
dacr$acr<-trimws(dacr$acr,which = c("both"),whitespace = "[ ]")
head(dacr)
#combine dataframe
c<-left_join(pacr,dacr,by="acr")
head(c)
#filter
c<-filter(c,c$pident>39,c$length>70)
c$cctype<-trimws(c$cctype,which = c("both"),whitespace = "[ ]")
c<-arrange(c,c$plasmid,c$pident)
head(c)
c2<-c[duplicated(c$plasmid),]
c1<-distinct(c,plasmid,.keep_all = T)
head(c1)

##??????Acr???????ֲ?????
table(c1$acrfamily)

##????anti-IE??Acr??��????
d<-filter(c,c$cctype!="I-E")

#?ϲ???��????
p<-read.table("plasmidfKPplasmid/star_out/plasmidfinder.tsv,header = T)[,2:6]
colnames(p)<-c("ptype","identity","overlap","length_total","plasmid")
p$plasmid<-gsub("ref[|]","",p$plasmid)
p$plasmid<-gsub("[|]","",p$plasmid)
head(p)
c3<-select(c,plasmid,acr,pident,length,protein,acrfamily,cctype)
head(c3)
result<-left_join(c3,p,by="plasmid")
head(result)

head(c3)
c4<-select(c3,plasmid,acrfamily)
m=length(c4[,1])-1
c4$g<-1
for (i in 1:m) {
  if (c4[i+1,1]==c4[i,1]) {
    c4[i+1,3]=c4[i,3]+1
    i=i+1
  }else{
    c4[i+1,3]=1
    i=i+1
  }
}
head(c4)
colnames(c4)<-c("plasmid","value","variab,e")
c4<-dcast(c4,plasmid~variable)
head(c4)
c5<-na.omit(c4)
head(c5)



####????????��????ͳ??
head(result)
test<-resultresultp<-select(result,plasmid,ptype)
p[duplicated(resultp$plasmid),]
#
resultp<-selresultp<-select(result,plasmid,ptype)na(resultp)]<-0
head(resultp)
resultp<-distinct(resultp,plasmid,ptype,.keep_all = T)

m=length(resultp[,1])-1
resultp$g<-1
for (i in 1:m) {
  if (resultp[i+1,1]==resultp[i,1]) {
    resultp[i+1,3]=resultp[i,3]+1
    i=i+1
  }else{
    resultp[i+1,3]=1
    i=i+1
  }
}

colnames(resultp)<-c("plasmid","value","variable")
resultp<-dcast(resultp,plasmid~variable)
#
head(resultp)

#???????ӵ???��??��
k=0
for (i in 1:763) {
  if (resultp[i,2]==0) {
    k=k+1
    i=i+1
    
  }
  
}
k#k=69
  }

}
k#k=69
  }
<-763-k
k1sum(!is.na(resultp$`2`)) #318
sum(!is.na(resultp$`3`)) #18

resultp$ptype<-paste(resultp$`1`,resultp$`2`,
                     resultp$`3`,sep = "_")
resultp$ptype<-gsub("_NA","",resultp$ptype)
head(resultp)
as.data.frame(table(resultp$ptype))

#??????Acr????��??
head(p)
head(result)
nap<-select(result,plasmid,ptype)
nap<-distinct(nap,plasmid,ptype)
ap<-select(p,plasmid,ptype)
ap<-distinct(ap,plasmid,ptype)

head(nap)
head(ap)

na<-anti_join(ap,nap,by="plasmid")
head(na)
na<-arrange(na,plasmid,ptype)

m=length(na[,1])-1
m
na$g<-1
for (i in 1:m) {
  if (na[i+1,1]==na[i,1]) {
    na[i+1,3]=na[i,3]+1
    i=i+1
  }else{
    na[i+1,3]=1
    i=i+1
  }
}
head(na)

colnames(na)<-c("plasmid","value","variable")
na1<-dcast(na,plasmid~variable,value.var = "value")
head(na1)

sum(!is.na(na1$`5`)) #1
sum(!is.na(na1$`4`))-sum(!is.na(na1$`5`)) #2
sum(!is.na(na1$`3`))-2#30
sum(!is.na(na1$`2`))-30#265
length(na1$`1`)-26

a<-as.data.frame(table(na$value))
colnames(a)<-c("ptype","Freq")
a$group<-"na"
b<-as.data.frame(table(nap$ptype))
colnames(b)<-c("ptype","Freq")
b$group<-"p"
b
a
all<-rbind(a,b)
head(all)
new<-data.frame(c("None","None"),
                c(69,1252),
                c("p","na"))
new
names(new)<-c("ptype","Freq","group")
all<-rbind(all,new)

all

all$ptype<-gsub("\\(.{1,20}\\)","",all$ptype)
all$ptype<-gsub("^FIA","IncFIA",all$ptype)
all$ptype<-gsub("^FII","IncFII",all$ptype)
all$ptype<-gsub("^Col.{1,5}","Col",all$ptype)
all$ptype<-gsub("pXuzhou21","other",all$ptype)
all$ptype<-gsub("RepA","other",all$ptype)
all$ptype<-gsub("repA","other",all$ptype)
all<-aggregate(Freq~ptype+group,sum,data = all)
all
all1<-dcast(all,ptype~group,value.var = "Freq")
all1[is.na(all1)]<-0
all1$na<-all1$na/3092*100
all1$p<-all1$p/763*100
all1<-melt(all1,id.vars = "ptype",value.name = "Freq",variable.name = "group")
all1$group<-gsub("na","Acr(-)",all1$group)
all1$group<-gsub("p","Acr(+)",all1$group)
all1
compare_means(Freq~group,data = all1)
mycop<-list(c("Acr(-)","Acr(+)"))

ggplot(all1) +
  aes(x = group, y = Freq) +
  geom_violin(adjust = 3L, boxplot=# "area",trim = T,
              aes(fill = group),show.legend = F) +
  scale_fill_manual(vascale_fill_manual(values = c("#009DAE","#FF7272"))+
  tle("Plasmid differences\Repliconen Acr(+) and Acr(-) groups")+
  theme(plot.title = element_text(hjust = 0.5))
#.y.   group1 group2       p  p.adj p.format p.signif method
#<chr> <chr>  <chr>    <dbl>  <dbl> <chr>    <chr>    <chr>
# Freq  Acr(-) Acr(+) 0.00000000500 0.000000005 5e-09    ****     Wilcoxon


allp1<-all1
allp1$p###########复制子类型显著性差异图bp1t
bp1<-ggplot(all1) +
  aes(x = group,y = Freq,fill = group) +
  geom_boxplot(show.legend = F) +
  theme_test()+
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize=0.1, binwidth = 2,show.legend = F)+
  stat_compare_means(comparisons = mycop,label = "p.signif")+
  scale_fill_manual(values = c("#009DAE","#FF7272"))+
  stat_boxplot(geom = "errorbar",width=0.15)+
  ggtitle("Replicon differences")+
  theme(plot.title = element_text(hjust = 0.5))
bp1







ype[allp1$Freq<5]<-"other"

#绘制质粒分型比较柱图
library(ggplot2)

ap1<-ggplot(allp1) +
     aes(x = ptype, y = Freq, fill = ptype) +
     geom_col() +
     theme_bw() +
     facet_wrap(vars(group),nrow = 2)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  xlab("Plasmid replicon type")+
  ylab("Frequency(%)")+
  scale_y_continuous(expand = c(0,0))

ap1
ap1<-ap1+scale_fill_npg()
ap1
##????????��????ͳ??

head(na1)
head(resultp)
head(result)
pdb<-read.csv("2022pdb.csv",header = T)[,1:2]
colnames(pdb)<-c("plasmid","Lenth")
head(pdb)
test<-sample(pdb$Lenth,5000,replace = T)
shapiro.test(test)

acrp<-resultp[,c(1,5)]
head(test)
head(acrp)
acrp<-left_join(acrp,pdb,by="plasmid")[-2]
head(acrp)
acrp$group<-"Acr(+)"
sum(is.na(acrp$Lenth))

allp<-read.table("../../毕业论文第二部分-质粒上Acr特征研究/本章图表/allplasmid.txt",sep = "\t")
colnames(allp)<-c("plasmid","complete")
table(allp$complete)
nap<-anti_join(allp,acrp,by="plasmid")
head(allp)
nap<-left_join(nap,pdb,by="plasmid")[,-2]
head(nap)
nap$group<-"Acr(-)"
head(acrp)

test<-rbind(acrp,nap)
head(test)
test<-test[sample(nrow(test),5000,replace = T),]
test$n<-c(1:5000)

pp1<-ggplot(test) +
    aes(x = n, y = Lenth,colour = factor(group)) +
    geom_point(shape = "circle", size = 1.8,alpha=0.5) +
    theme_test()+
    scale_colour_manual(values = c("#009DAE","#FF7272"))+
    geom_smooth()+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  xlab("")+
  ylab("Length")

pp1

#
head(acrp)

a<-data.frame()
acrp$Log10<-floor(log10(acrp$Lenth))
count(acrp,Log10)

pdb$Log10<-floor(log10(pdb$Lenth))
head(pdb)
count(pdb,Log10)



for (i in 1:length(acrp$Lenth)) {
  if (log(acrp[i,2])<8&log(acrp[i,2])>7) {
    a[i,1]=acrp[i,1]
    a[i,2]=7
    i=i+1
  }else if (log(acrp[i,2])<9&log(acrp[i,2])>8) {
    a[i,1]=acrp[i,1]
    a[i,2]=8
    i=i+1

  }else if (log(acrp[i,2])<10&log(acrp[i,2])>9) {
    a[i,1]=acrp[i,1]
    a[i,2]=9
    i=i+1

  }else if (log(acrp[i,2])<11&log(acrp[i,2])>10) {
    a[i,1]=acrp[i,1]
    a[i,2]=10
    i=i+1

  }else if (log(acrp[i,2])<12&log(acrp[i,2])>11) {
    a[i,1]=acrp[i,1]
    a[i,2]=11
    i=i+1

  }else if (log(acrp[i,2])<13&log(acrp[i,2])>12) {
    a[i,1]=acrp[i,1]
    a[i,2]=12
    i=i+1

  }else if (log(acrp[i,2])<14&log(acrp[i,2])>13) {
    a[i,1]=acrp[i,1]
    a[i,2]=13
    i=i+1

  }else if (log(acrp[i,2])<15&log(acrp[i,2])>14) {
    a[i,1]=acrp[i,1]
    a[i,2]=14
    i=i+1

  }

}

a1<-a
colnames(a1)<-c("plasmid","Log(length)")
head(a1)
a2<-as.data.frame(table(a1$`Log(length)`))
colnames(a2)<-c("Log(length)","Freq")
a2
a2$group<-c("Acr(+)")

a<-data.frame()
for (i in 1:length(nap$Lenth)) {
  if (log(nap[i,2])<8&log(nap[i,2])>7) {
    a[i,1]=nap[i,1]
    a[i,2]=7
    i=i+1
  }else if (log(nap[i,2])<9&log(nap[i,2])>8) {
    a[i,1]=nap[i,1]
    a[i,2]=8
    i=i+1
  }else if (log(nap[i,2])<10&log(nap[i,2])>9) {
    a[i,1]=nap[i,1]
    a[i,2]=9
    i=i+1
  }else if (log(nap[i,2])<11&log(nap[i,2])>10) {
    a[i,1]=nap[i,1]
    a[i,2]=10
    i=i+1
  }else if (log(nap[i,2])<12&log(nap[i,2])>11) {
    a[i,1]=nap[i,1]
    a[i,2]=11
    i=i+1
  }else if (log(nap[i,2])<13&log(nap[i,2])>12) {
    a[i,1]=nap[i,1]
    a[i,2]=12
    i=i+1
  }else if (log(nap[i,2])<14&log(nap[i,2])>13) {
    a[i,1]=nap[i,1]
    a[i,2]=13
    i=i+1
  }else if (log(nap[i,2])<15&log(nap[i,2])>14) {
    a[i,1]=nap[i,1]
    a[i,2]=14
    i=i+1
  }

}
head(a)
b<-a
colnames(b)<-c("plasmid","Log(length)")
head(b)
b<-as.data.frame(table(b$`Log(length)`))
colnames(b)<-c("Log(length)","Freq")
b
b$group<-c("Acr(-)")
sig<-rbind(a2,b)
sig<-dcast(sig,`Log(length)`~group,value.var = "Freq")
sig[is.na(sig)]<-0
sig<-melt(sig,id.vars = "Log(length)",value.name = "Freq",variable.name = "group")
sig
compare_means(Freq~group,data = sig,method = "t.test")

#> compare_means(Freq~group,data = sig)
# A tibble: 1 x 8
#.y.   group1 group2      p p.adj p.format p.signif method
#<chr> <chr>  <chr>   <dbl> <dbl> <chr>    <chr>    <chr>
#  1 Freq  Acr(-) Acr(+) 0.0284 0.028 0.028    *        Wilcoxon
#> compare_means(Freq~group,data = sig,method = "t.test")
# A tibble: 1 x 8
#.y.   group1 group2      p p.adj p.format p.signif method
#<chr> <chr>  <chr>   <dbl> <dbl> <chr>    <chr>    <chr>
#  1 Freq  Acr(-) Acr(+) 0.0460 0.046 0.046    *        T-test
mycop<-list(c("Acr(+)","Acr(-)"))


ggplot(sig) +
  aes(bp2<-x = group,y = Freq,fill = group) +
  geom_boxplot(show.legend = F) +
  theme_test()+
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize=6, binwidth = 2,show.legend = F)+
  stat_compare_means(comparisons = mycop,label = "p.signif")+
  scale_fill_manual(values = c("#009DAE","#FF7272"))+
  stat_boxplot(geom = "errorbar",width=0.15)+
  ggtitle("Plasmid Length diffe between Acr(+) anle = element_text(hjust = 0.5))



##????????��?bp2ƶ?

??ͳ??
head(p)
acrp<-resultp[,c(1,5)]
head(acrp)
mob<-read.table("Klebsiella_pneumoniae_plasmid_plasmids_classification_sum.txt",
                  fill = T)[,1:2]
head(mob)
colnames(mob)<-c("plasmid","mobility")
mobacr<-left_join(acrp,mob,by="plasmid")
head(mobacr)

length(mobacr$plasmid)
sum(is.na(mobacr))
a<-as.data.frame(table(mobacr$mobility))
a$p<-a$Freq/763*100
a

namob<-anti_join(mob,acrp,by="plasmid")
head(namob)
namob<-filter(namob,namob$mobility=="Conj"|namob$mobility=="mob_unconj"|namob$mobility=="unmob")
namob<-filter(namob,namob$plasmid!="NC_006625.1")
b<-as.data.frame(table(namob$mobility))
b$p<-b$Freq/3091*100
b$group<-c("Acr(-)")
a$group<-c("Acr(+)")

allmob<-rbind(a,b)
allmob
#绘制质粒分型比较柱图
allp2<-allmob
allp2$test="test"
allp2$test[allp2$Var1=="mob_unconj"]<-"mob"
allp2$test[allp2$Var1=="Conj"]<-"conj"
allp2$test[allp2$Var1=="unmob"]<-"unmob"

ap2<-ggplot(allp2) +
  aes(x = test, y = p, fill = test) +
  geom_col(width = 0.33) +
  theme_bw() +
  facet_wrap(vars(group),nrow = 2)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  xlab("Plasmid mobility type")+
  ylab("Frequency(%)")+
  scale_y_continuous(expand = c(0,0))

ap2
ap2<-ap2+scale_fill_npg()
ap2

#拼接图形
plot_grid(ap1,ap2, labels=c('b','c'), ncol=2, nrow=1)


compare_means(Freq~group,data =allmob,method = "t.test")
mycop<-list(c("Acr(+)","Acr(-)"))
t.test(a$p,b$p)

a
b
ks.test(a$p,b$p)
shapiro.test(a$p)
shapiro.test(b$p)

ggplot(allmob) +
  abp3<-es(x = group,y = Freq,fill = group) +
  geom_boxplot(show.legend = F) +
  theme_test()+
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize=6, binwidth = 2,show.legend = F)+
  stat_compare_means(comparisons = mycop,label = "p.signif",method = "t.test")+
  scale_fill_manual(values = c("#009DAE","#FF7272"))+
  stat_boxplot(geom = "errorbar",width=0.15)+
  ggtitle("Plasmid mobility difMn between Acr(+) anle = element_text(hjust = 0.5))

#> compare_means(Frbp3e

#################################拼接差异性图形
plot_grid(bp2,bp1,bp3, labels=c('d','e','f'), ncol=3, nrow=1)

q~group,data =allmob,method = "t.test")
# A tibble: 1 x 8
#.y.   group1 group2      p p.adj p.format p.signif method
#<chr> <chr>  <chr>   <dbl> <dbl> <chr>    <chr>    <chr>
#  1 Freq  Acr(+) Acr(-) 0.0377 0.038 0.038    *        T-test

########??ҩ????
arg<-read.table("Klebsiella_pneumoniae_plasmid_plasmids_classification_sum.txt",
                fill = T)
arg<-filter(arg,arg$V2=="Conj"|arg$V2=="mob_unconj"|arg$V2=="unmob")
arg<-filter(arg,arg$V1!="NC_006625.1")
head(arg)

colnames(arg)<-c("plasmid","mobility","ARG")

#acrp<-resultp[,c(1,5)]
acrp<-resultp
head(acrp)
Acrarg<-left_join(acrp,arg,by="plasmid")
head(Acrarg)
head(acr)
acr<-c1[,c(1,6)]
Acrarg<-left_join(Acrarg,acr,by="plasmid")
head(Acrarg)

write.csv(Acrarg,"acrarg_test.csv")
test<-filter(Acrarg,Acrarg$acrfamily==c("AcrIF11"))

head(argresult)
##????ͼ
Acrarg<-left_join(acrp,betaL,by="plasmid")
head(Acrarg)
Acrarg<-na.omit(Acrarg)
write.csv(Acrarg,"node.csv")
###
arg1<-read.csv("KP_arg.csv")
head(ar../../毕业论文第二部分-质粒上Acr特征研究/本章图表/KP_arg.csvarg1[,1:2]

head(arg2)
arg2<-arrange(arg2,plasmid,ARG)
arg2$ARG<-trimws(arg2$ARG,which = c("both"),whitespace = "[ ]")
unique(arg2$ARG)


NDM<-filter(arg2,ARG=="blaNDM-1"|ARG=="blaNDM-2"|ARG=="blaNDM-3"|ARG=="blaNDM-4"|ARG=="blaNDM-5"|ARG=="blaNDM-6"|ARG=="blaNDM-7"|ARG=="blaNDM-8"|ARG=="blaNDM-9")
NDM$group<-c("NDM")
CMY<-filter(arg2,ARG=="blaCMY-16"|ARG=="blaCMY-2"|ARG=="blaCMY-36"|ARG=="blaCMY-4"|ARG=="blaCMY-54"|ARG=="blaCMY-6"|ARG=="blaCMY-62"|ARG=="blaCMY-8")
CMY$group<-c("CMY")
CTXM<-filter(arg2,AR#G=="blaCTX-M-1"|ARG=="blaCTX-M-104"|ARG=="blaCTX-M-106"|ARG=="blaCTX-M-125"|ARG=="blaCTX-M-14"|ARG=="blaCTX-M-14b"|ARG=="blaCTX-M-15"|ARG=="blaCTX-M-17"|ARG=="blaCTX-M-2"|ARG=="blaCTX-M-24"|ARG=="blaCTX-M-25"|ARG=="blaCTX-M-27"|ARG=="blaCTX-M-3"|ARG=="blaCTX-M-38"||ARG=="blaCTX-M-55"|ARG=="blaCTX-M-59"|ARG=="blaCTX-M-62"|ARG=="blaCTX-M-63"|ARG=="blaCTX-M-64"|ARG=="blaCTX-M-65"|ARG=="blaCTX-M-8"|ARG=="blaCTX-M-90")
CTXM$group<-c("CTXM"CTXM<-arg2%>%filter(grepl('blaCTX-M',ARG))
head(CTXM)
)
DHA<-filter(arg2,ARG=="blaDHA-1"|ARG=="blaDHA-7")
DHA$group<-c("DHA")
IMP<-filter(arg2,ARG=="blaIMP-1"|ARG=="blaIMP-11"|ARG=="blaIMP-26"|ARG=="blaIMP-34"|ARG=="blaIMP-38"|ARG=="blaIMP-4"|ARG=="blaIMP-6"|ARG=="blaIMP-8")
IMP$group<-c("IMP")
OXA<-filter(arg2,ARG=="blaOXA-1"|ARG=="blaOXA-10"|ARG=="blaOXA-181"|ARG=="blaOXA-2"|ARG=="blaOXA-204"|ARG=="blaOXA-21"|ARG=="blaOXA-232"|ARG=="blaOXA-244"|ARG=="blaOXA-427"|ARG=="blaOXA-48"|ARG=="blaOXA-9")
OXA$group<-c("OXA")
KPC<-filter(arg2,ARG=="blaKPC-14"|ARG=="blaKPC-2"|ARG=="blaKPC-3"|ARG=="blaKPC-4"|ARG=="blaKPC-5"|ARG=="blaKPC-6")
KPC$group<-c("KPC")
SHV<-filter(arg2,ARG=="blaSHV-1"|ARG=="blaSHV-11"|ARG=="blaSHV-12"|ARG=="blaSHV-2"|ARG=="blaSHV-28"|ARG=="blaSHV-30"|ARG=="blaSHV-31"|ARG=="blaSHV-5"|ARG=="blaSHV-66")
SHV$group<-c("SHV")
TEM<-filter(arg2,ARG=="blaTEM-15"|ARG=="blaTEM-1A"|ARG=="blaTEM-1B"|ARG=="blaTEM-1C"|ARG=="blaTEM-1D"|ARG=="blaTEM-2"|ARG=="blaTEM-20"|ARG=="blaTEM-206"|ARG=="blaTEM-30"|ARG=="blaTEM-33"|ARG=="blaTEM-94")
TEM$group<-c("TEM")
VIM<-filter(arg2,ARG=="blaVIM-1"|ARG=="blaVIM-19"|ARG=="blaVIM-27"|ARG=="blaVIM-4")
VIM$group<-c("VIM")
VEB<-filter(arg2,ARG=="blaVEB-1"|ARG=="blaVEB-3"|ARG=="blaVEB-8")
VEB$group<-c("VEB")
OTH<-filter(arg2,ARG=="blaBEL-1"|ARG=="blaBKC-1"|ARG=="blaCARB-2"|ARG=="blaFOX-7"|ARG=="blaGES-5"|ARG=="blaPER-1"|ARG=="blaPER-2"|ARG=="blaSCO-1"|ARG=="blaSFO-1"|ARG=="blaSIM-1")
OTH$group<-c("OtherbetaL")
betaL<-rbind(CMY,CTXM,DHA,IMP,OXA,KPC,NDM,SHV,TEM,VIM,VEB,OTH)
betaL$AMR<-"beta-lactamase"
head(betaL)
write.csv(betaL,"betaL.csv")

Amino<-filter(arg2,ARG=="aac(3)-IIa"|ARG=="aac(3)-IId"|ARG=="aac(3)-IVa"|ARG=="aac(3)-Ia"|
                ARG=="aac(3)-Ib"|ARG=="aac(6')-29b"|ARG=="aac(6')-33"|
                ARG=="aac(6')-IIc"| ARG=="aac(6')-Ia"|ARG=="aac(6')-Ian"|
                ARG=="aac(6')-Ib"|ARG=="aac(6')-Ib-11"|ARG=="aac(6')-Ib-Hangzhou"|
                ARG=="aac(6')-Ib-cr"|ARG=="aac(6')-Ib3"|ARG=="aac(6')-Il"|
                ARG=="aac(6')-Im"|ARG=="aac(6')-aph(2'')"|
                ARG=="ant(2'')-Ia"|ARG=="ant(3'')-Ia"|
                ARG=="aph(2'')-Ib"|ARG=="aph(3'')-Ib"|ARG=="aph(3')-IIa"|
                ARG=="aph(3')-Ia"|ARG=="aph(3')-VI"|ARG=="aph(3')-VIa"|
                ARG=="aph(3')-VIb"|ARG=="aph(3')-XV"|ARG=="aph(4)-Ia"|
                ARG=="aph(6)-Id"|ARG=="aph(7'')-Ia"|
                ARG=="aadA1"|ARG=="aadA11"|ARG=="aadA15"|ARG=="aadA16"|
                ARG=="aadA1b"|ARG=="aadA2"|ARG=="aadA22"|
                ARG=="aadA24"|ARG=="aadA5"|ARG=="aadA6"|ARG=="aadA8b"|ARG=="rmtB"|
                ARG== "rmtC" |ARG=="rmtD" |ARG== "rmtD2"|ARG== "rmtG" |ARG== "rmtf")
head(Amino)
Amino$AMR<-"aminoglycoside"

Diami<-filter(arg2, ARG=="dfrA1"|ARG=="dfrA10" |ARG=="dfrA12" |ARG=="dfrA14"
              |ARG=="dfrA15"  |ARG=="dfrA15b"|ARG=="dfrA16"
              |ARG=="dfrA17"  |ARG=="dfrA18" |ARG== "dfrA22"
              |ARG=="dfrA23"  |ARG=="dfrA25" |ARG=="dfrA27"
              |ARG=="dfrA30"  |ARG=="dfrA5" |ARG== "dfrA7"
              |ARG=="dfrA8" |ARG=="dfrB1")
Diami$AMR<-"diaminopyrimidine"

macrolide<-filter(arg2,ARG=="ere(A)"|ARG=="ere(B)" |ARG== "mef(B)"|ARG=="mph(A)"  |ARG=="mph(E)")
macrolide$AMR<-"macrolide"

MLS<-filter(arg2,ARG=="erm(42)" |ARG== "erm(B)"|ARG=="erm(T)")
MLS$AMR<-"MLS"

fosfomycin<-filter(arg2,ARG=="fosA3"|ARG=="fosA5" |ARG== "fosA6" |ARG=="fosE")
fosfomycin$AMR<-"fosfomycin"

lincosamide<-filter(arg2,ARG=="lnu(F)"|ARG=="lnu(G)")
lincosamide$AMR<-"lincosamide"

colistin<-filter(arg2,ARG=="mcr-1"
                 |ARG=="mcr-1.2"  |ARG== "mcr-1.3" |ARG=="mcr-3.1"
                 |ARG=="mcr-3.11" |ARG== "mcr-3.21"|ARG=="mcr-3.22"
                 |ARG=="mcr-3.5" |ARG=="mcr-7.1" |ARG== "mcr-8")
colistin$AMR<-"colistin"


multidrug<-filter(arg2,ARG=="msr(E)" |ARG== "oqxA" |ARG== "oqxB")
multidrug$AMR<-"multidrug"

fluoroquinolone<-filter(arg2,ARG=="qepA1" |ARG=="qepA3" |ARG=="qnrA1"
                        | ARG =="qnrA6" |ARG== "qnrA7" |ARG=="qnrB1"
                        |ARG=="qnrB19" |ARG=="qnrB2" |ARG== "qnrB4"
                        | ARG == "qnrB52"|ARG=="qnrB6"|ARG== "qnrB81"
                        |ARG=="qnrB9"|ARG== "qnrD1" |ARG== "qnrE1"
                        | ARG == "qnrS1" |ARG=="qnrS2"|ARG== "qnrS3"
                        | ARG == "qnrVC4")

fluoroquinolone$AMR<-"fluoroquinolone"

sulfonamide<-filter(arg2,ARG=="sul1"|ARG=="sul2"|ARG=="sul3")
sulfonamide$AMR<-"sulfonamide"

tetracycline<-filter(arg2,ARG=="tet(A)"|ARG=="tet(B)"|ARG=="tet(D)" |ARG=="tet(G)"|ARG=="tet(M)")
tetracycline$AMR<-"tetracycline"

chloramphenicol<-filter(arg2,ARG=="catA1" |ARG=="catB2"|ARG=="catB4")
chloramphenicol$AMR<-"chloramphenicol"

phenicol<-filter(arg2,ARG=="cmlA1" |ARG=="floR")
phenicol$AMR<-"phenicol"

beta_lactamase<-filter(arg2,ARG=="blaBEL-1" |ARG=="blaBKC-1" |ARG=="blaCARB-2"
                       |ARG=="blaCMY-16"  |ARG== "blaCMY-2" |ARG=="blaCMY-36"
                       |ARG=="blaCMY-4" |ARG== "blaCMY-54" |ARG== "blaCMY-6"
                       |ARG=="blaCMY-62" |ARG=="blaCMY-8" |ARG== "blaCTX-M-1"
                       |ARG=="blaCTX-M-104" |ARG=="blaCTX-M-106"|ARG== "blaCTX-M-125"
                       |ARG=="blaCTX-M-14" |ARG== "blaCTX-M-14b" |ARG== "blaCTX-M-15"
                       |ARG=="blaCTX-M-17" |ARG== "blaCTX-M-2" |ARG=="blaCTX-M-24"
                       |ARG=="blaCTX-M-25" |ARG== "blaCTX-M-27" |ARG== "blaCTX-M-3"
                       |ARG=="blaCTX-M-38"  |ARG== "blaCTX-M-55"|ARG== "blaCTX-M-59"
                       |ARG=="blaCTX-M-62" |ARG== "blaCTX-M-63" |ARG== "blaCTX-M-64"
                       |ARG=="blaCTX-M-65" |ARG== "blaCTX-M-8"|ARG== "blaCTX-M-90"
                       |ARG=="blaDHA-1"|ARG== "blaDHA-7"|ARG== "blaFOX-7"
                       |ARG=="blaGES-5" |ARG== "blaIMP-1" |ARG== "blaIMP-11"
                       |ARG=="blaIMP-26"|ARG== "blaIMP-34" |ARG== "blaIMP-38"
                       |ARG=="blaIMP-4" |ARG== "blaIMP-6" |ARG== "blaIMP-8"
                       |ARG=="blaKPC-14" |ARG=="blaKPC-2" |ARG=="blaKPC-3"
                       |ARG=="blaKPC-4" |ARG== "blaKPC-5"|ARG== "blaKPC-6"
                       |ARG=="blaNDM-1" |ARG== "blaNDM-4" |ARG=="blaNDM-5"
                       | ARG == "blaNDM-6" |ARG== "blaNDM-7"|ARG== "blaNDM-9"
                       | ARG == "blaOXA-1" |ARG=="blaOXA-10" |ARG=="blaOXA-181"
                       | ARG == "blaOXA-2" |ARG== "blaOXA-204" |ARG== "blaOXA-21"
                       | ARG == "blaOXA-232"|ARG== "blaOXA-244" |ARG=="blaOXA-427"
                       |ARG=="blaOXA-48"|ARG== "blaOXA-9"|ARG== "blaPER-1"
                       |ARG=="blaPER-2"|ARG== "blaSCO-1" |ARG=="blaSFO-1"
                       |ARG=="blaSHV-1" |ARG== "blaSHV-11" |ARG== "blaSHV-12"
                       |ARG=="blaSHV-2" |ARG== "blaSHV-28" |ARG== "blaSHV-30"
                       |ARG=="blaSHV-31" |ARG== "blaSHV-5" |ARG== "blaSHV-66"
                       |ARG=="blaSIM-1" |ARG== "blaTEM-15" |ARG=="blaTEM-1A"
                       |ARG=="blaTEM-1B" |ARG== "blaTEM-1C" |ARG== "blaTEM-1D"
                       |ARG=="blaTEM-2"|ARG== "blaTEM-20"|ARG== "blaTEM-206"
                       |ARG=="blaTEM-30" |ARG== "blaTEM-33" |ARG=="blaTEM-94"
                       |ARG=="blaVEB-1" |ARG=="blaVEB-3" |ARG== "blaVEB-8"
                       |ARG=="blaVIM-1" |ARG== "blaVIM-19" |ARG== "blaVIM-27"
                       |ARG=="blaVIM-4")

beta_lactamase$AMR<-"beta_lactamase"
head(beta_lactamase)

other<-filter(arg2,ARG=="ARR-2" |ARG=="ARR-3"|ARG=="strA")
other$AMR<-"other"

head(arg1)
head(arg2)
unique(arg2$ARG)
length(unique(allARG$ARG))

allARG<-rbind(beta_lactamase,Amino,Diami,macrolide,MLS,fosfomycin,lincosamide,
              colistin,multidrug,fluoroquinolone,sulfonamide,tetracycline,phenicol,chloramphenicol,other)

as.data.frame(table(allARG$AMR))
head(allARG)
length(unique(allARG$plasmid))
#
acrARG<-left_join(achead(acrp)
rp,allARG,by="plasmid")[,c(1,3,4)]
head(acr[,c(1,6,7)]filter(acrARG,length(unique(acrARG$plasmid))
acrARG$acr<-"Acr(+)"
count(acrARG,AMR)
##绘制Acr有无耐药基因比例图




AMR=="beta_lactamase")
length(unique(test$plasmid))
head(test)

head(betaL)
test<-left_join(acrp,betaL,by="plasmid")
head(test)
test<-na[,c(1,6,7)].omit(test)
length(unique(test$plasmid))
as.data.frame(table(test$group))

head(test)


arrange()
########
##对于Acr(+)
lable<-c("antimicrobial-resistant","sensitive")
data<-c(623,150)hpied1<-data.frame(lable=lable,data=data)
pie1<-ggplot(data = pied1,mapping = aes(x="",y=data,fill=lable))+
  geom_bar(stat = 'identity', position = 'stack')+
  coord_polar(theta = 'y')+
  scale_fill_npg()+
  theme(legend.position = "none")
pie1le2<-c("antimicrobial-resistant","sensitive")
data2<-c(1541,1550)
pied2<-data.frame(lable=lable2,data=data2)

label_value <- paste('(', round(pied2$data/sum(pied2$data) * 100, 1), '%)', sep = '')

pie2<-ggplot(data = pied2,mapping = aes(x="",y=data,fill=lable))+
  geom_bar(stat = 'identity', position = 'stack')+
  coord_polar(theta = 'y')+
  scale_fill_npg()
pie2

#################################拼接图形
plot_grid(pie1,pie2, labels=c('a','b'), ncol=2, nrow=1)

######将耐药结果匹配在3854个质粒上
#所有质粒
head(allp)
#Acr(+)质粒的耐药
head(acrARG)
length(unique(acrARG$plasmid))
#####计算占比
acrARG<-distinct(acrARG,plasmid,AMR,.keep_all = T)
df1<-count(acrARG,acr,AMR)
df1$Frequency<-df1$n/763
head(df1)

#Acr(-)质粒的耐药
##所有耐药质粒
head(allARG)
length(unique(allARG$plasmid))
##获取Acr(-)的耐药
###全是含有耐药的质粒
naacrARG<-anti_join(allARG,acrARG,by="plasmid")
head(naacrARG)
length(unique(naacrARG$plasmid))
###补充不含有耐药的质粒
naacr<-anti_join(allp,acrARG,by="plasmid")
head(naacr)
naacrARG<-left_join(naacr,naacrARG,by="plasmid")
naacrARG$acr<-"Acr(-)"
head(naacrARG)
length(unique(naacrARG$plasmid))
naacrARG<-naacrARG[,-2]
#####计算占比
naacrARG<-distinct(naacrARG,plasmid,AMR,.keep_all = T)
df2<-count(naacrARG,acr,AMR)
df2$Frequency<-df2$n/3091
head(df2)
####合并
allpARG<-rbind(df1,df2)
head(allpARG)

####绘制大类耐药基因比较柱图
library(ggplot2)

p1<-ggplot(allpARG) +
 aes(x = AMR, y = Frequency*100, fill = acr) +
 geom_col(position = "dodge") +
  scale_fill_manual(values = c("#009DAE","#FF7272"))+
 theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  xlab("AMR")+
  ylab("Frequency(%)")+
  scale_y_continuous(expand = c(0,0))



#####################################耐药基因小类β-内酰胺分类比较柱图
head(betaL)
length(unique(betaL$plasmid))
##Acr(+)
head(acrARG)
acrARG1<-acrARG[,c(1,4)]%>%distinct(plasmid,acr)
acrARG1<-left_join(acrARG1,betaL,by='plasmid')[,c(1,2,4)]%>%distinct(plasmid,group,.keep_all = TRUE)
head(acrARG1)
####计算频率
df1<-count(acrARG1,acr,group)
df1$Frequency<-df1$n/763
head(df1)

##Acr(-)
head(naacrARG)
naacrARG1<-naacrARG[,c(1,4)]%>%distinct(plasmid,acr)
naacrARG1<-left_join(naacrARG1,betaL,by='plasmid')[,c(1,2,4)]%>%distinct(plasmid,group,.keep_all = TRUE)
head(naacrARG1)
####计算频率
df2<-count(naacrARG1,acr,group)
df2$Frequency<-df2$n/3091
head(df2)

####合并
allpARG1<-rbind(df1,df2)
head(allpARG1)

####绘制beta-内酰胺耐药基因比较柱图
library(ggplot2)

p2<-ggplot(allpARG1) +
  aes(x = group, y = Frequency*100, fill = acr) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#009DAE","#FF7272"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  xlab("β-lactam")+
  ylab("Frequency(%)")+
  scale_y_continuous(expand = c(0,0))


#################################拼接图形
plot_grid(p1,p2, labels=c('c','d'), ncol=1, nrow=2)





ead(acrp)
head(allARG)
naARG<-anti_join(allARG,acrp,by="plasmid")
head(naARG)
table(naARG$AMR)


#as.data.frame(ͳ??С??
test<-n)alength(unique(naARG$plasmid)).omit(naARG)
unique(test$AMR)

test<-arrange(test,plasmid)[,c(1,2)]
head(test)
length(unique(test$ARG))

test<-na.omit(naARG)
test<-arrange(test,plasmid)[,c(1,3)]
head(test)
length(unique(test$ARG))

test$n<-1
m=length(test$plasmid)-1
for (i in 1:m) {
  if (test[i,1]==test[i+1,1]) {
    test[i+1,3]=test[i,3]+1
    i=i+1

  }else{
    i=i+1
  }

}
test<-dcast(test,plasmid~n,value.var = "ARG")
sum(!is.na(test$`1`))
sum(!is.na(test$`2`))
sum(!is.na(test$`3`))

head(test)

####ͳ?ƴ???
test<-na.omit(naARG)
test<-arrange(test,plasmid)[,c(1,3)]

head(test)

test<-filter(test,AMR=="beta_lactamase")
length(unique(test$plasmid))
test<-filter(test,AMR=="aminoglycoside")
length(unique(test$plasmid))
test<-filter(test,AMR=="sulfonamide")
length(unique(test$plasmid))

####ͳ??Acr(-)BetaL
head(naARG)
test<-left_join(naARG,betaL,by="plasmid")
head(test)
test<-test[,c(1,5)]
test<-distinct(test,plasmid,group)
naacrarg<-arrange(as.data.frame(table(test$group)),Freq)
naacrarg$group<-"Acr(-)"




as.data.frame(table(acrARG$AMR))
test<-filter(acrARG,AMR=="sulfonamide")
length(unique(test$plasmid))
#
acrARG<-left_join(acrp,allARG,by="plasmid")[,c(1,3)]

head(acrARG)
as.data.frame(table(acrARG$ARG))





arg2<-arrange(arg2,plasmid)
head(arg2)
arg2$n<-1
m=length(arg2$plasmid)-1
for (i in 1:m) {
  if (arg2[i,1]==arg2[i+1,1]) {
    arg2[i+1,3]=arg2[i,3]+1
    i=i+1

  }else{
    i=i+1
  }

}

#??????ҩ????????
arg2<-dcast(arg2,plasmid~n,value.var = "ARG")
head(arg2)
sum(!is.na(arg2$`2`))

head(acrARG)
arg32<-left_join(acrp,acrARG,by="plasmid")[,c(1,4)]
head(arg32)
head(betaL)
as.data.frame(table(arg32$ARG))
as.data.frame(table(betaL$ARG))
unique(as.data.frame(table(acrARG$ARG))$Var1)

arg32$n<-1
arg32<-arrange(arg32,plasmid)
head(arg32)

m=length(arg32$plasmid)-1
for (i in 1:m) {
  if (arg32[i,1]==arg32[i+1,1]) {
    arg32[i+1,3]=arg32[i,3]+1
    i=i+1

  }else{
    i=i+1
  }
value.var = "AMR")
head(arg32)

arg3<-left_join(acrp,arg2,by="plasmid")
head(arg3)
arg3<-left_join(arg3,acr,by="plasmid")

test<-filter(arg32,AMR=="sulfonamide")
head(test)
length(unique(test$plasmid))


#write.csv(arg3,"test.csv")
sum(!is.na(arg3$`1`))
sum(!is.na(arg3$`2`))


#????beta_L??ҩ????
head(betaL)
table(betaL$group)
















#bacteria group of IE-CC

kp<-read.csv("E:/Anti-CRISPR-Plasmid/part1/statistics/??ͼ????/KPcc.csv")
kp<-filter(kp,kp$Subtype_probability>0.6)
kp<-select(kp,Contig,Subtype)
kp<-kp%>%distinct(Contig,Subtype)

#kp2<-kp%>%distinct(Subtype,.keep_all = T)
#kp3<-kp%>%distinct(Contig,.keep_all = T)

kp1<-kp%>%dcast(Contig~Subtype)
colnames(kp1)<-c("Contig","c1","c2","c3")
kp1<-unite(kp1,c1,c2,c3,sep = ",")
kp1$c1<-gsub("NA,","",kp1$c1)
kp1$c1<-gsub(",NA","",kp1$c1)
kp1<-filter(kp1,kp$c1!="NA")
kp<-kp1
head(kp)

kparg<-read.csv("E:/Anti-CRISPR-Plasmid/part2/KP????/KP.csv")[,2:8]
colnames(kparg)<-c("Contig","1","ARG","Phenotype","Plasmid","Scheme","ST")
kparg<-kparg[,-2]
head(kparg)
kpresult<-left_join(kp,kparg,by="Contig")
allresult<-full_join(kp,kparg,by="Contig")
head(allresult)
as.data.frame(table(allresult$c1))
write.csv(allresult,"kp??????????ϸ??Ϣ.csv")

#??????NAֵ????
#??????dat[is.na(dat)] <- 0
allresult[is.na(allresult[,2]), 2] <- "None"


##three cc groups of kp
#ST
ieresult<-filter(kpresult,kpresult$c1=="I-E")
ccresult<-filter(kpresult,kpresult$c1!="I-E")
naresult<-anti_join(kparg,kp,by="Contig")
naresult$c1<-c("None")

iest<-ieresult%>%select(Contig,c1,ST)
ccst<-ccresult%>%select(Contig,c1,ST)
nast<-naresult%>%select(Contig,c1,ST)

head(nast)
st<-rbind(iest,ccst,nast)

############################ARGs
##############total ARGs
head(allresult)
argresult<-allresult%>%select(Contig,c1,Phenotype)
head(argresult0)
argresult0<-separate(argresult,Phenotype,
                     into = c("A1","A2","A3","A4","A5","A6","A7","A8",
                              "A9","A10","A11","A12","A13","A14","A15","A16","A17","A18","A19","A20"),
                     sep = ",")

arg0<-filter(argresult0,A1=="Sensitive")[,c(1,2,3)]
argresult0<-filter(argresult0,A1!="Sensitive")

arg<-melt(argresult0,id.vars = c("Contig","c1"))
arg$value<- trimws(arg$value,which = c("both"),
                 whitespace = "[ ]")
argt<-arg[!is.na(arg$value),]
head(argt)
argt[is.na(argt)]<-"None"
argt$c1<-gsub("NA","None",argt$c1)

argt<-select(argt,Contig,c1,value)
head(argt)
colnames(argt)<-c("Contig","CCtype","Phenotype")

arg1<-aggregate(argt,CCtype)

arg1
compare_means(n~CCtype,data = arg1)
mycop<-list(c("I-E","I-E,IV-A3"),
            c("I-E","IV-A3"),
            c("I-E","None"),
            c("I-E,IV-A3","IV-A3"),
            c("I-E,IV-A3","None"),
            c("IV-A3","None"))

ggplot(arg1) +
  aes(x = CCtype, y = log(n), fill = CCtype) +
  geom_boxplot(shape = "circle",show.legend = F) +
  scale_fill_hue(direction = 1) +
  theme_minimal()+
  stat_compare_means(comparisons = mycop)+
  ggtitle("Antibiotic Phenotype differences between different types\n
          of CRISPR-Cas systems")+
  theme(plot.title = element_text(hjust = 0.5))










a<-arg[complete.cases(arg),]
a$value<-gsub("\\[.{3,20}\\]","",a$value)
#ȥ??
a$value<- trimws(a$value,which = c("both"),
                     whitespace = "[ ]")

colnames(a)<-c("Contig","c1","v1","Phenotype")
a<-distinct(a,Contig,c1,Phenotype,.keep_all = T)

head(a)
a$c1<-gsub("NA","None",a$c1)
a<-filter(a,Phenotype!="Sensitive"&Phenotype!="lincomycin")

ggplot(a) +
  aes(x = c1, fill = Phenotype) +
  geom_bar() +
  scale_fill_manual(values = c("#30A9DE","#EFDC05","#E53A40",
                               "#58C9B9","#F17F42","#9055A2",
                               "#E71D36","#4F86C6","#F8F8FF",
                               "#791E94","#f9c00c","#fab1ce",
                               "#f349eb","#00dffc","#008c9e",
                               "#C65146","#EC6A5C","#C5C6B6",
                               "#3ac569","#8c9184","#3a5134",
                               "#9145B6","#ff7f0e"))+
  theme_test()

ggplot(a) +
  aes(x = Phenotype, fill = c1) +
#  geom_bar() +
#  scale_fill_hue(direction = 1) +
  geom_bar(position="fill") +
  scale_fill_manual(values=c("#aec7e8","#ff7f0e","#ffbb78","#1f77b4")) +
   theme_test()+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))



ggplot(a) +
  aes(x = c1, fill = Phenotype) +
  geom_bar() +
  scale_fill_hue(direction = 1) +
  theme_test()
###################heatmap
argresult0$c1<-gsub("NA","None",argresult0$c1)
sum<-count(argresult0,c1)
rownames(sum)<-sum$c1
sum<-sum[,-1]
sum
a1<-count(a,c1,Phenotype)
head(a1)
a1<-dcast(a1,c1~Phenotype,value.var = 'n')
a1[is.na(a1)]<-0
rownames(a1)<-a1$c1
a1<-a1[,-1]
a1
a2<-cbind(a1,sum)
a2

pa<-pheatmap(a2,cluster_rows = F,scale = "row",
             cluster_cols = F,cellheight =12,cellwidth = 24,
             angle_col = 45,fontsize = 10)


#??��??ȡ???ȳ????ݿ?>>??ҩ??????��
for (i in 1:19) {
  assign(paste("arg",i,sep=""),
         filter(argresult0,!is.na(argresult0[,i+2])&is.na(argresult0[,i+3]))[,1:(i+2)])
}

arg20<-filter(argresult0,!is.na(A20))


##########################################################
###############################plasmids
head(allresult)
presult<-allresult%>%select(Contig,c1,Plasmid)
p0<-separate(presult,Plasmid,
                     into = c("A1","A2","A3","A4","A5","A6",
                              "A7","A8","A9"),
                     sep = ",")
head(p0)
p<-melt(p0,id.vars = c("Contig","c1"))
p$value<- trimws(p$value,which = c("both"),
                     whitespace = "[ ]")
p<-distinct(p,Contig,c1,value,.keep_all = T)

pt<-p[!is.na(p$value),]
pt[is.na(pt)]<-"None"
head(pt)
pt<-select(pt,Contig,c1,value)
colnames(pt)<-c("Contig","CCtype","Plasmid")
pt$CCtype<-gsub("NA","None",pt$CCtype)
########
colnames(p)<-c("Contig","c1","v1","Plasmid")
p$Plasmid<-gsub("\\(.{1,20}\\)","",p$Plasmid)
p$Plasmid<-gsub("^FIA","IncFIA",p$Plasmid)
p$Plasmid<-gsub("^FII","IncFII",p$Plasmid)
p$Plasmid<-gsub("^Col.{1,5}","Col",p$Plasmid)
p$Plasmid<-gsub("pXuzhou21","other",p$Plasmid)
p$Plasmid<-gsub("RepA","other",p$Plasmid)
p$Plasmid<-gsub("repA","other",p$Plasmid)

p<-distinct(p,Contig,c1,Plasmid,.keep_all = T)

head(p)
p$c1<-gsub("NA","None",p$c1)
unique(p$c1)

ggplot(p) +
  aes(x = c1, fill = Plasmid) +
  geom_bar() +
  scale_fill_hue(direction = -1) +
  theme_test()

ggplot(p) +
  aes(x = Plasmid, fill = c1) +
 # geom_bar() +
#scale_fill_hue(direction = 1) +
  geom_bar(position="fill") +
  scale_fill_manual(values=c("#aec7e8","#ff7f0e","#ffbb78","#1f77b4")) +
   theme_test()+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))
#################
library(ggpubr)
p1<-count(pt,CCtype,Plasmid)
p1
compare_means(n~CCtype,data = p1)
mycop<-list(c("I-E","I-E,IV-A3"),
            c("I-E","IV-A3"),
            c("I-E","None"),
            c("I-E,IV-A3","IV-A3"),
            c("I-E,IV-A3","None"),
            c("IV-A3","None"))

ggplot(p1) +
  aes(x = CCtype, y = log(n), fill = CCtype) +
  geom_boxplot(shape = "circle",show.legend = F) +
  scale_fill_hue(direction = 1) +
  theme_minimal()+
  stat_compare_means(comparisons = mycop,hide.ns = T,label = "p.signif")+
  ggtitle("Plasmid differences between different types\n
          of CRISPR-Cas systems")+
  theme(plot.title = element_text(hjust = 0.5))

p1<-dcast(p1,Plasmid~c1,value.var = 'n')
p1[is.na(p1)]<-0
rownames(p1)<-p1$Plasmid
p1<-p1[,-1]
p1
p1<-as.matrix(p1)
p1


p2<-cbind(p1,sum)
p2
pp<-pheatmap(p2,scale = "row",cluster_rows = F,
             cluster_cols = F,cellheight =12,
             fontsize = 10,
             cellwidth = 21,angle_col = 45 )
########################
#######################antibiotic gene subtype
head(allresult)
presult<-allresult%>%select(Contig,c1,ARG)
p0<-separate(presult,ARG,
             into = c("A1","A2","A3","A4","A5","A6","A7","A8","A9",
                      "A10","A12","A13","A14","A15",
                      "A16","A17","A18","A19","A20",
                      "A21","A22","A23","A24","A25",
                      "A26","A27","A28","A29","A30",
                      "A31","A32","A33","A34","A35",
                      "A36","A37","A38","A39","A40",
                      "A41","A42","A43","A44","A45",
                      "A46","A47","A48","A49","A50"),
             sep = ",")
p<-melt(p0,id.vars = c("Contig","c1"))
p$value<- trimws(p$value,which = c("both"),whitespace = "[ ]")
p<-distinct(p,Contig,c1,value,.keep_all = T)
p$c1<-gsub("NA","None",p$c1)

colnames(p)<-c("Contig","c1","v1","ARG")

p<-distinct(p,Contig,c1,ARG,.keep_all = T)
p[is.na(p)]<-"None"
head(p)

ggplot(p) +
  aes(x = c1, fill = ARG) +
  geom_bar() +
  scale_fill_hue(direction = -1) +
  theme_test()

ggplot(p) +
  aes(x = ARG, fill = c1) +
  geom_bar(position="fill") +
  scale_fill_manual(values=c("#aec7e8","#ff7f0e","#ffbb78","#1f77b4")) +
  theme_test()+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))

########################

ppt<-p[!is.na(p$value),]
ppt[is.na(ppt)]<-"None"
ppt$c1<-gsub("NA","None",ppt$c1)
head(ppt)
ppt<-select(ppt,Contig,c1,value)
colnames(ppt)<-c("Contig","CCtype","ARG")

pp1<-count(ppt,CCtype,ARG)
pp1
compare_means(n~CCtype,data = pp1)
mycop<-list(c("I-E","I-E,IV-A3"),
            c("I-E","IV-A3"),
            c("I-E","None"),
            c("I-E,IV-A3","IV-A3"),
            c("I-E,IV-A3","None"),
            c("IV-A3","None"))

ggplot(pp1) +
  aes(x = CCtype, y = log(n), fill = CCtype) +
  geom_boxplot(show.legend = F) +
  stat_compare_means(comparisons = mycop)+
  ggtitle("Antibiotic resistence gene differences between different types\n
          of CRISPR-Cas systems")+
  theme(plot.title = element_text(hjust = 0.5))+
  #stat_boxplot(geom = "errorbar",width=0.15)+
  theme_minimal()
################ST
setwd("E:/Anti-CRISPR-Plasmid/part2/ST")
st<-read.csv("KP_ST_CC.csv")
head(st)
sta<-count(st,ST)
stb<-count(st,ST,group)%>%dcast(ST~group)
stb[is.na(stb)]<-0
head(stb)
head(sta)
stb<-left_join(sta,stb,by="ST")%>%arrange(desc(n))
write.csv(stb,"ST??ϸ??.csv")

