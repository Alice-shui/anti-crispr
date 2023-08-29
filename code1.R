library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)

set.seed(123)
setwd("D:/wechatpublicstation/川大")
getwd()
# 读入表格
## 读入菌株信息
a<-read.table("KP菌株所有详细信息.csv",sep = ",",header = T)
head(a)
## 读入质粒上Acr信息
b<-read.table("比对结果/kpplasmid_136vfacr")
colnames(b)<-c("plasmid","acr","ident","length","mismatch","gapopen","qs","qe","ss","se","ev","bs")
b<-select(b,plasmid,acr,ident,length,ev,bs)
b1<-b%>%filter(ident>40 & ev<1e-5)
head(b1)
## 读入Acr背景信息
c<-read.table("dacr3.csv",sep = ",",header = T)
head(c)



##读入自靶向信息
e<-read.table("selftarge_resultall.csv",sep = ",",header = T)%>%select(name1,V2,type)
colnames(e)<-c("KP","plasmid","type")
e$KP<-gsub(".fna","",e$KP)
e$target="targeted"
head(e)


#将所有信息链接起来
##添加菌株自靶向信息
all1<-left_join(a,e,by="KP")
##添加菌株靶向质粒的acr携带情况
all2<-left_join(all1,b,by="plasmid")
##添加acr的其他信息
all3<-left_join(all2,c,by="acr")
all<-select(all3,KP,CC,ST,plasmid,target,type,ident,ev,bs,acrfamily,cctype,source,Phenotype)
all$count<-1

length(unique(all$KP))

all_t_i40<-filter(all,type=="plasmid"& ident>40 & ev<1e-5)
all_t<-filter(all,type=="plasmid")

length(unique(all_t_i40$KP))
length(unique(all_t$KP))

length(unique(all_t$plasmid))
length(unique(all_t_i40$plasmid))

unique(all_t_i40$ST)
unique(all_t$ST)

length(filter(a,ST==11))
length(unique(filter(all_t_i40,ST==147))$KP)


#绘制菌株CC与自靶向环图
library(circlize)
category<-c("Klebsiella pneumoniae","CRISPR-Cas Carrying strains","Self-target strains")
percent<-c(100,29.89,11.49)
color<-c("#18929E","#55CAA6","#9ECA45")
#初始化
circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 100))
#添加图形
circos.track(
  ylim = c(0.5, length(percent)+0.5), track.height = 0.8,
  bg.border = NA,
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    # 添加圆形中心线
    circos.segments(rep(xlim[1], 3), 1:3,
                    rep(xlim[2], 3), 1:3,
                    col = "#CCCCCC")
    # 添加 3 个圆形矩形
    circos.rect(rep(0, 3),
                1:3 - 0.45,
                percent,
                1:3 + 0.45,
                col = color,
                border = "white")
    # 添加文本信息
    circos.text(
      rep(xlim[1], 3),
      1:3,
      paste(category, " - ", percent, "%"),
      facing = "downward",
      adj = c(1.05, 0.5),
      cex = 0.8
    )
    # 添加轴信息
    breaks = seq(0, 75, by = 5)
    circos.axis(
      h = "top",
      major.at = breaks,
      labels = paste0(breaks, "%"),
      labels.cex = 0.6
    )
  })
circos.clear()

##质粒的复制子类型


##绘制CRISPR-Cas类型柱状图
c1<-c("I-E","IV-A3","I-E&IV-A3","None")
c2<-c(248,29,30,720)
df<-data.frame(c1,c2)
head(df)
colnames(df)<-c("Type","Counts")

library(ggplot2)

df$Type<-factor(df$Type,levels = c("None","I-E&IV-A3","IV-A3","I-E"))

p2<-ggplot(df) +
  aes(x = Type, y = Counts) +
  geom_col(fill = "#55CAA6",width = 0.5) +
  ylim(0,800)+
  coord_flip() +
  theme_test()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1))+
  scale_y_continuous(expand = c(0,0))

p2
##含有CC的菌株
c1<-c("I-E","IV-A3","I-E&IV-A3")
c2<-c(248,29,30)
df<-data.frame(c1,c2)
head(df)
colnames(df)<-c("Type","Counts")

library(ggplot2)

df$Type<-factor(df$Type,levels = c("I-E&IV-A3","IV-A3","I-E"))

p3<-ggplot(df) +
  aes(x = Type, y = Counts) +
  geom_col(fill = "#55CAA6",width = 0.5) +
  coord_flip() +
  theme_test()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1))+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)))

p3


##自靶向菌株中含有CC的类型
c1<-c("I-E","IV-A3","I-E&IV-A3")
c2<-c(63,27,27)
df<-data.frame(c1,c2)
head(df)
colnames(df)<-c("Type","Counts")

library(ggplot2)

df$Type<-factor(df$Type,levels = c("I-E&IV-A3","IV-A3","I-E"))

p4<-ggplot(df) +
  aes(x = Type, y = Counts) +
  geom_col(fill = "#9ECA45",width = 0.5) +
  coord_flip() +
  theme_test()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1))+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)))

p4
#拼接图形
plot_grid(p3,p4, labels=c('b','c'), ncol=1, nrow=2)


#读入质粒类型数据
pt<-read.table("results_tab.csv",sep=",",header = T)
head(pt)
length(unique(pt$Accession.number))
pt<-distinct(pt,Accession.number,Plasmid)
head(pt)

#########################part2,3854个质粒的详细分析，分组比较
pacr<-read.csv("序列/part2/pacr.csv")[,1:3]
head(pacr)
head(c)
pacr1<-left_join(pacr,c,by="acr")%>%distinct(plasmid,acr,.keep_all = TRUE)%>%filter(pident>=39.5)
head(pacr1)

##
read.table("序列/part2/肺克质粒plasmidfinder结果.tsv",sep="\t")





