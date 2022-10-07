################frequency figure###############
library(dplyr)
library(tidyverse)
library(ggplot2)
library(fBasics)
library(car)
library(MASS)
library(ggpubr)
library(ggthemes)
######bac and its 0_2##########################
bac.otu.0_2.freq.tab <- read.csv("bac.otu.0_2.freq.0.01.csv")
head(bac.otu.0_2.freq.tab)
its.otu.0_2.freq.tab <- read.csv("its.otu.0_2.freq.0.01.csv")
head(its.otu.0_2.freq.tab)
bac.its.0_2.freq.tab <- full_join(bac.otu.0_2.freq.tab,its.otu.0_2.freq.tab)#different number of rows combined in the figure.
bac.its.0_2.freq.tab <- bac.its.0_2.freq.tab[2:9]
head(bac.its.0_2.freq.tab)
dim(bac.its.0_2.freq.tab)
colnames(bac.its.0_2.freq.tab)
#bac.its.0_2.freq.tab$freq<- apply(bac.its.0_2.freq.tab,1,function(x)sum(bac.its.0_2.freq.tab$bac.otu.T2W.cor>0&bac.its.0_2.freq.tab$bac.otu.T2W.cor<=0.1))
head(bac.its.0_2.freq.tab)
cor.matrix[1,2] <- sum(bac.its.0_2.freq.tab$bac.otu.T2W.cor>cor.strength[1]&bac.its.0_2.freq.tab$bac.otu.T2W.cor<=cor.strength[2])


cor.strength <- as.numeric(list("-1","-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
cor.strength[2]
cor_freq_count <- function(file){
  cor.matrix <- matrix(nrow=20,ncol = 8)
  colnames(cor.matrix) <- c("bac.T2W.freq","bac.Control.freq","bac.T2.freq","bac.W.freq","its.T2W.freq","its.Control.freq","its.T2.freq","its.W.freq")
  for(b in 1:8){
  for (a in 1:20){
    cor.matrix[a,b]<- sum(file[,b]>cor.strength[a] & file[,b]<=cor.strength[a+1],na.rm=TRUE)
    }
  }
  return(cor.matrix)
  print(b)
}

bac.its.0_2.freq.count <- cor_freq_count(bac.its.0_2.freq.tab)

dim(bac.its.0_2.freq.count)
rownames(bac.its.0_2.freq.count) <- as.factor(c("-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
head(bac.its.0_2.freq.count)
# sumbycol <- apply(bac.its.0_2.freq.count,2,sum)#do a test to make sure count is correct
# sumbycol#correct!



#normalization
bac.its.0_2.freq.count.r <- apply(bac.its.0_2.freq.count,2,function(x)x/sum(x))#calculate relative frequency
head(bac.its.0_2.freq.count.r)

#organize the matrix for plotting
treat <- c("T2W","Control","T2","W")
domain <- c(rep("Bacteria",4),rep("Fungi",4))
range.level <- c(rep((c("-0.9")),8),rep((c("-0.8")),8),rep((c("-0.7")),8),rep((c("-0.6")),8),
                          rep((c("-0.5")),8),rep((c("-0.4")),8),rep((c("-0.3")),8),rep((c("-0.2")),8),
                          rep((c("-0.1")),8),rep((c("0")),8),rep((c("0.1")),8),rep((c("0.2")),8),rep((c("0.3")),8),
                          rep((c("0.4")),8),rep((c("0.5")),8),rep((c("0.6")),8),rep((c("0.7")),8),rep((c("0.8")),8),
                          rep((c("0.9")),8),rep((c("1")),8))
dim(bac.its.0_2.freq.count.r.trans)
bac.its.0_2.freq.count.r.trans <- pivot_longer(as.data.frame(bac.its.0_2.freq.count.r),
                                               cols=c("bac.T2W.freq","bac.Control.freq","bac.T2.freq","bac.W.freq",
                                                      "its.T2W.freq","its.Control.freq","its.T2.freq","its.W.freq")) %>% 
  mutate(Range=range.level, Treatment = c(rep(treat,40)),Domain=c(rep(domain,20)))
head(bac.its.0_2.freq.count.r.trans)
#View(bac.its.0_2.freq.count.r.trans)

bac.its.0_2.freq.count.r.trans$Range <- factor(bac.its.0_2.freq.count.r.trans$Range,levels =c("-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
tiff('bac.its.0_2.freq.count.r.tiff', units="in", width=10, height=5, res=300)
ggplot(bac.its.0_2.freq.count.r.trans, aes(x = Range, y = value, fill = Domain))+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.3)+
  scale_fill_manual(values = c("#f5a21d","aquamarine4"),name="Treatment")+
  labs(title="V0_2",x="Correlation",y="Frequency")+
  facet_wrap(~Treatment)+
  theme_bw()
dev.off()

#strong, moderate and weak category
head(bac.its.0_2.freq.count.r)
head(bac.its.0_2.freq.tab)
hist(bac.its.0_2.freq.tab$bac.otu.T2W.cor)

#overall stronger links
normalTest(sqrt(abs(bac.its.0_2.freq.tab))$bac.otu.T2W.cor,"da")
t.test(abs(bac.its.0_2.freq.tab)$bac.otu.T2W.cor,abs(bac.its.0_2.freq.tab)$its.otu.T2W.cor,paired = F)#fungi
t.test(abs(bac.its.0_2.freq.tab)$bac.otu.T2.cor,abs(bac.its.0_2.freq.tab)$its.otu.T2.cor,paired = F)#bac
t.test(abs(bac.its.0_2.freq.tab)$bac.otu.W.cor,abs(bac.its.0_2.freq.tab)$its.otu.W.cor,paired = F)#fungi
t.test(abs(bac.its.0_2.freq.tab)$bac.otu.Control.cor,abs(bac.its.0_2.freq.tab)$its.otu.Control.cor,paired = F)#bac


######bac and its 2_10##########################
bac.otu.2_10.freq.tab <- read.csv("bac.otu.2_10.freq.0.01.csv")
head(bac.otu.2_10.freq.tab)
its.otu.2_10.freq.tab <- read.csv("its.otu.2_10.freq.0.01.csv")
head(its.otu.2_10.freq.tab)
bac.its.2_10.freq.tab <- full_join(bac.otu.2_10.freq.tab,its.otu.2_10.freq.tab)#different number of rows combined in the figure.
bac.its.2_10.freq.tab <- bac.its.2_10.freq.tab[2:9]
head(bac.its.2_10.freq.tab)
dim(bac.its.2_10.freq.tab)
colnames(bac.its.2_10.freq.tab)
#bac.its.2_10.freq.tab$freq<- apply(bac.its.2_10.freq.tab,1,function(x)sum(bac.its.2_10.freq.tab$bac.otu.T2W.cor>0&bac.its.2_10.freq.tab$bac.otu.T2W.cor<=0.1))
head(bac.its.2_10.freq.tab)
#cor.matrix[1,2] <- sum(bac.its.2_10.freq.tab$bac.otu.T2W.cor>cor.strength[1]&bac.its.2_10.freq.tab$bac.otu.T2W.cor<=cor.strength[2])


# cor.strength <- as.numeric(list("-1","-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
# cor.strength[2]
# cor_freq_count <- function(file){
#   cor.matrix <- matrix(nrow=20,ncol = 8)
#   colnames(cor.matrix) <- c("bac.T2W.freq","bac.Control.freq","bac.T2.freq","bac.W.freq","its.T2W.freq","its.Control.freq","its.T2.freq","its.W.freq")
#   for(b in 1:8){
#     for (a in 1:20){
#       cor.matrix[a,b]<- sum(file[,b]>cor.strength[a] & file[,b]<=cor.strength[a+1],na.rm=TRUE)
#     }
#   }
#   return(cor.matrix)
#   print(b)
# }

bac.its.2_10.freq.count <- cor_freq_count(bac.its.2_10.freq.tab)
dim(bac.its.2_10.freq.count)
rownames(bac.its.2_10.freq.count) <- as.factor(c("-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
head(bac.its.2_10.freq.count)
sumbycol <- apply(bac.its.2_10.freq.count,2,sum)#do a test to make sure count is correct
sumbycol#correct!

#normalization
bac.its.2_10.freq.count.r <- apply(bac.its.2_10.freq.count,2,function(x)x/sum(x))#calculate relative frequency
head(bac.its.2_10.freq.count.r)

#organize the matrix for plotting
treat <- c("T2W","Control","T2","W")
domain <- c(rep("Bacteria",4),rep("Fungi",4))
range.level <- c(rep((c("-0.9")),8),rep((c("-0.8")),8),rep((c("-0.7")),8),rep((c("-0.6")),8),
                 rep((c("-0.5")),8),rep((c("-0.4")),8),rep((c("-0.3")),8),rep((c("-0.2")),8),
                 rep((c("-0.1")),8),rep((c("0")),8),rep((c("0.1")),8),rep((c("0.2")),8),rep((c("0.3")),8),
                 rep((c("0.4")),8),rep((c("0.5")),8),rep((c("0.6")),8),rep((c("0.7")),8),rep((c("0.8")),8),
                 rep((c("0.9")),8),rep((c("1")),8))
dim(bac.its.2_10.freq.count.r.trans)
bac.its.2_10.freq.count.r.trans <- pivot_longer(as.data.frame(bac.its.2_10.freq.count.r),
                                               cols=c("bac.T2W.freq","bac.Control.freq","bac.T2.freq","bac.W.freq",
                                                      "its.T2W.freq","its.Control.freq","its.T2.freq","its.W.freq")) %>% 
  mutate(Range=range.level, Treatment = c(rep(treat,40)),Domain=c(rep(domain,20)))
head(bac.its.2_10.freq.count.r.trans)
#View(bac.its.2_10.freq.count.r.trans)

bac.its.2_10.freq.count.r.trans$Range <- factor(bac.its.2_10.freq.count.r.trans$Range,levels =c("-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
tiff('bac.its.2_10.freq.count.r.tiff', units="in", width=10, height=5, res=300)
ggplot(bac.its.2_10.freq.count.r.trans, aes(x = Range, y = value, fill = Domain))+
  geom_col(alpha = 0.5, position = 'dodge')+
  ylim(0,0.3)+
  scale_fill_manual(values = c("#f5a21d","aquamarine4"),name="Treatment")+
  labs(title="V2_10",x="Correlation",y="Frequency")+
  facet_wrap(~Treatment)+
  theme_bw()
dev.off()

#overall stronger links
t.test(abs(bac.its.2_10.freq.tab)$bac.otu.T2W.cor,abs(bac.its.2_10.freq.tab)$its.otu.T2W.cor,paired = F)#fungi
t.test(abs(bac.its.2_10.freq.tab)$bac.otu.T2.cor,abs(bac.its.2_10.freq.tab)$its.otu.T2.cor,paired = F)#bac
t.test(abs(bac.its.2_10.freq.tab)$bac.otu.W.cor,abs(bac.its.2_10.freq.tab)$its.otu.W.cor,paired = F)#bac
t.test(abs(bac.its.2_10.freq.tab)$bac.otu.Control.cor,abs(bac.its.2_10.freq.tab)$its.otu.Control.cor,paired = F)#bac

######bac and its 20_30##########################
bac.otu.20_30.freq.tab <- read.csv("bac.otu.20_30.freq.0.01.csv")
head(bac.otu.20_30.freq.tab)
its.otu.20_30.freq.tab <- read.csv("its.otu.20_30.freq.0.01.csv")
bac.its.20_30.freq.tab <- full_join(bac.otu.20_30.freq.tab,its.otu.20_30.freq.tab)#different number of rows combined in the figure.
bac.its.20_30.freq.tab <- bac.its.20_30.freq.tab[2:9]
head(bac.its.20_30.freq.tab)
dim(bac.its.20_30.freq.tab)
colnames(bac.its.20_30.freq.tab)
#bac.its.20_30.freq.tab$freq<- apply(bac.its.20_30.freq.tab,1,function(x)sum(bac.its.20_30.freq.tab$bac.otu.T2W.cor>0&bac.its.20_30.freq.tab$bac.otu.T2W.cor<=0.1))
head(bac.its.20_30.freq.tab)
#cor.matrix[1,2] <- sum(bac.its.20_30.freq.tab$bac.otu.T2W.cor>cor.strength[1]&bac.its.20_30.freq.tab$bac.otu.T2W.cor<=cor.strength[2])


# cor.strength <- as.numeric(list("-1","-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
# cor.strength[2]
# cor_freq_count <- function(file){
#   cor.matrix <- matrix(nrow=20,ncol = 8)
#   colnames(cor.matrix) <- c("bac.T2W.freq","bac.Control.freq","bac.T2.freq","bac.W.freq","its.T2W.freq","its.Control.freq","its.T2.freq","its.W.freq")
#   for(b in 1:8){
#     for (a in 1:20){
#       cor.matrix[a,b]<- sum(file[,b]>cor.strength[a] & file[,b]<=cor.strength[a+1],na.rm=TRUE)
#     }
#   }
#   return(cor.matrix)
#   print(b)
# }

bac.its.20_30.freq.count <- cor_freq_count(bac.its.20_30.freq.tab)
dim(bac.its.20_30.freq.count)
rownames(bac.its.20_30.freq.count) <- as.factor(c("-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
head(bac.its.20_30.freq.count)
sumbycol <- apply(bac.its.20_30.freq.count,2,sum)#do a test to make sure count is correct
sumbycol#correct!

#normalization
bac.its.20_30.freq.count.r <- apply(bac.its.20_30.freq.count,2,function(x)x/sum(x))#calculate relative frequency
head(bac.its.20_30.freq.count.r)

#organize the matrix for plotting
treat <- c("T2W","Control","T2","W")
domain <- c(rep("Bacteria",4),rep("Fungi",4))
range.level <- c(rep((c("-0.9")),8),rep((c("-0.8")),8),rep((c("-0.7")),8),rep((c("-0.6")),8),
                 rep((c("-0.5")),8),rep((c("-0.4")),8),rep((c("-0.3")),8),rep((c("-0.2")),8),
                 rep((c("-0.1")),8),rep((c("0")),8),rep((c("0.1")),8),rep((c("0.2")),8),rep((c("0.3")),8),
                 rep((c("0.4")),8),rep((c("0.5")),8),rep((c("0.6")),8),rep((c("0.7")),8),rep((c("0.8")),8),
                 rep((c("0.9")),8),rep((c("1")),8))
dim(bac.its.20_30.freq.count.r.trans)
bac.its.20_30.freq.count.r.trans.raw <- pivot_longer(as.data.frame(bac.its.20_30.freq.count.r),
                                                cols=c("bac.T2W.freq","bac.Control.freq","bac.T2.freq","bac.W.freq",
                                                       "its.T2W.freq","its.Control.freq","its.T2.freq","its.W.freq")) %>% 
  mutate(Range=range.level, Treatment = c(rep(treat,40)),Domain=c(rep(domain,20))) 
  
bac.its.20_30.freq.count.r.trans <- subset(bac.its.20_30.freq.count.r.trans.raw,
                                           !(bac.its.20_30.freq.count.r.trans.raw$Treatment=="T2"&bac.its.20_30.freq.count.r.trans.raw$Domain=="Fungi")) 
head(bac.its.20_30.freq.count.r.trans)
View(bac.its.20_30.freq.count.r.trans)

bac.its.20_30.freq.count.r.trans$Range <- factor(bac.its.20_30.freq.count.r.trans$Range,levels =c("-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
tiff('bac.its.20_30.freq.count.r.tiff', units="in", width=10, height=5, res=300)
ggplot(bac.its.20_30.freq.count.r.trans, aes(x = Range, y = value, fill = Domain))+
  geom_col(alpha = 0.5, position = position_dodge(preserve = "single"))+
  ylim(0,0.3)+
  scale_fill_manual(values = c("#f5a21d","aquamarine4"),name="Treatment")+
  labs(title="V20_30",x="Correlation",y="Frequency")+
  facet_wrap(~Treatment)+
  theme_bw()
dev.off()

#overall stronger links
t.test(abs(bac.its.20_30.freq.tab)$bac.otu.T2W.cor,abs(bac.its.20_30.freq.tab)$its.otu.T2W.cor,paired = F)#bac
#t.test(abs(bac.its.20_30.freq.tab)$bac.otu.T2.cor,abs(bac.its.20_30.freq.tab)$its.otu.T2.cor,paired = F)#bac
t.test(abs(bac.its.20_30.freq.tab)$bac.otu.W.cor,abs(bac.its.20_30.freq.tab)$its.otu.W.cor,paired = F)#bac
t.test(abs(bac.its.20_30.freq.tab)$bac.otu.Control.cor,abs(bac.its.20_30.freq.tab)$its.otu.Control.cor,paired = F)#fungi


#####################statistical test##################
#chisquare http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
#z.test http://www.sthda.com/english/wiki/two-proportions-z-test-in-r
?prop.trend.test()
bac.its.0_2.freq.count
colnames(bac.its.0_2.freq.count)
rownames(bac.its.0_2.freq.count)
######0-2
###control
bac.its.0_2.freq.count.control <- bac.its.0_2.freq.count[,c(2,6)]
head(bac.its.0_2.freq.count.control)
##strong negative
bac.its.0_2.freq.count.control.sum.connection <-apply(bac.its.0_2.freq.count.control,2,sum) 
bac.its.0_2.freq.count.control.sum.connection
bac.its.0_2.freq.count.control.strong.negative <-apply(bac.its.0_2.freq.count.control[c(1:3),],2,sum)#-1 to -0.7
bac.its.0_2.freq.count.control.strong.negative
bac.its.0_2.freq.count.control.strong.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.control.strong.negative,
          n=bac.its.0_2.freq.count.control.sum.connection)
bac.its.0_2.freq.count.control.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.0_2.freq.count.control.sum.connection <-apply(bac.its.0_2.freq.count.control,2,sum) 
bac.its.0_2.freq.count.control.sum.connection
bac.its.0_2.freq.count.control.moderate.negative <-apply(bac.its.0_2.freq.count.control[c(4:7),],2,sum)#-7 to -0.3
bac.its.0_2.freq.count.control.moderate.negative
bac.its.0_2.freq.count.control.moderate.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.control.moderate.negative,
                                                                   n=bac.its.0_2.freq.count.control.sum.connection)
bac.its.0_2.freq.count.control.moderate.negative.z.test

##weak negative
bac.its.0_2.freq.count.control.Weak.negative <-apply(bac.its.0_2.freq.count.control[c(8:10),],2,sum)#
bac.its.0_2.freq.count.control.Weak.negative
bac.its.0_2.freq.count.control.Weak.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.control.Weak.negative,
                                                                   n=c(93553,7490))
bac.its.0_2.freq.count.control.Weak.negative.z.test$

##strong positive
bac.its.0_2.freq.count.control.sum.connection <-apply(bac.its.0_2.freq.count.control,2,sum) 
bac.its.0_2.freq.count.control.sum.connection
bac.its.0_2.freq.count.control.strong.positive <-apply(bac.its.0_2.freq.count.control[c(18:20),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.control.strong.positive
bac.its.0_2.freq.count.control.strong.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.control.strong.positive,
                                                                   n=bac.its.0_2.freq.count.control.sum.connection)
bac.its.0_2.freq.count.control.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.0_2.freq.count.control.sum.connection <-apply(bac.its.0_2.freq.count.control,2,sum) 
bac.its.0_2.freq.count.control.sum.connection
bac.its.0_2.freq.count.control.moderate.positive <-apply(bac.its.0_2.freq.count.control[c(14:17),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.control.moderate.positive
bac.its.0_2.freq.count.control.moderate.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.control.moderate.positive,
                                                                   n=bac.its.0_2.freq.count.control.sum.connection)
bac.its.0_2.freq.count.control.moderate.positive.z.test
##weak positive
bac.its.0_2.freq.count.control.Weak.positive <-apply(bac.its.0_2.freq.count.control[c(11:13),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.control.Weak.positive
bac.its.0_2.freq.count.control.Weak.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.control.Weak.positive,
                                                                 n=bac.its.0_2.freq.count.control.sum.connection)
bac.its.0_2.freq.count.control.Weak.positive.z.test

###T2W

bac.its.0_2.freq.count.T2W <- bac.its.0_2.freq.count[,c(1,5)]
##strong negative
bac.its.0_2.freq.count.T2W.sum.connection <-apply(bac.its.0_2.freq.count.T2W,2,sum) 
bac.its.0_2.freq.count.T2W.sum.connection
bac.its.0_2.freq.count.T2W.strong.negative <-apply(bac.its.0_2.freq.count.T2W[c(1:3),],2,sum)#-1 to -0.7
bac.its.0_2.freq.count.T2W.strong.negative
bac.its.0_2.freq.count.T2W.strong.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.T2W.strong.negative,
                                                                   n=bac.its.0_2.freq.count.T2W.sum.connection)
bac.its.0_2.freq.count.T2W.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.0_2.freq.count.T2W.sum.connection <-apply(bac.its.0_2.freq.count.T2W,2,sum) 
bac.its.0_2.freq.count.T2W.sum.connection
bac.its.0_2.freq.count.T2W.moderate.negative <-apply(bac.its.0_2.freq.count.T2W[c(4:7),],2,sum)#-7 to -0.3
bac.its.0_2.freq.count.T2W.moderate.negative
bac.its.0_2.freq.count.T2W.moderate.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.T2W.moderate.negative,
                                                                     n=bac.its.0_2.freq.count.T2W.sum.connection)
bac.its.0_2.freq.count.T2W.moderate.negative.z.test

##weak negative
bac.its.0_2.freq.count.T2W.Weak.negative <-apply(bac.its.0_2.freq.count.T2W[c(8:10),],2,sum)#
bac.its.0_2.freq.count.T2W.Weak.negative
bac.its.0_2.freq.count.T2W.Weak.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.T2W.Weak.negative,
                                                                 n=c(93553,7490))
bac.its.0_2.freq.count.T2W.Weak.negative.z.test

##strong positive
bac.its.0_2.freq.count.T2W.sum.connection <-apply(bac.its.0_2.freq.count.T2W,2,sum) 
bac.its.0_2.freq.count.T2W.sum.connection
bac.its.0_2.freq.count.T2W.strong.positive <-apply(bac.its.0_2.freq.count.T2W[c(18:20),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.T2W.strong.positive
bac.its.0_2.freq.count.T2W.strong.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.T2W.strong.positive,
                                                                   n=bac.its.0_2.freq.count.T2W.sum.connection)
bac.its.0_2.freq.count.T2W.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.0_2.freq.count.T2W.sum.connection <-apply(bac.its.0_2.freq.count.T2W,2,sum) 
bac.its.0_2.freq.count.T2W.sum.connection
bac.its.0_2.freq.count.T2W.moderate.positive <-apply(bac.its.0_2.freq.count.T2W[c(14:17),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.T2W.moderate.positive
bac.its.0_2.freq.count.T2W.moderate.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.T2W.moderate.positive,
                                                                     n=bac.its.0_2.freq.count.T2W.sum.connection)
bac.its.0_2.freq.count.T2W.moderate.positive.z.test
##weak positive
bac.its.0_2.freq.count.T2W.Weak.positive <-apply(bac.its.0_2.freq.count.T2W[c(11:13),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.T2W.Weak.positive
bac.its.0_2.freq.count.T2W.Weak.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.T2W.Weak.positive,
                                                                 n=bac.its.0_2.freq.count.T2W.sum.connection)
bac.its.0_2.freq.count.T2W.Weak.positive.z.test




###T2
#strong negative
bac.its.0_2.freq.count.T2 <- bac.its.0_2.freq.count[,c(3,7)]
##strong negative
bac.its.0_2.freq.count.T2.sum.connection <-apply(bac.its.0_2.freq.count.T2,2,sum) 
bac.its.0_2.freq.count.T2.sum.connection
bac.its.0_2.freq.count.T2.strong.negative <-apply(bac.its.0_2.freq.count.T2[c(1:3),],2,sum)#-1 to -0.7
bac.its.0_2.freq.count.T2.strong.negative
bac.its.0_2.freq.count.T2.strong.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.T2.strong.negative,
                                                                   n=bac.its.0_2.freq.count.T2.sum.connection)
bac.its.0_2.freq.count.T2.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.0_2.freq.count.T2.sum.connection <-apply(bac.its.0_2.freq.count.T2,2,sum) 
bac.its.0_2.freq.count.T2.sum.connection
bac.its.0_2.freq.count.T2.moderate.negative <-apply(bac.its.0_2.freq.count.T2[c(4:7),],2,sum)#-7 to -0.3
bac.its.0_2.freq.count.T2.moderate.negative
bac.its.0_2.freq.count.T2.moderate.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.T2.moderate.negative,
                                                                     n=bac.its.0_2.freq.count.T2.sum.connection)
bac.its.0_2.freq.count.T2.moderate.negative.z.test

##weak negative
bac.its.0_2.freq.count.T2.Weak.negative <-apply(bac.its.0_2.freq.count.T2[c(8:10),],2,sum)#
bac.its.0_2.freq.count.T2.Weak.negative
bac.its.0_2.freq.count.T2.Weak.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.T2.Weak.negative,
                                                                 n=c(93553,7490))
bac.its.0_2.freq.count.T2.Weak.negative.z.test

##strong positive
bac.its.0_2.freq.count.T2.sum.connection <-apply(bac.its.0_2.freq.count.T2,2,sum) 
bac.its.0_2.freq.count.T2.sum.connection
bac.its.0_2.freq.count.T2.strong.positive <-apply(bac.its.0_2.freq.count.T2[c(18:20),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.T2.strong.positive
bac.its.0_2.freq.count.T2.strong.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.T2.strong.positive,
                                                                   n=bac.its.0_2.freq.count.T2.sum.connection)
bac.its.0_2.freq.count.T2.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.0_2.freq.count.T2.sum.connection <-apply(bac.its.0_2.freq.count.T2,2,sum) 
bac.its.0_2.freq.count.T2.sum.connection
bac.its.0_2.freq.count.T2.moderate.positive <-apply(bac.its.0_2.freq.count.T2[c(14:17),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.T2.moderate.positive
bac.its.0_2.freq.count.T2.moderate.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.T2.moderate.positive,
                                                                     n=bac.its.0_2.freq.count.T2.sum.connection)
bac.its.0_2.freq.count.T2.moderate.positive.z.test$p.value
##weak positive
bac.its.0_2.freq.count.T2.Weak.positive <-apply(bac.its.0_2.freq.count.T2[c(11:13),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.T2.Weak.positive
bac.its.0_2.freq.count.T2.Weak.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.T2.Weak.positive,
                                                                 n=bac.its.0_2.freq.count.T2.sum.connection)
bac.its.0_2.freq.count.T2.Weak.positive.z.test

##W
bac.its.0_2.freq.count.W <- bac.its.0_2.freq.count[,c(4,8)]
##strong negative
bac.its.0_2.freq.count.W.sum.connection <-apply(bac.its.0_2.freq.count.W,2,sum) 
bac.its.0_2.freq.count.W.sum.connection
bac.its.0_2.freq.count.W.strong.negative <-apply(bac.its.0_2.freq.count.W[c(1:3),],2,sum)#-1 to -0.7
bac.its.0_2.freq.count.W.strong.negative
bac.its.0_2.freq.count.W.strong.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.W.strong.negative,
                                                                   n=bac.its.0_2.freq.count.W.sum.connection)
bac.its.0_2.freq.count.W.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.0_2.freq.count.W.sum.connection <-apply(bac.its.0_2.freq.count.W,2,sum) 
bac.its.0_2.freq.count.W.sum.connection
bac.its.0_2.freq.count.W.moderate.negative <-apply(bac.its.0_2.freq.count.W[c(4:7),],2,sum)#-7 to -0.3
bac.its.0_2.freq.count.W.moderate.negative
bac.its.0_2.freq.count.W.moderate.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.W.moderate.negative,
                                                                     n=bac.its.0_2.freq.count.W.sum.connection)
bac.its.0_2.freq.count.W.moderate.negative.z.test

##weak negative
bac.its.0_2.freq.count.W.Weak.negative <-apply(bac.its.0_2.freq.count.W[c(8:10),],2,sum)#
bac.its.0_2.freq.count.W.Weak.negative
bac.its.0_2.freq.count.W.Weak.negative.z.test <- prop.test(x=bac.its.0_2.freq.count.W.Weak.negative,
                                                                 n=c(93553,7490))
bac.its.0_2.freq.count.W.Weak.negative.z.test

##strong positive
bac.its.0_2.freq.count.W.sum.connection <-apply(bac.its.0_2.freq.count.W,2,sum) 
bac.its.0_2.freq.count.W.sum.connection
bac.its.0_2.freq.count.W.strong.positive <-apply(bac.its.0_2.freq.count.W[c(18:20),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.W.strong.positive
bac.its.0_2.freq.count.W.strong.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.W.strong.positive,
                                                                   n=bac.its.0_2.freq.count.W.sum.connection)
bac.its.0_2.freq.count.W.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.0_2.freq.count.W.sum.connection <-apply(bac.its.0_2.freq.count.W,2,sum) 
bac.its.0_2.freq.count.W.sum.connection
bac.its.0_2.freq.count.W.moderate.positive <-apply(bac.its.0_2.freq.count.W[c(14:17),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.W.moderate.positive
bac.its.0_2.freq.count.W.moderate.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.W.moderate.positive,
                                                                     n=bac.its.0_2.freq.count.W.sum.connection)
bac.its.0_2.freq.count.W.moderate.positive.z.test
##weak positive
bac.its.0_2.freq.count.W.Weak.positive <-apply(bac.its.0_2.freq.count.W[c(11:13),],2,sum)#-1 to -0.6
bac.its.0_2.freq.count.W.Weak.positive
bac.its.0_2.freq.count.W.Weak.positive.z.test <- prop.test(x=bac.its.0_2.freq.count.W.Weak.positive,
                                                                 n=bac.its.0_2.freq.count.W.sum.connection)
bac.its.0_2.freq.count.W.Weak.positive.z.test

####################
######2_10
###control
bac.its.2_10.freq.count.control <- bac.its.2_10.freq.count[,c(2,6)]
head(bac.its.2_10.freq.count.control)
##strong negative
bac.its.2_10.freq.count.control.sum.connection <-apply(bac.its.2_10.freq.count.control,2,sum) 
bac.its.2_10.freq.count.control.sum.connection
bac.its.2_10.freq.count.control.strong.negative <-apply(bac.its.2_10.freq.count.control[c(1:3),],2,sum)#-1 to -0.7
bac.its.2_10.freq.count.control.strong.negative
bac.its.2_10.freq.count.control.strong.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.control.strong.negative,
                                                                   n=bac.its.2_10.freq.count.control.sum.connection)
bac.its.2_10.freq.count.control.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.2_10.freq.count.control.sum.connection <-apply(bac.its.2_10.freq.count.control,2,sum) 
bac.its.2_10.freq.count.control.sum.connection
bac.its.2_10.freq.count.control.moderate.negative <-apply(bac.its.2_10.freq.count.control[c(4:7),],2,sum)#-7 to -0.3
bac.its.2_10.freq.count.control.moderate.negative
bac.its.2_10.freq.count.control.moderate.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.control.moderate.negative,
                                                                     n=bac.its.2_10.freq.count.control.sum.connection)
bac.its.2_10.freq.count.control.moderate.negative.z.test

##weak negative
bac.its.2_10.freq.count.control.Weak.negative <-apply(bac.its.2_10.freq.count.control[c(8:10),],2,sum)#
bac.its.2_10.freq.count.control.Weak.negative
bac.its.2_10.freq.count.control.Weak.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.control.Weak.negative,
                                                                 n=c(93553,7490))
bac.its.2_10.freq.count.control.Weak.negative.z.test

##strong positive
bac.its.2_10.freq.count.control.sum.connection <-apply(bac.its.2_10.freq.count.control,2,sum) 
bac.its.2_10.freq.count.control.sum.connection
bac.its.2_10.freq.count.control.strong.positive <-apply(bac.its.2_10.freq.count.control[c(18:20),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.control.strong.positive
bac.its.2_10.freq.count.control.strong.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.control.strong.positive,
                                                                   n=bac.its.2_10.freq.count.control.sum.connection)
bac.its.2_10.freq.count.control.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.2_10.freq.count.control.sum.connection <-apply(bac.its.2_10.freq.count.control,2,sum) 
bac.its.2_10.freq.count.control.sum.connection
bac.its.2_10.freq.count.control.moderate.positive <-apply(bac.its.2_10.freq.count.control[c(14:17),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.control.moderate.positive
bac.its.2_10.freq.count.control.moderate.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.control.moderate.positive,
                                                                     n=bac.its.2_10.freq.count.control.sum.connection)
bac.its.2_10.freq.count.control.moderate.positive.z.test
##weak positive
bac.its.2_10.freq.count.control.Weak.positive <-apply(bac.its.2_10.freq.count.control[c(11:13),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.control.Weak.positive
bac.its.2_10.freq.count.control.Weak.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.control.Weak.positive,
                                                                 n=bac.its.2_10.freq.count.control.sum.connection)
bac.its.2_10.freq.count.control.Weak.positive.z.test

###T2W

bac.its.2_10.freq.count.T2W <- bac.its.2_10.freq.count[,c(1,5)]
##strong negative
bac.its.2_10.freq.count.T2W.sum.connection <-apply(bac.its.2_10.freq.count.T2W,2,sum) 
bac.its.2_10.freq.count.T2W.sum.connection
bac.its.2_10.freq.count.T2W.strong.negative <-apply(bac.its.2_10.freq.count.T2W[c(1:3),],2,sum)#-1 to -0.7
bac.its.2_10.freq.count.T2W.strong.negative
bac.its.2_10.freq.count.T2W.strong.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.T2W.strong.negative,
                                                               n=bac.its.2_10.freq.count.T2W.sum.connection)
bac.its.2_10.freq.count.T2W.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.2_10.freq.count.T2W.sum.connection <-apply(bac.its.2_10.freq.count.T2W,2,sum) 
bac.its.2_10.freq.count.T2W.sum.connection
bac.its.2_10.freq.count.T2W.moderate.negative <-apply(bac.its.2_10.freq.count.T2W[c(4:7),],2,sum)#-7 to -0.3
bac.its.2_10.freq.count.T2W.moderate.negative
bac.its.2_10.freq.count.T2W.moderate.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.T2W.moderate.negative,
                                                                 n=bac.its.2_10.freq.count.T2W.sum.connection)
bac.its.2_10.freq.count.T2W.moderate.negative.z.test

##weak negative
bac.its.2_10.freq.count.T2W.Weak.negative <-apply(bac.its.2_10.freq.count.T2W[c(8:10),],2,sum)#
bac.its.2_10.freq.count.T2W.Weak.negative
bac.its.2_10.freq.count.T2W.Weak.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.T2W.Weak.negative,
                                                             n=c(93553,7490))
bac.its.2_10.freq.count.T2W.Weak.negative.z.test

##strong positive
bac.its.2_10.freq.count.T2W.sum.connection <-apply(bac.its.2_10.freq.count.T2W,2,sum) 
bac.its.2_10.freq.count.T2W.sum.connection
bac.its.2_10.freq.count.T2W.strong.positive <-apply(bac.its.2_10.freq.count.T2W[c(18:20),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.T2W.strong.positive
bac.its.2_10.freq.count.T2W.strong.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.T2W.strong.positive,
                                                               n=bac.its.2_10.freq.count.T2W.sum.connection)
bac.its.2_10.freq.count.T2W.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.2_10.freq.count.T2W.sum.connection <-apply(bac.its.2_10.freq.count.T2W,2,sum) 
bac.its.2_10.freq.count.T2W.sum.connection
bac.its.2_10.freq.count.T2W.moderate.positive <-apply(bac.its.2_10.freq.count.T2W[c(14:17),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.T2W.moderate.positive
bac.its.2_10.freq.count.T2W.moderate.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.T2W.moderate.positive,
                                                                 n=bac.its.2_10.freq.count.T2W.sum.connection)
bac.its.2_10.freq.count.T2W.moderate.positive.z.test
##weak positive
bac.its.2_10.freq.count.T2W.Weak.positive <-apply(bac.its.2_10.freq.count.T2W[c(11:13),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.T2W.Weak.positive
bac.its.2_10.freq.count.T2W.Weak.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.T2W.Weak.positive,
                                                             n=bac.its.2_10.freq.count.T2W.sum.connection)
bac.its.2_10.freq.count.T2W.Weak.positive.z.test




###T2
#strong negative
bac.its.2_10.freq.count.T2 <- bac.its.2_10.freq.count[,c(3,7)]
##strong negative
bac.its.2_10.freq.count.T2.sum.connection <-apply(bac.its.2_10.freq.count.T2,2,sum) 
bac.its.2_10.freq.count.T2.sum.connection
bac.its.2_10.freq.count.T2.strong.negative <-apply(bac.its.2_10.freq.count.T2[c(1:3),],2,sum)#-1 to -0.7
bac.its.2_10.freq.count.T2.strong.negative
bac.its.2_10.freq.count.T2.strong.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.T2.strong.negative,
                                                              n=bac.its.2_10.freq.count.T2.sum.connection)
bac.its.2_10.freq.count.T2.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.2_10.freq.count.T2.sum.connection <-apply(bac.its.2_10.freq.count.T2,2,sum) 
bac.its.2_10.freq.count.T2.sum.connection
bac.its.2_10.freq.count.T2.moderate.negative <-apply(bac.its.2_10.freq.count.T2[c(4:7),],2,sum)#-7 to -0.3
bac.its.2_10.freq.count.T2.moderate.negative
bac.its.2_10.freq.count.T2.moderate.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.T2.moderate.negative,
                                                                n=bac.its.2_10.freq.count.T2.sum.connection)
bac.its.2_10.freq.count.T2.moderate.negative.z.test

##weak negative
bac.its.2_10.freq.count.T2.Weak.negative <-apply(bac.its.2_10.freq.count.T2[c(8:10),],2,sum)#
bac.its.2_10.freq.count.T2.Weak.negative
bac.its.2_10.freq.count.T2.Weak.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.T2.Weak.negative,
                                                            n=c(93553,7490))
bac.its.2_10.freq.count.T2.Weak.negative.z.test

##strong positive
bac.its.2_10.freq.count.T2.sum.connection <-apply(bac.its.2_10.freq.count.T2,2,sum) 
bac.its.2_10.freq.count.T2.sum.connection
bac.its.2_10.freq.count.T2.strong.positive <-apply(bac.its.2_10.freq.count.T2[c(18:20),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.T2.strong.positive
bac.its.2_10.freq.count.T2.strong.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.T2.strong.positive,
                                                              n=bac.its.2_10.freq.count.T2.sum.connection)
bac.its.2_10.freq.count.T2.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.2_10.freq.count.T2.sum.connection <-apply(bac.its.2_10.freq.count.T2,2,sum) 
bac.its.2_10.freq.count.T2.sum.connection
bac.its.2_10.freq.count.T2.moderate.positive <-apply(bac.its.2_10.freq.count.T2[c(14:17),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.T2.moderate.positive
bac.its.2_10.freq.count.T2.moderate.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.T2.moderate.positive,
                                                                n=bac.its.2_10.freq.count.T2.sum.connection)
bac.its.2_10.freq.count.T2.moderate.positive.z.test
##weak positive
bac.its.2_10.freq.count.T2.Weak.positive <-apply(bac.its.2_10.freq.count.T2[c(11:13),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.T2.Weak.positive
bac.its.2_10.freq.count.T2.Weak.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.T2.Weak.positive,
                                                            n=bac.its.2_10.freq.count.T2.sum.connection)
bac.its.2_10.freq.count.T2.Weak.positive.z.test

##W
bac.its.2_10.freq.count.W <- bac.its.2_10.freq.count[,c(4,8)]
##strong negative
bac.its.2_10.freq.count.W.sum.connection <-apply(bac.its.2_10.freq.count.W,2,sum) 
bac.its.2_10.freq.count.W.sum.connection
bac.its.2_10.freq.count.W.strong.negative <-apply(bac.its.2_10.freq.count.W[c(1:3),],2,sum)#-1 to -0.7
bac.its.2_10.freq.count.W.strong.negative
bac.its.2_10.freq.count.W.strong.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.W.strong.negative,
                                                             n=bac.its.2_10.freq.count.W.sum.connection)
bac.its.2_10.freq.count.W.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.2_10.freq.count.W.sum.connection <-apply(bac.its.2_10.freq.count.W,2,sum) 
bac.its.2_10.freq.count.W.sum.connection
bac.its.2_10.freq.count.W.moderate.negative <-apply(bac.its.2_10.freq.count.W[c(4:7),],2,sum)#-7 to -0.3
bac.its.2_10.freq.count.W.moderate.negative
bac.its.2_10.freq.count.W.moderate.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.W.moderate.negative,
                                                               n=bac.its.2_10.freq.count.W.sum.connection)
bac.its.2_10.freq.count.W.moderate.negative.z.test

##weak negative
bac.its.2_10.freq.count.W.Weak.negative <-apply(bac.its.2_10.freq.count.W[c(8:10),],2,sum)#
bac.its.2_10.freq.count.W.Weak.negative
bac.its.2_10.freq.count.W.Weak.negative.z.test <- prop.test(x=bac.its.2_10.freq.count.W.Weak.negative,
                                                           n=c(93553,7490))
bac.its.2_10.freq.count.W.Weak.negative.z.test

##strong positive
bac.its.2_10.freq.count.W.sum.connection <-apply(bac.its.2_10.freq.count.W,2,sum) 
bac.its.2_10.freq.count.W.sum.connection
bac.its.2_10.freq.count.W.strong.positive <-apply(bac.its.2_10.freq.count.W[c(18:20),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.W.strong.positive
bac.its.2_10.freq.count.W.strong.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.W.strong.positive,
                                                             n=bac.its.2_10.freq.count.W.sum.connection)
bac.its.2_10.freq.count.W.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.2_10.freq.count.W.sum.connection <-apply(bac.its.2_10.freq.count.W,2,sum) 
bac.its.2_10.freq.count.W.sum.connection
bac.its.2_10.freq.count.W.moderate.positive <-apply(bac.its.2_10.freq.count.W[c(14:17),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.W.moderate.positive
bac.its.2_10.freq.count.W.moderate.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.W.moderate.positive,
                                                               n=bac.its.2_10.freq.count.W.sum.connection)
bac.its.2_10.freq.count.W.moderate.positive.z.test
##weak positive
bac.its.2_10.freq.count.W.Weak.positive <-apply(bac.its.2_10.freq.count.W[c(11:13),],2,sum)#-1 to -0.6
bac.its.2_10.freq.count.W.Weak.positive
bac.its.2_10.freq.count.W.Weak.positive.z.test <- prop.test(x=bac.its.2_10.freq.count.W.Weak.positive,
                                                           n=bac.its.2_10.freq.count.W.sum.connection)
bac.its.2_10.freq.count.W.Weak.positive.z.test

#################

######20_30
###control
bac.its.20_30.freq.count.control <- bac.its.20_30.freq.count[,c(2,6)]
head(bac.its.20_30.freq.count.control)
##strong negative
bac.its.20_30.freq.count.control.sum.connection <-apply(bac.its.20_30.freq.count.control,2,sum) 
bac.its.20_30.freq.count.control.sum.connection
bac.its.20_30.freq.count.control.strong.negative <-apply(bac.its.20_30.freq.count.control[c(1:3),],2,sum)#-1 to -0.7
bac.its.20_30.freq.count.control.strong.negative
bac.its.20_30.freq.count.control.strong.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.control.strong.negative,
                                                                   n=bac.its.20_30.freq.count.control.sum.connection)
bac.its.20_30.freq.count.control.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.20_30.freq.count.control.sum.connection <-apply(bac.its.20_30.freq.count.control,2,sum) 
bac.its.20_30.freq.count.control.sum.connection
bac.its.20_30.freq.count.control.moderate.negative <-apply(bac.its.20_30.freq.count.control[c(4:7),],2,sum)#-7 to -0.3
bac.its.20_30.freq.count.control.moderate.negative
bac.its.20_30.freq.count.control.moderate.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.control.moderate.negative,
                                                                     n=bac.its.20_30.freq.count.control.sum.connection)
bac.its.20_30.freq.count.control.moderate.negative.z.test

##weak negative
bac.its.20_30.freq.count.control.Weak.negative <-apply(bac.its.20_30.freq.count.control[c(8:10),],2,sum)#
bac.its.20_30.freq.count.control.Weak.negative
bac.its.20_30.freq.count.control.Weak.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.control.Weak.negative,
                                                                 n=c(93553,7490))
bac.its.20_30.freq.count.control.Weak.negative.z.test

##strong positive
bac.its.20_30.freq.count.control.sum.connection <-apply(bac.its.20_30.freq.count.control,2,sum) 
bac.its.20_30.freq.count.control.sum.connection
bac.its.20_30.freq.count.control.strong.positive <-apply(bac.its.20_30.freq.count.control[c(18:20),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.control.strong.positive
bac.its.20_30.freq.count.control.strong.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.control.strong.positive,
                                                                   n=bac.its.20_30.freq.count.control.sum.connection)
bac.its.20_30.freq.count.control.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.20_30.freq.count.control.sum.connection <-apply(bac.its.20_30.freq.count.control,2,sum) 
bac.its.20_30.freq.count.control.sum.connection
bac.its.20_30.freq.count.control.moderate.positive <-apply(bac.its.20_30.freq.count.control[c(14:17),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.control.moderate.positive
bac.its.20_30.freq.count.control.moderate.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.control.moderate.positive,
                                                                     n=bac.its.20_30.freq.count.control.sum.connection)
bac.its.20_30.freq.count.control.moderate.positive.z.test
##weak positive
bac.its.20_30.freq.count.control.Weak.positive <-apply(bac.its.20_30.freq.count.control[c(11:13),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.control.Weak.positive
bac.its.20_30.freq.count.control.Weak.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.control.Weak.positive,
                                                                 n=bac.its.20_30.freq.count.control.sum.connection)
bac.its.20_30.freq.count.control.Weak.positive.z.test

###T2W

bac.its.20_30.freq.count.T2W <- bac.its.20_30.freq.count[,c(1,5)]
##strong negative
bac.its.20_30.freq.count.T2W.sum.connection <-apply(bac.its.20_30.freq.count.T2W,2,sum) 
bac.its.20_30.freq.count.T2W.sum.connection
bac.its.20_30.freq.count.T2W.strong.negative <-apply(bac.its.20_30.freq.count.T2W[c(1:3),],2,sum)#-1 to -0.7
bac.its.20_30.freq.count.T2W.strong.negative
bac.its.20_30.freq.count.T2W.strong.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.T2W.strong.negative,
                                                               n=bac.its.20_30.freq.count.T2W.sum.connection)
bac.its.20_30.freq.count.T2W.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.20_30.freq.count.T2W.sum.connection <-apply(bac.its.20_30.freq.count.T2W,2,sum) 
bac.its.20_30.freq.count.T2W.sum.connection
bac.its.20_30.freq.count.T2W.moderate.negative <-apply(bac.its.20_30.freq.count.T2W[c(4:7),],2,sum)#-7 to -0.3
bac.its.20_30.freq.count.T2W.moderate.negative
bac.its.20_30.freq.count.T2W.moderate.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.T2W.moderate.negative,
                                                                 n=bac.its.20_30.freq.count.T2W.sum.connection)
bac.its.20_30.freq.count.T2W.moderate.negative.z.test

##weak negative
bac.its.20_30.freq.count.T2W.Weak.negative <-apply(bac.its.20_30.freq.count.T2W[c(8:10),],2,sum)#
bac.its.20_30.freq.count.T2W.Weak.negative
bac.its.20_30.freq.count.T2W.Weak.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.T2W.Weak.negative,
                                                             n=c(93553,7490))
bac.its.20_30.freq.count.T2W.Weak.negative.z.test

##strong positive
bac.its.20_30.freq.count.T2W.sum.connection <-apply(bac.its.20_30.freq.count.T2W,2,sum) 
bac.its.20_30.freq.count.T2W.sum.connection
bac.its.20_30.freq.count.T2W.strong.positive <-apply(bac.its.20_30.freq.count.T2W[c(18:20),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.T2W.strong.positive
bac.its.20_30.freq.count.T2W.strong.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.T2W.strong.positive,
                                                               n=bac.its.20_30.freq.count.T2W.sum.connection)
bac.its.20_30.freq.count.T2W.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.20_30.freq.count.T2W.sum.connection <-apply(bac.its.20_30.freq.count.T2W,2,sum) 
bac.its.20_30.freq.count.T2W.sum.connection
bac.its.20_30.freq.count.T2W.moderate.positive <-apply(bac.its.20_30.freq.count.T2W[c(14:17),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.T2W.moderate.positive
bac.its.20_30.freq.count.T2W.moderate.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.T2W.moderate.positive,
                                                                 n=bac.its.20_30.freq.count.T2W.sum.connection)
bac.its.20_30.freq.count.T2W.moderate.positive.z.test
##weak positive
bac.its.20_30.freq.count.T2W.Weak.positive <-apply(bac.its.20_30.freq.count.T2W[c(11:13),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.T2W.Weak.positive
bac.its.20_30.freq.count.T2W.Weak.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.T2W.Weak.positive,
                                                             n=bac.its.20_30.freq.count.T2W.sum.connection)
bac.its.20_30.freq.count.T2W.Weak.positive.z.test




###T2
#strong negative
bac.its.20_30.freq.count.T2 <- bac.its.20_30.freq.count[,c(3,7)]
##strong negative
bac.its.20_30.freq.count.T2.sum.connection <-apply(bac.its.20_30.freq.count.T2,2,sum) 
bac.its.20_30.freq.count.T2.sum.connection
bac.its.20_30.freq.count.T2.strong.negative <-apply(bac.its.20_30.freq.count.T2[c(1:3),],2,sum)#-1 to -0.7
bac.its.20_30.freq.count.T2.strong.negative
bac.its.20_30.freq.count.T2.strong.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.T2.strong.negative,
                                                              n=bac.its.20_30.freq.count.T2.sum.connection)
bac.its.20_30.freq.count.T2.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.20_30.freq.count.T2.sum.connection <-apply(bac.its.20_30.freq.count.T2,2,sum) 
bac.its.20_30.freq.count.T2.sum.connection
bac.its.20_30.freq.count.T2.moderate.negative <-apply(bac.its.20_30.freq.count.T2[c(4:7),],2,sum)#-7 to -0.3
bac.its.20_30.freq.count.T2.moderate.negative
bac.its.20_30.freq.count.T2.moderate.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.T2.moderate.negative,
                                                                n=bac.its.20_30.freq.count.T2.sum.connection)
bac.its.20_30.freq.count.T2.moderate.negative.z.test

##weak negative
bac.its.20_30.freq.count.T2.Weak.negative <-apply(bac.its.20_30.freq.count.T2[c(8:10),],2,sum)#
bac.its.20_30.freq.count.T2.Weak.negative
bac.its.20_30.freq.count.T2.Weak.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.T2.Weak.negative,
                                                            n=c(93553,7490))
bac.its.20_30.freq.count.T2.Weak.negative.z.test

##strong positive
bac.its.20_30.freq.count.T2.sum.connection <-apply(bac.its.20_30.freq.count.T2,2,sum) 
bac.its.20_30.freq.count.T2.sum.connection
bac.its.20_30.freq.count.T2.strong.positive <-apply(bac.its.20_30.freq.count.T2[c(18:20),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.T2.strong.positive
bac.its.20_30.freq.count.T2.strong.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.T2.strong.positive,
                                                              n=bac.its.20_30.freq.count.T2.sum.connection)
bac.its.20_30.freq.count.T2.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.20_30.freq.count.T2.sum.connection <-apply(bac.its.20_30.freq.count.T2,2,sum) 
bac.its.20_30.freq.count.T2.sum.connection
bac.its.20_30.freq.count.T2.moderate.positive <-apply(bac.its.20_30.freq.count.T2[c(14:17),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.T2.moderate.positive
bac.its.20_30.freq.count.T2.moderate.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.T2.moderate.positive,
                                                                n=bac.its.20_30.freq.count.T2.sum.connection)
bac.its.20_30.freq.count.T2.moderate.positive.z.test
##weak positive
bac.its.20_30.freq.count.T2.Weak.positive <-apply(bac.its.20_30.freq.count.T2[c(11:13),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.T2.Weak.positive
bac.its.20_30.freq.count.T2.Weak.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.T2.Weak.positive,
                                                            n=bac.its.20_30.freq.count.T2.sum.connection)
bac.its.20_30.freq.count.T2.Weak.positive.z.test

##W
bac.its.20_30.freq.count.W <- bac.its.20_30.freq.count[,c(4,8)]
##strong negative
bac.its.20_30.freq.count.W.sum.connection <-apply(bac.its.20_30.freq.count.W,2,sum) 
bac.its.20_30.freq.count.W.sum.connection
bac.its.20_30.freq.count.W.strong.negative <-apply(bac.its.20_30.freq.count.W[c(1:3),],2,sum)#-1 to -0.7
bac.its.20_30.freq.count.W.strong.negative
bac.its.20_30.freq.count.W.strong.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.W.strong.negative,
                                                             n=bac.its.20_30.freq.count.W.sum.connection)
bac.its.20_30.freq.count.W.strong.negative.z.test
#p-value = 4.956e-14
# prop 1      prop 2 
# 0.002608147 0.007610147 
##moderate negative
bac.its.20_30.freq.count.W.sum.connection <-apply(bac.its.20_30.freq.count.W,2,sum) 
bac.its.20_30.freq.count.W.sum.connection
bac.its.20_30.freq.count.W.moderate.negative <-apply(bac.its.20_30.freq.count.W[c(4:7),],2,sum)#-7 to -0.3
bac.its.20_30.freq.count.W.moderate.negative
bac.its.20_30.freq.count.W.moderate.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.W.moderate.negative,
                                                               n=bac.its.20_30.freq.count.W.sum.connection)
bac.its.20_30.freq.count.W.moderate.negative.z.test

##weak negative
bac.its.20_30.freq.count.W.Weak.negative <-apply(bac.its.20_30.freq.count.W[c(8:10),],2,sum)#
bac.its.20_30.freq.count.W.Weak.negative
bac.its.20_30.freq.count.W.Weak.negative.z.test <- prop.test(x=bac.its.20_30.freq.count.W.Weak.negative,
                                                           n=c(93553,7490))
bac.its.20_30.freq.count.W.Weak.negative.z.test

##strong positive
bac.its.20_30.freq.count.W.sum.connection <-apply(bac.its.20_30.freq.count.W,2,sum) 
bac.its.20_30.freq.count.W.sum.connection
bac.its.20_30.freq.count.W.strong.positive <-apply(bac.its.20_30.freq.count.W[c(18:20),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.W.strong.positive
bac.its.20_30.freq.count.W.strong.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.W.strong.positive,
                                                             n=bac.its.20_30.freq.count.W.sum.connection)
bac.its.20_30.freq.count.W.strong.positive.z.test
#p-value = 1.242e-07
# sample estimates:
#   prop 1      prop 2 
# 0.004029801 0.008277704
##moderate positive
bac.its.20_30.freq.count.W.sum.connection <-apply(bac.its.20_30.freq.count.W,2,sum) 
bac.its.20_30.freq.count.W.sum.connection
bac.its.20_30.freq.count.W.moderate.positive <-apply(bac.its.20_30.freq.count.W[c(14:17),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.W.moderate.positive
bac.its.20_30.freq.count.W.moderate.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.W.moderate.positive,
                                                               n=bac.its.20_30.freq.count.W.sum.connection)
bac.its.20_30.freq.count.W.moderate.positive.z.test
##weak positive
bac.its.20_30.freq.count.W.Weak.positive <-apply(bac.its.20_30.freq.count.W[c(11:13),],2,sum)#-1 to -0.6
bac.its.20_30.freq.count.W.Weak.positive
bac.its.20_30.freq.count.W.Weak.positive.z.test <- prop.test(x=bac.its.20_30.freq.count.W.Weak.positive,
                                                           n=bac.its.20_30.freq.count.W.sum.connection)
bac.its.20_30.freq.count.W.Weak.positive.z.test



##########################################script for connection statistical test###################################
# head(bac.its.0_2.freq.count.r.trans)
# str(bac.its.0_2.freq.count.r.trans)
# #plot(bac.its.0_2.freq.count.r.trans$value~bac.its.0_2.freq.count.r.trans$Range,data=bac.its.0_2.freq.count.r.trans)
# length(bac.its.0_2.freq.count.r.trans$value)
# length(bac.its.0_2.freq.count.r.trans$Domain)
# sub.positive.strong <- subset(bac.its.0_2.freq.count.r.trans, Range=="-0.8"|Range=="-0.9"|Range=="-0.7"|Range=="-0.6")
# sub.positive.strong.control <- subset(sub.positive.strong, Treatment=="T2W")
# plot(sub.positive.strong.control$value~sub.positive.strong.control$Range,data=sub.positive.strong.control)
# 
# #view(sub.positive.strong)
# mod.test <- lm((value)^(1/2)~Domain,data = sub.positive.strong.control)
# Anova(mod.test)
# test.mod<-glm(value~Domain,family="poisson",data=sub.positive.strong)
# summary(test.mod)
# pchisq(0.018342,24,lower.tail=FALSE)
# Anova(test.mod)
# # normalTest(as.data.frame(bac.its.0_2.freq.count)$its.Control.freq,"da")
# 
# 
# 
# 
# #STEP1
# #The null hypothesis of the Chi-Square test is that no relationship exists on the categorical variables in the population
# #chisqure test
# bac.its.0_2.freq.count
# colnames(bac.its.0_2.freq.count)
# rownames(bac.its.0_2.freq.count)
# ###control
# bac.its.0_2.freq.count.control <- bac.its.0_2.freq.count[,c(2,6)]
# bac.its.0_2.freq.count.T2W <- bac.its.0_2.freq.count[,c(1,5)]
# bac.its.0_2.freq.count.T2 <- bac.its.0_2.freq.count[,c(3,7)]
# bac.its.0_2.freq.count.W <- bac.its.0_2.freq.count[,c(4,8)]
# 
# #control
# #head(bac.its.0_2.freq.count.control)
# # bac.its.0_2.freq.count.control.chi <- chisq.test(t(bac.its.0_2.freq.count.control),simulate.p.value = TRUE)
# # bac.its.0_2.freq.count.control.chi.contri <- 100*bac.its.0_2.freq.count.control.chi$residuals^2/bac.its.0_2.freq.count.control.chi$statistic
# # round(bac.its.0_2.freq.count.control.chi.contri, 3)
# # corrplot(bac.its.0_2.freq.count.control.chi.contri, is.cor = FALSE)
# # corrplot(bac.its.0_2.freq.count.control.chi$residuals, is.cor = FALSE)
# 
# # #T2W
# # head(bac.its.0_2.freq.count.T2W)
# # bac.its.0_2.freq.count.T2W.chi <- chisq.test(t(bac.its.0_2.freq.count.T2W),simulate.p.value = TRUE)
# # bac.its.0_2.freq.count.T2W.chi.contri <- 100*bac.its.0_2.freq.count.T2W.chi$residuals^2/bac.its.0_2.freq.count.T2W.chi$statistic
# # round(bac.its.0_2.freq.count.T2W.chi.contri, 3)
# # corrplot(bac.its.0_2.freq.count.T2W.chi.contri, is.cor = FALSE)
# # #T2
# # head(bac.its.0_2.freq.count.T2)
# # bac.its.0_2.freq.count.T2.chi <- chisq.test(t(bac.its.0_2.freq.count.T2),simulate.p.value = TRUE)
# # bac.its.0_2.freq.count.T2.chi.contri <- 100*bac.its.0_2.freq.count.T2.chi$residuals^2/bac.its.0_2.freq.count.T2.chi$statistic
# # round(bac.its.0_2.freq.count.T2.chi.contri, 3)
# # corrplot(bac.its.0_2.freq.count.T2.chi.contri, is.cor = FALSE)
# # #W
# # head(bac.its.0_2.freq.count.W)
# # bac.its.0_2.freq.count.W.chi <- chisq.test(t(bac.its.0_2.freq.count.W),simulate.p.value = TRUE)
# # bac.its.0_2.freq.count.W.chi.contri <- 100*bac.its.0_2.freq.count.W.chi$residuals^2/bac.its.0_2.freq.count.W.chi$statistic
# # round(bac.its.0_2.freq.count.W.chi.contri, 3)
# # corrplot(bac.its.0_2.freq.count.W.chi.contri, is.cor = FALSE)
# 
# #he row and the column variables are statistically significantly associated (p-value = 0)
# # chi$method
# # chi$expected
# #chisq.test(bac.its.0_2.freq.count.control,simulate.p.value = TRUE)
# ##strong negative
# # bac.its.0_2.freq.count.control.strong.neg <- bac.its.0_2.freq.count.control[c(1,2,3,4),]#-0.9 to -0.6
# # head(bac.its.0_2.freq.count.control.strong.neg)
# # #chisq.test(bac.its.0_2.freq.count.control.strong.neg,simulate.p.value = TRUE)
# # #X-squared = 108.73, df = NA, p-value =0.0004998
# # ##weak negative
# # bac.its.0_2.freq.count.control.weak.neg <- bac.its.0_2.freq.count.control[c(5,6,7,8,9,10),]#-0.5 to 0
# # head(bac.its.0_2.freq.count.control.weak.neg)
# # chisq.test(bac.its.0_2.freq.count.control.weak.neg,simulate.p.value = TRUE)
# # #X-squared = 108.73, df = NA, p-value =0.0004998
# # 
# # ##strong positive
# # bac.its.0_2.freq.count.control.strong.pos <- bac.its.0_2.freq.count.control[c(20,19,18,17),]#-0.9 to -0.6
# # head(bac.its.0_2.freq.count.control.strong.pos)
# # chisq.test(bac.its.0_2.freq.count.control.strong.pos,simulate.p.value = TRUE)
# # 
# # ##weak positive
# # bac.its.0_2.freq.count.control.weak.pos <- bac.its.0_2.freq.count.control[c(16,15,14,13,12,11),]#-0.9 to -0.6
# # head(bac.its.0_2.freq.count.control.weak.pos)
# # chisq.test(bac.its.0_2.freq.count.control.weak.pos,simulate.p.value = TRUE)
# 
# 
# #0_2
# #strong
# head(bac.its.0_2.freq.count.r)
# colnames(bac.its.0_2.freq.count.r)
# rownames(bac.its.0_2.freq.count.r)
# ###control
# 
# bac.its.0_2.negative.strong <- subset(bac.its.0_2.freq.count.r.trans, Range=="-0.8"|Range=="-0.9"|Range=="-0.7"|Range=="-0.6")
# bac.its.0_2.negative.strong.W <- subset(bac.its.0_2.negative.strong,Treatment=="W")
# normalTest(bac.its.0_2.negative.strong.W$value,"sw")
# t.test(bac.its.0_2.negative.strong.W)
# 
# 
# 
# hist(bac.its.0_2.negative.strong.W$value)
# bac.its.0_2.negative.strong.W.mod <- lm(value~Domain,data =bac.its.0_2.negative.strong.W)
# normalTest(bac.its.0_2.negative.strong.W.mod$residuals,"sw")
# hist(bac.its.0_2.negative.strong.W.mod$residuals)
# Anova(bac.its.0_2.negative.strong.W.mod)
# 
# bac.its.0_2.negative.strong.T2W <- subset(bac.its.0_2.negative.strong,Treatment=="T2W")
# bac.its.0_2.negative.strong.T2W.mod <- lm(value~Domain,data =bac.its.0_2.negative.strong.T2W)
# normalTest(bac.its.0_2.negative.strong.T2W.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.strong.T2W.mod)
# plot(bac.its.0_2.negative.strong.T2W.mod)
# bac.its.0_2.negative.strong.Control <- subset(bac.its.0_2.negative.strong,Treatment=="Control")
# bac.its.0_2.negative.strong.Control.mod <- lm(value~Domain,data =bac.its.0_2.negative.strong.Control)
# normalTest(bac.its.0_2.negative.strong.Control.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.strong.Control.mod)
# 
# bac.its.0_2.negative.strong.T2 <- subset(bac.its.0_2.negative.strong,Treatment=="T2") #p=0.02057 * significant bacteria signficant greater than fungi
# bac.its.0_2.negative.strong.T2.mod <- lm(value~Domain,data =bac.its.0_2.negative.strong.T2)
# summary(bac.its.0_2.negative.strong.T2.mod)
# normalTest(bac.its.0_2.negative.strong.T2.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.strong.T2.mod)
# 
# #positive
# bac.its.0_2.positive.strong <- subset(bac.its.0_2.freq.count.r.trans, Range=="1"|Range=="0.9"|Range=="0.8"|Range=="0.7")
# 
# bac.its.0_2.positive.strong.W <- subset(bac.its.0_2.positive.strong,Treatment=="W")
# bac.its.0_2.positive.strong.W.mod <- lm(value~Domain,data =bac.its.0_2.positive.strong.W)
# normalTest(bac.its.0_2.positive.strong.W.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.strong.W.mod)
# 
# bac.its.0_2.positive.strong.T2W <- subset(bac.its.0_2.positive.strong,Treatment=="T2W")
# bac.its.0_2.positive.strong.T2W.mod <- lm(value~Domain,data =bac.its.0_2.positive.strong.T2W)
# normalTest(bac.its.0_2.positive.strong.T2W.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.strong.T2W.mod)
# 
# bac.its.0_2.positive.strong.Control <- subset(bac.its.0_2.positive.strong,Treatment=="Control")
# bac.its.0_2.positive.strong.Control.mod <- lm(value~Domain,data =bac.its.0_2.positive.strong.Control)
# normalTest(bac.its.0_2.positive.strong.Control.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.strong.Control.mod)
# 
# bac.its.0_2.positive.strong.T2 <- subset(bac.its.0_2.positive.strong,Treatment=="T2")
# bac.its.0_2.positive.strong.T2.mod <- lm(value~Domain,data =bac.its.0_2.positive.strong.T2)
# summary(bac.its.0_2.positive.strong.T2.mod)
# normalTest(bac.its.0_2.positive.strong.T2.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.strong.T2.mod)
# 
# #moderate
# bac.its.0_2.negative.moderate <- subset(bac.its.0_2.freq.count.r.trans, Range=="-0.5"|Range=="-0.4"|Range=="-0.3"|Range=="-0.2"|Range=="-0.1"|Range=="0")
# bac.its.0_2.negative.moderate.W <- subset(bac.its.0_2.negative.moderate,Treatment=="W")
# bac.its.0_2.negative.moderate.W.mod <- lm(value~Domain,data =bac.its.0_2.negative.moderate.W)
# normalTest(bac.its.0_2.negative.moderate.W.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.moderate.W.mod)
# 
# bac.its.0_2.negative.moderate.T2W <- subset(bac.its.0_2.negative.moderate,Treatment=="T2W")
# bac.its.0_2.negative.moderate.T2W.mod <- lm(value~Domain,data =bac.its.0_2.negative.moderate.T2W)
# normalTest(bac.its.0_2.negative.moderate.T2W.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.moderate.T2W.mod)
# 
# bac.its.0_2.negative.moderate.Control <- subset(bac.its.0_2.negative.moderate,Treatment=="Control")
# bac.its.0_2.negative.moderate.Control.mod <- lm(value~Domain,data =bac.its.0_2.negative.moderate.Control)
# normalTest(bac.its.0_2.negative.moderate.Control.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.moderate.Control.mod)
# 
# bac.its.0_2.negative.moderate.T2 <- subset(bac.its.0_2.negative.moderate,Treatment=="T2")
# bac.its.0_2.negative.moderate.T2.mod <- lm(value~Domain,data =bac.its.0_2.negative.moderate.T2)
# normalTest(bac.its.0_2.negative.moderate.T2.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.moderate.T2.mod)
# 
# #positive
# bac.its.0_2.positive.moderate <- subset(bac.its.0_2.freq.count.r.trans, Range=="0.6"|Range=="0.5"|Range=="0.4"|Range=="0.3"|Range=="0.2"|Range=="0.1")
# bac.its.0_2.positive.moderate.W <- subset(bac.its.0_2.positive.moderate,Treatment=="W")
# bac.its.0_2.positive.moderate.W.mod <- lm(value~Domain,data =bac.its.0_2.positive.moderate.W)
# normalTest(bac.its.0_2.positive.moderate.W.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.moderate.W.mod)
# 
# bac.its.0_2.positive.moderate.T2W <- subset(bac.its.0_2.positive.moderate,Treatment=="T2W")
# bac.its.0_2.positive.moderate.T2W.mod <- lm(value~Domain,data =bac.its.0_2.positive.moderate.T2W)#p=0.04202 *
# normalTest(bac.its.0_2.positive.moderate.T2W.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.moderate.T2W.mod)
# 
# bac.its.0_2.positive.moderate.Control <- subset(bac.its.0_2.positive.moderate,Treatment=="Control")
# bac.its.0_2.positive.moderate.Control.mod <- lm(value~Domain,data =bac.its.0_2.positive.moderate.Control)
# normalTest(bac.its.0_2.positive.moderate.Control.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.moderate.Control.mod)
# 
# bac.its.0_2.positive.moderate.T2 <- subset(bac.its.0_2.positive.moderate,Treatment=="T2")
# bac.its.0_2.positive.moderate.T2.mod <- lm(value~Domain,data =bac.its.0_2.positive.moderate.T2)
# normalTest(bac.its.0_2.positive.moderate.T2.mod$residuals,"sw")
# anova(bac.its.0_2.positive.moderate.T2.mod)
# 
# ##weak
# #negative
# bac.its.0_2.negative.weak <- subset(bac.its.0_2.freq.count.r.trans, Range=="-0.2"|Range=="-0.1"|Range=="0")
# bac.its.0_2.negative.weak.W <- subset(bac.its.0_2.negative.weak,Treatment=="W")
# bac.its.0_2.negative.weak.W.mod <- lm(value~Domain,data =bac.its.0_2.negative.weak.W)
# normalTest(bac.its.0_2.negative.weak.W.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.weak.W.mod)
# 
# bac.its.0_2.negative.weak.T2W <- subset(bac.its.0_2.negative.weak,Treatment=="T2W")
# bac.its.0_2.negative.weak.T2W.mod <- lm(value~Domain,data =bac.its.0_2.negative.weak.T2W)
# normalTest(bac.its.0_2.negative.weak.T2W.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.weak.T2W.mod)
# 
# bac.its.0_2.negative.weak.Control <- subset(bac.its.0_2.negative.weak,Treatment=="Control")
# bac.its.0_2.negative.weak.Control.mod <- lm(value~Domain,data =bac.its.0_2.negative.weak.Control)
# normalTest(bac.its.0_2.negative.weak.Control.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.weak.Control.mod)
# 
# bac.its.0_2.negative.weak.T2 <- subset(bac.its.0_2.negative.weak,Treatment=="T2")
# bac.its.0_2.negative.weak.T2.mod <- lm(value~Domain,data =bac.its.0_2.negative.weak.T2)
# normalTest(bac.its.0_2.negative.weak.T2.mod$residuals,"sw")
# Anova(bac.its.0_2.negative.weak.T2.mod)
# 
# #positive
# bac.its.0_2.positive.weak <- subset(bac.its.0_2.freq.count.r.trans, Range=="0.3"|Range=="0.2"|Range=="0.1")
# bac.its.0_2.positive.weak.W <- subset(bac.its.0_2.positive.weak,Treatment=="W")
# bac.its.0_2.positive.weak.W.mod <- lm(value~Domain,data =bac.its.0_2.positive.weak.W)
# normalTest(bac.its.0_2.positive.weak.W.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.weak.W.mod)
# 
# bac.its.0_2.positive.weak.T2W <- subset(bac.its.0_2.positive.weak,Treatment=="T2W")
# bac.its.0_2.positive.weak.T2W.mod <- lm(value~Domain,data =bac.its.0_2.positive.weak.T2W)
# normalTest(bac.its.0_2.positive.weak.T2W.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.weak.T2W.mod)
# 
# bac.its.0_2.positive.weak.Control <- subset(bac.its.0_2.positive.weak,Treatment=="Control")
# bac.its.0_2.positive.weak.Control.mod <- lm(value~Domain,data =bac.its.0_2.positive.weak.Control)
# normalTest(bac.its.0_2.positive.weak.Control.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.weak.Control.mod)
# 
# bac.its.0_2.positive.weak.T2 <- subset(bac.its.0_2.positive.weak,Treatment=="T2")
# bac.its.0_2.positive.weak.T2.mod <- lm(value~Domain,data =bac.its.0_2.positive.weak.T2)
# normalTest(bac.its.0_2.positive.weak.T2.mod$residuals,"sw")
# Anova(bac.its.0_2.positive.weak.T2.mod)
# # 
# # 
# # 
# # 
# # head(bac.its.0_2.freq.count.r.trans)
# # bac.its.0_2.freq.count.r.trans.Control <- subset(bac.its.0_2.freq.count.r.trans,Treatment=="Control")
# # head(bac.its.0_2.freq.count.r.trans.Control)
# # bac.its.0_2.freq.count.r.trans.Control.bac <- subset(bac.its.0_2.freq.count.r.trans.Control,Domain=="Bacteria"&Range=="1")
# # bac.its.0_2.freq.count.r.trans.Control.fungi <- subset(bac.its.0_2.freq.count.r.trans.Control,Domain=="Fungi")
# # 
# # normalTest((bac.its.0_2.freq.count.r.trans.Control.bac$value)^(1/2),"da")
# # normalTest((bac.its.0_2.freq.count.r.trans.Control.fungi$value)^(1/2),"da")
# # 
# # bac.its.0_2.freq.count.r.trans.aov.mod <- lm((value)^(1/2)~Domain,data = bac.its.0_2.freq.count.r.trans.Control)
# # summary(bac.its.0_2.freq.count.r.trans.aov.mod)
# # anova(bac.its.0_2.freq.count.r.trans.aov.mod)
# # 
# # 
# # bac.its.0_2.freq.count.r.trans
# # bac.its.2_10.freq.count.r.trans
# # bac.its.20_30.freq.count.r.trans
#0_2
bac.its.0_2.freq.count.control.strong.negative.z.test
bac.its.0_2.freq.count.T2.strong.negative.z.test
bac.its.0_2.freq.count.W.strong.negative.z.test
bac.its.0_2.freq.count.T2W.strong.negative.z.test
bac.its.0_2.freq.count.control.moderate.negative.z.test
bac.its.0_2.freq.count.T2.moderate.negative.z.test
bac.its.0_2.freq.count.W.moderate.negative.z.test
bac.its.0_2.freq.count.T2W.moderate.negative.z.test
bac.its.0_2.freq.count.control.Weak.negative.z.test
bac.its.0_2.freq.count.T2.Weak.negative.z.test
bac.its.0_2.freq.count.W.Weak.negative.z.test
bac.its.0_2.freq.count.T2W.Weak.negative.z.test

bac.its.0_2.freq.count.control.strong.positive.z.test
bac.its.0_2.freq.count.T2.strong.positive.z.test
bac.its.0_2.freq.count.W.strong.positive.z.test
bac.its.0_2.freq.count.T2W.strong.positive.z.test
bac.its.0_2.freq.count.control.moderate.positive.z.test
bac.its.0_2.freq.count.T2.moderate.positive.z.test
bac.its.0_2.freq.count.W.moderate.positive.z.test
bac.its.0_2.freq.count.T2W.moderate.positive.z.test
bac.its.0_2.freq.count.control.Weak.positive.z.test
bac.its.0_2.freq.count.T2.Weak.positive.z.test
bac.its.0_2.freq.count.W.Weak.positive.z.test
bac.its.0_2.freq.count.T2W.Weak.positive.z.test

#2_10
bac.its.2_10.freq.count.control.strong.negative.z.test
bac.its.2_10.freq.count.T2.strong.negative.z.test
bac.its.2_10.freq.count.W.strong.negative.z.test
bac.its.2_10.freq.count.T2W.strong.negative.z.test
bac.its.2_10.freq.count.control.moderate.negative.z.test
bac.its.2_10.freq.count.T2.moderate.negative.z.test
bac.its.2_10.freq.count.W.moderate.negative.z.test
bac.its.2_10.freq.count.T2W.moderate.negative.z.test
bac.its.2_10.freq.count.control.Weak.negative.z.test
bac.its.2_10.freq.count.T2.Weak.negative.z.test
bac.its.2_10.freq.count.W.Weak.negative.z.test
bac.its.2_10.freq.count.T2W.Weak.negative.z.test

bac.its.2_10.freq.count.control.strong.positive.z.test
bac.its.2_10.freq.count.T2.strong.positive.z.test
bac.its.2_10.freq.count.W.strong.positive.z.test
bac.its.2_10.freq.count.T2W.strong.positive.z.test
bac.its.2_10.freq.count.control.moderate.positive.z.test
bac.its.2_10.freq.count.T2.moderate.positive.z.test
bac.its.2_10.freq.count.W.moderate.positive.z.test
bac.its.2_10.freq.count.T2W.moderate.positive.z.test
bac.its.2_10.freq.count.control.Weak.positive.z.test
bac.its.2_10.freq.count.T2.Weak.positive.z.test
bac.its.2_10.freq.count.W.Weak.positive.z.test
bac.its.2_10.freq.count.T2W.Weak.positive.z.test

#20_30
bac.its.20_30.freq.count.control.strong.negative.z.test
bac.its.20_30.freq.count.T2.strong.negative.z.test
bac.its.20_30.freq.count.W.strong.negative.z.test
bac.its.20_30.freq.count.T2W.strong.negative.z.test
bac.its.20_30.freq.count.control.moderate.negative.z.test
bac.its.20_30.freq.count.T2.moderate.negative.z.test
bac.its.20_30.freq.count.W.moderate.negative.z.test
bac.its.20_30.freq.count.T2W.moderate.negative.z.test
bac.its.20_30.freq.count.control.Weak.negative.z.test
bac.its.20_30.freq.count.T2.Weak.negative.z.test
bac.its.20_30.freq.count.W.Weak.negative.z.test
bac.its.20_30.freq.count.T2W.Weak.negative.z.test

bac.its.20_30.freq.count.control.strong.positive.z.test
bac.its.20_30.freq.count.T2.strong.positive.z.test
bac.its.20_30.freq.count.W.strong.positive.z.test
bac.its.20_30.freq.count.T2W.strong.positive.z.test
bac.its.20_30.freq.count.control.moderate.positive.z.test
bac.its.20_30.freq.count.T2.moderate.positive.z.test
bac.its.20_30.freq.count.W.moderate.positive.z.test
bac.its.20_30.freq.count.T2W.moderate.positive.z.test
bac.its.20_30.freq.count.control.Weak.positive.z.test
bac.its.20_30.freq.count.T2.Weak.positive.z.test
bac.its.20_30.freq.count.W.Weak.positive.z.test
bac.its.20_30.freq.count.T2W.Weak.positive.z.test

