#**************************#
# Edited by Qiao Jiaying
# Revised by Xzm
# Last updated: 2018/1/26
# n = 210 （A=105，B=105)
#**************************#

options(digits=3)

# loading packages
library(plyr)
library(dplyr)
library(tidyverse)

setwd("C:/Users/xzm/Desktop/runforlove")

# read data
df <- read_csv("followup/aa-0202.csv")


#定量数据检验
qtest <- function(x, y = df.bl){
  a <- 'normal'
  b <- 'not normal'
  options(digits = 3)
  o1 <- shapiro.test(x)
  o2 <- summary(x)
  o4 <- aggregate(x ~ group, summary, data = y)  
  o3 <- tapply(x, y$group, shapiro.test)
  ifelse(o3$B$p.value < 0.1 | o3$A$p.value < 0.1, o5 <- wilcox.test(x ~ group, y), o5 <- t.test(x ~ group, y))
  #  outcome <- list(b, a)
  return(o5$p.value)
  #  return(outcome)
}

attach(df)
# select variables and range by q1_id
df.sub1 <- df %>% select(q1_id,group,pssscore,poscope,negcope,cesdscore,MET,
                         gsesscore,qolscore,stigmascore,insti,obs) %>% arrange(q1_id)
names(df)
# calculate the change of pss, poscope and negcope, ces-d, MET, gses, qol, stigma
df.sub1.d <- ddply(df.sub1, .(q1_id), mutate, d_pss = pssscore[2] - pssscore[1], 
                   d_poscope = poscope[2] - poscope[1], 
                   d_negcope = negcope[2] - negcope[1], 
                   d_cesd = cesdscore[2] - cesdscore[1], 
                   d_MET = MET[2] - MET[1], 
                   d_gsesscore = gsesscore[2] - gsesscore[1], 
                   d_qolscore = qolscore[2] - qolscore[1], 
                   d_stigmascore = stigmascore[2] - stigmascore[1], 
                   d_insti = insti[2] - insti[1])

# baseline data
df.bl <- filter(df.sub1.d, obs == 1)
# follow data
df.fl <- filter(df.sub1.d, obs == 2)
# baseline of group A 
df.bl.gpA <- df.bl %>% filter(group == 'A')
# baseline of group B
df.bl.gpB <- df.bl %>% filter(group == 'B')
# followup of group A 
df.fl.gpA <- df.fl %>% filter(group == 'A')
# followup of group B
df.fl.gpB <- df.fl %>% filter(group == 'B')


# baseline statistical description 
df.bldes <- df.bl %>% ddply(.(group), summary)
# statistical description after three mongths
df.fldes <- df.fl %>% ddply(.(group), summary)


library(stringr)
#基线##############################################################################################

#提取数值
for (i in 4:ncol(df.bldes)){
  u1<- str_split_fixed(df.bldes[,i], ":",2)
  df.bldes[,i] <- str_trim(u1[,2])
}

#删去多余列
qdes <- df.bldes[,-3]

#定义变量
names(qdes)[2] <- c('var')

#提取分组
u2 <- str_split_fixed(qdes$var, ":",2)
qdes$var <- str_trim(u2[,1])

###分组
qa <- qdes %>% filter(group == 'A')
qb <- qdes %>% filter(group == 'B')

#合并中位数四分位数
Median.A <- str_c(qa[3,], ' (', qa[2,], ', ', qa[5,], ')', collapse = NULL)
Median.B <- str_c(qb[3,], ' (', qb[2,], ', ', qb[5,], ')', collapse = NULL)

#提取标准差，保留两位小数
usa <- round(apply(df.bl.gpA, 2, sd), 2)
usb <- round(apply(df.bl.gpB, 2, sd), 2)

#合并均数标准差
Mean.A <- str_c(qa[4,], ' ± ', usa, collapse = NULL)
Mean.B <- str_c(qb[4,], ' ± ', usb, collapse = NULL)

Variable <- names(qa) %>% str_trim()

#合成表格
bls<- as.data.frame(cbind(Variable, Median.A, Median.B, Mean.A, Mean.B), stringsAsFactors = FALSE)
blsum <- bls[3:11,]

#随访##############################################################################################


#提取数值
for (i in 4:ncol(df.fldes)){
  u1<- str_split_fixed(df.fldes[,i], ":",2)
  df.fldes[,i] <- str_trim(u1[,2])
}

#删去多余列
qdes <- df.fldes[,-3]

#定义变量
names(qdes)[2] <- c('var')

#提取分组
u2 <- str_split_fixed(qdes$var, ":", 2)
qdes$var <- str_trim(u2[, 1])

###分组
qa <- qdes %>% filter(group == 'A')
qb <- qdes %>% filter(group == 'B')

#合并中位数四分位数
Median.A<- str_c(qa[3,],' (',qa[2,], ', ', qa[5,],')', collapse = NULL)
Median.B<- str_c(qb[3,],' (',qb[2,], ', ', qb[5,],')', collapse = NULL)

#提取标准差，保留两位小数
usa <- round(apply(df.fl.gpA, 2, sd), 2)
usb <- round(apply(df.fl.gpB, 2, sd), 2)

#合并均数标准差
Mean.A <- str_c(qa[4,],' ± ', usa, collapse = NULL)
Mean.B <- str_c(qb[4,],' ± ', usb, collapse = NULL)

Variable <- names(qa) %>% str_trim

#合成表格
fls <- as.data.frame(cbind(Variable, Median.A, Median.B, Mean.A, Mean.B), stringsAsFactors = FALSE)

#生成随访描述
flsum <- fls[3:11, ]

#生成差值描述
dsum <-  fls[13:21, ]

###############################################Baseline test#####
#Variable to be tested
a <- select(df.bl, pssscore, poscope, negcope, cesdscore, MET, gsesscore, qolscore, stigmascore, insti)
#outcome output
u <- list()
for(i in 1:ncol(a)){
  #print(names(a[i]))
  u1 <- a[, i] %>% qtest()
  #u <- list(u,u1)
  #u[i] <- u1
  u <- c(u, u1)
  names(u)[i] <- names(a[i])
}

#blsum
p.a <- reshape2::melt(u)
names(p.a) <- c('P vales', 'variable')
blsum.p <- merge(blsum, p.a, by.x = 'Variable', by.y = 'variable')

############################################Follow-up test#####
#Variable to be tested
a <- select(df.fl, pssscore, poscope, negcope, cesdscore, MET, gsesscore, qolscore, stigmascore, insti)
#outcome output
u <- list()
for(i in 1:ncol(a)){
  #  print(names(a[i]))
  u1 <- a[,i] %>% qtest(y=df.fl)
  #  u <- list(u,u1)
  u <- c(u, u1)
  names(u)[i] <- names(a[i])
}
#flsum
p.a <- reshape2::melt(u)
names(p.a) <- c('P vales', 'variable')
flsum.p <- merge(flsum, p.a, by.x = 'Variable', by.y = 'variable')

###############################################d test#####
#Variable to be tested
a <- df.bl[,c('d_pss','d_poscope','d_negcope','d_cesd','d_MET',
              'd_gsesscore','d_qolscore','d_stigmascore','d_insti')]

#outcome output
u <- list()
for(i in 1:ncol(a)){
  #  print(names(a[i]))
  u1 <- a[,i] %>% qtest()
  #  u <- list(u,u1)
  u <- c(u, u1)
  names(u)[i] <- names(a[i])
}

#dsum
p.a <- reshape2::melt(u)
names(p.a) <- c('P vales', 'variable')
dsum.p <- merge(dsum, p.a, by.x = 'Variable', by.y = 'variable')


################################三线表导出到WORD##################################################
#加载包
library(ReporteRs)
library(magrittr)

#设定三线表vanilla.table
options( "ReporteRs-fontsize" = 11 )

ft_obj <- vanilla.table(flsum.p)
bt_obj <- vanilla.table(blsum.p)
d_obj <- vanilla.table(dsum.p)

#导出到word
doc = docx( )
doc = addFlexTable(doc, flextable = bt_obj )
doc = addFlexTable(doc, flextable = ft_obj )
doc = addFlexTable(doc, flextable = d_obj )
writeDoc(doc, file = "add_ft_ex.docx")
