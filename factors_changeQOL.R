glimpse(df.fl)
fit.dqol <- df.fl %$% 
  lm(d_qolscore ~ as.factor(group) + d_phq + d_pss + d_cesd + d_gsesscore + d_stigmascore + d_insti + d_poscope)
fit.dqol %>% summary()

df.fl %$% t.test(d_qolscore ~ group)

fit1.dqol <- df.fl %$% 
  lm(d_qolscore ~ as.factor(group)+ d_pss + d_cesd + d_gsesscore + d_poscope)
fit1.dqol %>% summary()

anova(fit.dqol, fit1.dqol)

# 标化系数
library(QuantPsyc)
lm.beta(fit1.dqol)
