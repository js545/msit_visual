# MSIT MIND Manuscript
# Load packages
library(ggpubr)

###############################################################
# Response Time Outliers by Condition
rt_df = read.csv('E:/Data/MSIT_MIND/group_RT_129_modified_v2.csv')

rt_cont = rt_df[which(rt_df$Condition == 'Control'),]
rt_mean = mean(rt_cont$Response_Time)
rt_stdv = sd(rt_cont$Response_Time)
rt_cont = subset(rt_cont, rt_cont$Response_Time > rt_mean - 2.5 * rt_stdv)
rt_cont = subset(rt_cont, rt_cont$Response_Time < rt_mean + 2.5 * rt_stdv)

rt_simon = rt_df[which(rt_df$Condition == 'Simon'),]
rt_mean = mean(rt_cont$Response_Time)
rt_stdv = sd(rt_cont$Response_Time)
rt_simon = subset(rt_simon, rt_simon$Response_Time > rt_mean - 2.5 * rt_stdv)
rt_simon = subset(rt_simon, rt_simon$Response_Time < rt_mean + 2.5 * rt_stdv)

rt_flan = rt_df[which(rt_df$Condition == 'Flanker'),]
rt_mean = mean(rt_cont$Response_Time)
rt_stdv = sd(rt_cont$Response_Time)
rt_flan = subset(rt_flan, rt_flan$Response_Time > rt_mean - 2.5 * rt_stdv)
rt_flan = subset(rt_flan, rt_flan$Response_Time < rt_mean + 2.5 * rt_stdv)

rt_ms = rt_df[which(rt_df$Condition == 'MultiSource'),]
rt_mean = mean(rt_cont$Response_Time)
rt_stdv = sd(rt_cont$Response_Time)
rt_ms = subset(rt_ms, rt_ms$Response_Time > rt_mean - 2.5 * rt_stdv)
rt_ms = subset(rt_ms, rt_ms$Response_Time < rt_mean + 2.5 * rt_stdv)

shapiro.test(rt_cont$Response_Time)
shapiro.test(rt_simon$Response_Time)
shapiro.test(rt_flan$Response_Time)
shapiro.test(rt_ms$Response_Time)

rt_df = rbind(rt_cont, rt_simon, rt_flan, rt_ms)

rt_df$Condition <- factor(rt_df$Condition, 
                          levels = c('Control', 'Simon', 'Flanker', 'MultiSource'),
                          labels = c('Control', 'Simon', 'Flanker', 'MultiSource'))

ggboxplot(rt_df, x='Group', y='Response_Time', color='Condition') + labs(y='Response Time (ms)')

###############################################################
# Accuracy Outliers by Condition

df = read.csv('E:/Data/MSIT_MIND/accuracy_129.csv')














