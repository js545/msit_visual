# MSIT MIND Manuscript
# Load packages
library(emmeans)
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
rt_mean = mean(rt_simon$Response_Time)
rt_stdv = sd(rt_simon$Response_Time)
rt_simon = subset(rt_simon, rt_simon$Response_Time > rt_mean - 2.5 * rt_stdv)
rt_simon = subset(rt_simon, rt_simon$Response_Time < rt_mean + 2.5 * rt_stdv)

rt_flan = rt_df[which(rt_df$Condition == 'Flanker'),]
rt_mean = mean(rt_flan$Response_Time)
rt_stdv = sd(rt_flan$Response_Time)
rt_flan = subset(rt_flan, rt_flan$Response_Time > rt_mean - 2.5 * rt_stdv)
rt_flan = subset(rt_flan, rt_flan$Response_Time < rt_mean + 2.5 * rt_stdv)

rt_ms = rt_df[which(rt_df$Condition == 'MultiSource'),]
rt_mean = mean(rt_ms$Response_Time)
rt_stdv = sd(rt_ms$Response_Time)
rt_ms = subset(rt_ms, rt_ms$Response_Time > rt_mean - 2.5 * rt_stdv)
rt_ms = subset(rt_ms, rt_ms$Response_Time < rt_mean + 2.5 * rt_stdv)

# Response Time Normality by Group and Condition

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
# Average RT by group and condition
rt_df = read.csv('E:/Data/MSIT_MIND/group_RT_129_modified_v2.csv')
df_con = rt_df[which(rt_df$Group == 'Control'),]
df_hiv = rt_df[which(rt_df$Group == 'HIV'),]

rt_cont = df_con[which(df_con$Condition == 'Control'),]
rt_mean = mean(rt_cont$Response_Time)
rt_stdv = sd(rt_cont$Response_Time)

rt_cont = df_hiv[which(df_hiv$Condition == 'Control'),]
rt_mean = mean(rt_cont$Response_Time)
rt_stdv = sd(rt_cont$Response_Time)

rt_simon = df_con[which(df_con$Condition == 'Simon'),]
rt_mean = mean(rt_simon$Response_Time)
rt_stdv = sd(rt_simon$Response_Time)

rt_simon = df_hiv[which(df_hiv$Condition == 'Simon'),]
rt_mean = mean(rt_simon$Response_Time)
rt_stdv = sd(rt_simon$Response_Time)

rt_flanker = df_con[which(df_con$Condition == 'Flanker'),]
rt_mean = mean(rt_flanker$Response_Time)
rt_stdv = sd(rt_flanker$Response_Time)

rt_flanker = df_hiv[which(df_hiv$Condition == 'Flanker'),]
rt_mean = mean(rt_flanker$Response_Time)
rt_stdv = sd(rt_flanker$Response_Time)

rt_ms = df_con[which(df_con$Condition == 'MultiSource'),]
rt_mean = mean(rt_ms$Response_Time)
rt_stdv = sd(rt_ms$Response_Time)

rt_ms = df_hiv[which(df_hiv$Condition == 'MultiSource'),]
rt_mean = mean(rt_ms$Response_Time)
rt_stdv = sd(rt_ms$Response_Time)

###############################################################
# Average accuracy by group and condition

df = read.csv('E:/Data/MSIT_MIND/accuracy_129.csv')

df_con = df[which(df$Group == 'Control'),]
df_hiv = df[which(df$Group == 'HIV'),]

rt_con_mean = mean(df_con$Control_Accuracy)
rt_con_stdv = sd(df_con$Control_Accuracy)
rt_sim_mean = mean(df_con$Simon_Accuracy)
rt_sim_stdv = sd(df_con$Simon_Accuracy)
rt_flanker_mean = mean(df_con$Flanker_Accuracy)
rt_flanker_stdv = sd(df_con$Flanker_Accuracy)
rt_ms_mean = mean(df_con$MultiSource_Accuracy)
rt_ms_stdv = sd(df_con$MultiSource_Accuracy)

rt_con_mean = mean(df_hiv$Control_Accuracy)
rt_con_stdv = sd(df_hiv$Control_Accuracy)
rt_sim_mean = mean(df_hiv$Simon_Accuracy)
rt_sim_stdv = sd(df_hiv$Simon_Accuracy)
rt_flanker_mean = mean(df_hiv$Flanker_Accuracy)
rt_flanker_stdv = sd(df_hiv$Flanker_Accuracy)
rt_ms_mean = mean(df_hiv$MultiSource_Accuracy)
rt_ms_stdv = sd(df_hiv$MultiSource_Accuracy)

###############################################################
# Accuracy Normality by Condition

df = read.csv('E:/Data/MSIT_MIND/accuracy_129.csv')

shapiro.test(df$Control_Accuracy)
shapiro.test(df$Simon_Accuracy)
shapiro.test(df$Flanker_Accuracy)
shapiro.test(df$MultiSource_Accuracy)

###############################################################
# Plot for Response Times

# Reorganize data for plot
rt_df = read.csv('E:/Data/MSIT_MIND/group_RT_129_modified_v2.csv')
rt_df$Condition <- factor(rt_df$Condition, 
                          levels = c('Control', 'Simon', 'Flanker', 'MultiSource'),
                          labels = c('Control', 'Simon', 'Flanker', 'MultiSource'))
tiff('E:/Data/MSIT_MIND/Manuscript/v1/Figures/response_times.png', units='in', width=5, height=5, res=300)
ggboxplot(rt_df, x='Group', y='Response_Time', color='Condition')
dev.off()

#############################################################################################






