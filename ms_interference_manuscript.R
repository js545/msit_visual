library(ggplot2)

###############################################################
# RT plot
# Reorganize RT data for plot
rt_df = read.csv('E:/Data/MSIT_MIND/group_RT_129_modified_v2.csv')

rt_df = rt_df[which(rt_df$Condition == 'Control' | rt_df$Condition == 'Simon' | rt_df$Condition == 'Flanker' | rt_df$Condition == 'MultiSource'),]

rt_df$Condition <- factor(rt_df$Condition, 
                          levels = c('Control', 'Simon', 'Flanker', 'MultiSource'),
                          labels = c('Congruent', 'Simon', 'Flanker', 'MultiSource'))
tiff('E:/Data/MSIT_MIND/Manuscript/v1/Figures/response_times.png', units='in', width=5, height=5, res=300)
# ggboxplot(rt_df, x='Group', y='Response_Time', color='Group', ylab='Response Time (ms)')
# ggbarplot(rt_df, x='Condition', y='Response_Time', fill='Group', color='Group',
#           ylab='Response Time (ms)', position = position_dodge(.9), add = 'mean_se', ylim=c(500,900))
ggboxplot(rt_df, x='Condition', y='Response_Time', color='Group', ylab='Response Time (ms)')
dev.off()


###############################################################
# Response Time Shapiro Wilk Test for Normality
rt_df = read.csv('E:/Data/MSIT_MIND/group_RT_129_modified_v2.csv')

rt_cont = rt_df[which(rt_df$Condition == 'Control'),]
rt_simon = rt_df[which(rt_df$Condition == 'Simon'),]
rt_flan = rt_df[which(rt_df$Condition == 'Flanker'),]
rt_ms = rt_df[which(rt_df$Condition == 'MultiSource'),]

# Response Time Normality by Group and Condition

shapiro.test(rt_cont$Response_Time)
shapiro.test(rt_simon$Response_Time)
shapiro.test(rt_flan$Response_Time)
shapiro.test(rt_ms$Response_Time)

###############################################################
# Accuracy Shapiro Wilk Test for Normality

df = read.csv('E:/Data/MSIT_MIND/accuracy_129_friedman.csv')

shapiro.test(df$Control_Accuracy)
shapiro.test(df$Simon_Accuracy)
shapiro.test(df$Flanker_Accuracy)
shapiro.test(df$MultiSource_Accuracy)

# Accuracy Plot

df = read.csv('E:/Data/MSIT_MIND/accuracy_129_friedman.csv')

df = df[which(df$Condition == 'Control' | df$Condition == 'Simon' | df$Condition == 'Flanker' | df$Condition == 'MultiSource'),]

df$Accuracy = df$Accuracy *100

df$Condition <- factor(df$Condition, 
                          levels = c('Control', 'Simon', 'Flanker', 'MultiSource'),
                          labels = c('Congruent', 'Simon', 'Flanker', 'MultiSource'))
# ggbarplot(df, x='Condition', y='Accuracy', fill='Group', color='Group',
          # ylab='Accuracy (%)', position = position_dodge(.9), add = 'mean_se', ylim = c(.8, 1))

tiff('E:/Data/MSIT_MIND/Manuscript/v1/Figures/accuracy.png', units='in', width=5, height=5, res=300)
ggboxplot(df, x='Condition', y='Accuracy', color='Group', ylab='Accuracy (% correct)')
dev.off()


###############################################################
# Additive vs. Multisource Figure

# Additive vs. MultiSource Alpha 2x2 Interaction

df = read.csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/additive_alpha_50_20_-3.csv')
df$Condition <- factor(df$Condition, 
                       levels = c('Additive', 'MultiSource'),
                       labels = c('Additive', 'MultiSource'))
ggbarplot(df, x='Condition', y='pseudo_tvalue', fill='Group', color='Group', position = position_dodge(.9), add='mean_sd')

ggboxplot(df, x='Condition', y='pseudo_tvalue')

ggplot(df, aes(x=pseudo_tvalue, color=Condition)) + 
  geom_histogram(binwidth=1)

###############################################################
# Misc

df = read.csv('E:/Data/MSIT_MIND/accuracy_129_additive_fixed.csv')

df$Condition <- factor(df$Condition, 
                       levels = c('Additive', 'Multisource'),
                       labels = c('Additive', 'Multisource'))

ggboxplot(df, x='Condition', y='Accuracy_int_per', add='jitter', shape='Condition') + scale_y_reverse()





