#############################################################################################
library(openxlsx)

demo_df = read.csv('E:/Data/MSIT_MIND/mind_msit_demographics_R_format.csv')

path = 'E:/Data/MSIT_MIND/MEG/Artifact_Scan_LogFiles/'
savename = 'E:/Data/MSIT_MIND/MEG/master_logfile.csv'

filenames_list = list.files(path=path, full.names=TRUE)

df_concat = read.csv(filenames_list[1])

for (i in 2:length(filenames_list)) {
  
  df_temp = read.csv(filenames_list[i])
  df_concat = rbind(df_concat, df_temp)
  
}

colnames(df_concat)[1] = 'MIND.ID.Updated'

df_concat = merge(demo_df, df_concat, by='MIND.ID.Updated')

#############################################################################################
# Correct Accepted Trials

mean_all = mean(df_concat$All_Correct_Total_Trials)
sd_all = sd(df_concat$All_Correct_Total_Trials)
mean_flanker = mean(df_concat$Identity_Total_Trials)
sd_flanker = sd(df_concat$Identity_Total_Trials)
mean_simon = mean(df_concat$Spatial_Total_Trials)
sd_simon = sd(df_concat$Spatial_Total_Trials)
mean_control = mean(df_concat$Control_Total_Trials)
sd_control = sd(df_concat$Control_Total_Trials)
mean_msit = mean(df_concat$Multi_Source_Total_Trials)
sd_msit = sd(df_concat$Multi_Source_Total_Trials)

threshold = mean_all - 2.5*sd_all
df_concat = subset(df_concat, df_concat$All_Correct_Total_Trials > threshold)

# Reorganize data

num_samples = dim(df_concat)[1]

num_correct = c(df_concat$Control_Total_Trials, df_concat$Spatial_Total_Trials,
                df_concat$Identity_Total_Trials, df_concat$Multi_Source_Total_Trials)
condition_label = c(rep(c('Control'), num_samples), rep(c('Spatial'), num_samples),
                    rep(c('Identity'), num_samples), rep(c('MSIT'), num_samples))
anova_df = as.data.frame(cbind(num_correct, condition_label))

anova_df[,1] = as.integer(as.character(anova_df[,1]))

# Run ANOVA

res.aov = aov(num_correct ~ condition_label, data=anova_df)
summary(res.aov)


#############################################################################################
# Response Times

rt_df = read.csv('E:/Data/MSIT_MIND/group_RT_original.csv')

mean_all = mean(rt_df$All_Correct)
sd_all = sd(rt_df$All_Correct)
mean_flanker = mean(rt_df$Flanker)
sd_flanker = sd(rt_df$Flanker)
mean_simon = mean(rt_df$Simon)
sd_simon = sd(rt_df$Simon)
mean_control = mean(rt_df$Control)
sd_control = sd(rt_df$Control)
mean_msit = mean(rt_df$MultiSource)
sd_msit = sd(rt_df$MultiSource)

threshold = mean_all - 2.5*sd_all
rt_df = subset(rt_df, rt_df$All_Correct > threshold)

# Reorganize data

num_samples = dim(rt_df)[1]

RT = c(rt_df$All_Correct, rt_df$Simon,
                rt_df$Flanker, rt_df$MultiSource)
condition_label = c(rep(c('Control'), num_samples), rep(c('Spatial'), num_samples),
                    rep(c('Identity'), num_samples), rep(c('MultiSource'), num_samples))
anova_df = as.data.frame(cbind(RT, condition_label))

anova_df[,1] = as.integer(as.character(anova_df[,1]))

# Run ANOVA

res.aov = aov(num_correct ~ condition_label, data=anova_df)
summary(res.aov)

emt1 = emmeans(res.aov, ~condition_label, var='RT')
test(emt1)
contrast(emt1, 'pairwise')
