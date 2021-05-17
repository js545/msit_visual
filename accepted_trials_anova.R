library(openxlsx)

path = 'E:/Data/MSIT_MIND/MEG/Artifact_Scan_LogFiles/'
savename = 'E:/Data/MSIT_MIND/MEG/master_logfile.csv'

filenames_list = list.files(path=path, full.names=TRUE)

df_concat = read.csv(filenames_list[1])

for (i in 2:length(filenames_list)) {
  
  df_temp = read.csv(filenames_list[i])
  df_concat = rbind(df_concat, df_temp)
  
}

mean_all = mean(df_concat$All_Correct_Accepted_Trials)
sd_all = sd(df_concat$All_Correct_Accepted_Trials)
mean_flanker = mean(df_concat$Identity_Accepted_Trials)
sd_flanker = sd(df_concat$Identity_Accepted_Trials)
mean_simon = mean(df_concat$Spatial_Accepted_Trials)
sd_simon = sd(df_concat$Spatial_Accepted_Trials)
mean_control = mean(df_concat$Control_Accepted_Trials)
sd_control = sd(df_concat$Control_Accepted_Trials)
mean_msit = mean(df_concat$Multi_Source_Accepted_Trials)
sd_msit = sd(df_concat$Multi_Source_Accepted_Trials)

threshold = mean_all - 2.5*sd_all
df_concat = subset(df_concat, df_concat$All_Correct_Accepted_Trials > threshold)

# Reorganize data

num_samples = dim(df_concat)[1]

num_correct = c(df_concat$Control_Accepted_Trials, df_concat$Spatial_Accepted_Trials,
                df_concat$Identity_Accepted_Trials, df_concat$Multi_Source_Accepted_Trials)
condition_label = c(rep(c('Control'), num_samples), rep(c('Spatial'), num_samples),
                    rep(c('Identity'), num_samples), rep(c('MSIT'), num_samples))
anova_df = as.data.frame(cbind(num_correct, condition_label))

anova_df[,1] = as.integer(as.character(anova_df[,1]))

# Run ANOVA

res.aov = aov(num_correct ~ condition_label, data=anova_df)
summary(res.aov)




