library(openxlsx)

path = 'E:/Data/MSIT_MIND/MEG/Artifact_Scan_LogFiles/'
savename = 'E:/Data/MSIT_MIND/MEG/master_logfile.csv'

filenames_list = list.files(path=path, full.names=TRUE)

df_concat = read.csv(filenames_list[1])

for (i in 2:length(filenames_list)) {
  
  df_temp = read.csv(filenames_list[i])
  df_concat = rbind(df_concat, df_temp)
  
}

mean_all_correct = mean(df_concat$All_Correct_Accepted_Trials)
mean_flanker_correct = mean(df_concat$Identity_Accepted_Trials)
mean_simon_correct = mean(df_concat$Spatial_Accepted_Trials)
mean_control_correct = mean(df_concat$Control_Accepted_Trials)
mean_msit_correct = mean(df_concat$Multi_Source_Accepted_Trials)

# Reorganize data