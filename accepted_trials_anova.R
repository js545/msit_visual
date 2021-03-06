#############################################################################################
# Create merged demographic / behavior csv
library(openxlsx)

# demo_df = read.csv('E:/Data/MSIT_MIND/mind_msit_demographics_R_format.csv')
# demo_df = read.csv('E:/Data/MSIT_MIND/MIND_Demographic_Data_v2_R_format.csv')
demo_df = read.csv('E:/Data/MSIT_MIND/msit_mind_demographics_150.csv')

path = 'E:/Data/MSIT_MIND/MEG/Artifact_Scan_LogFiles/'

filenames_list = list.files(path=path, full.names=TRUE)

df_concat = read.csv(filenames_list[1])
df_concat = subset(df_concat)

for (i in 2:length(filenames_list)) {
  
  df_temp = read.csv(filenames_list[i])
  df_temp = subset(df_temp)
  df_concat = rbind(df_concat, df_temp)
  
}

colnames(df_concat)[1] = 'MIND.ID.Updated'

df_concat = merge(demo_df, df_concat, by='MIND.ID.Updated')

# write.csv(df_concat, 'E:/Data/MSIT_MIND/merged_demo_behavioral_150.csv')

df_concat = read.csv('E:/Data/MSIT_MIND/merged_demo_behavioral_129.csv')

#############################################################################################
# Comparisons by HIV Status

df_HIV = df_concat[which(df_concat$Group == 'HIV+ Patient'),]
df_Control = df_concat[which(df_concat$Group == 'Control'),]

#No difference in accepted trials per condition across HIV/Control groups
t.test(df_Control$All_Correct_Accepted_Trials, df_HIV$All_Correct_Accepted_Trials, alternative='two.sided', var.equal=FALSE)
t.test(df_Control$Control_Accepted_Trials, df_HIV$Control_Accepted_Trials, alternative='two.sided', var.equal=FALSE)
t.test(df_Control$Spatial_Accepted_Trials, df_HIV$Spatial_Accepted_Trials, alternative='two.sided', var.equal=FALSE)
t.test(df_Control$Identity_Accepted_Trials, df_HIV$Identity_Accepted_Trials, alternative='two.sided', var.equal=FALSE)
t.test(df_Control$Multi_Source_Accepted_Trials, df_HIV$Multi_Source_Accepted_Trials, alternative='two.sided', var.equal=FALSE)
#############################################################################################
# Comparisons by Correct Accepted Trials

library(emmeans)
library(ggpubr)

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
# summary(res.aov)

emt1 = emmeans(res.aov, ~condition_label, var='num_correct')
test(emt1)
contrast(emt1, 'pairwise')


#############################################################################################
# Comparisons by Response Times

# Need to extract group RT with bad files removed
#  = read.csv('E:/Data/MSIT_MIND/group_RT_150.csv')
rt_df = read.csv('E:/Data/MSIT_MIND/group_RT_129.csv')
# rt_df = read.csv('E:/Data/MSIT_MIND/group_RT_129_modified_v2.csv')

rt_df = rt_df[which(rt_df$Condition == 'Control' | rt_df$Condition == 'Simon' | rt_df$Condition == 'Flanker' | rt_df$Condition == 'MultiSource'),]

# rt_df = rt_df[which(rt_df$Group == 'Control'),]
# rt_df = rt_df[which(rt_df$Group == 'HIV'),]

# shapiro.test(rt_df[which(rt_df$Condition == 'Control'),]$Response_Time)
# shapiro.test(rt_df[which(rt_df$Condition == 'Flanker'),]$Response_Time)
# shapiro.test(rt_df[which(rt_df$Condition == 'Simon'),]$Response_Time)
# shapiro.test(rt_df[which(rt_df$Condition == 'MultiSource'),]$Response_Time)

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

# threshold = mean_all - 2.5*sd_all
# rt_df = subset(rt_df, rt_df$All_Correct > threshold)
# threshold = mean_all + 2.5*sd_all
# rt_df = subset(rt_df, rt_df$All_Correct < threshold)

threshold = mean_control - 2.5*sd_control
rt_df = subset(rt_df, rt_df$Control > threshold)
threshold = mean_control + 2.5*sd_control
rt_df = subset(rt_df, rt_df$Control < threshold)
threshold = mean_simon - 2.5*sd_simon
rt_df = subset(rt_df, rt_df$Simon > threshold)
threshold = mean_simon + 2.5*sd_simon
rt_df = subset(rt_df, rt_df$Simon < threshold)
threshold = mean_flanker - 2.5*sd_flanker
rt_df = subset(rt_df, rt_df$Flanker > threshold)
threshold = mean_flanker + 2.5*sd_flanker
rt_df = subset(rt_df, rt_df$Flanker < threshold)
threshold = mean_msit - 2.5*sd_msit
rt_df = subset(rt_df, rt_df$MultiSource > threshold)
threshold = mean_msit + 2.5*sd_msit
rt_df = subset(rt_df, rt_df$MultiSource < threshold)



# Reorganize data

num_samples = dim(rt_df)[1]

RT = c(rt_df$Control, rt_df$Simon, rt_df$Flanker, rt_df$MultiSource)
condition_label = c(rep(c('Control'), num_samples), rep(c('Spatial'), num_samples),
                    rep(c('Identity'), num_samples), rep(c('MS'), num_samples))
anova_df = as.data.frame(cbind(RT, condition_label))

anova_df$condition_label <- factor(anova_df$condition_label , levels=c("Control", "Spatial", "Identity", "MS"))
anova_df[,1] = as.integer(as.character(anova_df[,1]))
boxplot(RT~condition_label, data=anova_df)

# Reorganize data 2.0 
rt_df$Condition <- factor(rt_df$Condition, 
                       levels = c('Control', 'Simon', 'Flanker', 'MultiSource'),
                       labels = c('Control', 'Simon', 'Flanker', 'MultiSource'))
ggboxplot(rt_df, x='Group', y='RT', color='Condition')



# Run ANOVA

res.aov = aov(RT ~ condition_label, data=anova_df)
summary(res.aov)

emt1 = emmeans(res.aov, ~condition_label, var='RT')
test(emt1)
contrast(emt1, 'pairwise')

#############################################################################################
# Superadditive RT Test

RT = c(rt_df$Superadd_Int, rt_df$MS_Int)
condition_label = c(rep(c('Additive'), num_samples), rep(c('MS'), num_samples))
rt_int = as.data.frame(cbind(RT, condition_label))

rt_int$condition_label <- factor(rt_int$condition_label , levels=c('Additive', 'MS'))
rt_int[,1] = as.integer(as.character(rt_int[,1]))
boxplot(RT~condition_label, data=rt_int)

rt_super = subset(rt_int, rt_int$condition_label == 'Additive')
rt_ms_raw = subset(rt_int, rt_int$condition_label == 'MS')

t.test(rt_super$RT, rt_ms_raw$RT, alternative='two.sided', var.equal=FALSE)


#############################################################################################
# Group ANOVA

df = read.csv('E:/Data/MSIT_MIND/group_anova_129.csv')

df$Condition <- factor(df$Condition, 
                       levels = c('Control', 'Simon', 'Flanker', 'MultiSource'),
                       labels = c('Control', 'Simon', 'Flanker', 'MultiSource'))
head(df)

table(df$Group, df$Condition)

df_Control = df
df_HIV = df_Control[which(df_Control$Group == 'HIV'),]
df_Healthy = df_Control[which(df_Control$Group == 'Control'),]

df_Simon = df[which(df$Condition == 'Control'),]

res.aov2 <- aov(RT ~ Group*Condition, data = df)
summary(res.aov2)

df = read.csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Alpha/Alpha_Interaction_-38_-8_-19_Composite_age_regressed.csv')

res.aov2 <- aov(age_residuals ~ Group*Condition, data=df)
summary(res.aov2)

df = read.csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Gamma/Gamma_Interaction_-26_-61_-43_age_regressed.csv')

res.aov2 <- aov(age_residuals ~ Group*Condition, data=df)
summary(res.aov2)
TukeyHSD(res.aov2)


# Not significant across groups
# Significant across Conditions

library('ggpubr')
ggboxplot(df, x='Condition', y='RT', color='Group')


################
# Temp

# Control
df = read.csv('E:/Data/MSIT_MIND/group_RT_129.csv')
df = df[which(df$Group == 'Control'),]
describe(df$Control)
describe(df$Simon)
describe(df$Flanker)
describe(df$MultiSource)
# HIV
df = read.csv('E:/Data/MSIT_MIND/group_RT_129.csv')
df = df[which(df$Group == 'HIV'),]
describe(df$Control)
describe(df$Simon)
describe(df$Flanker)
describe(df$MultiSource)




