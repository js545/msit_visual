# Load data

df = read.csv('E:/Data/MSIT_MIND/mind_msit_demographics_R_format.csv')

df_HIV = df[which(df$Group == 'HIV+ Patient'),]
df_Control = df[which(df$Group == 'Control'),]

t.test(df_Control$Age, df_HIV$Age, alternative='two.sided', var.equal=FALSE)
t.test(df_Control$Education, df_HIV$Education, alternative='two.sided', var.equal=FALSE)

# table(df$Group, df$Dominant.Hand)

sex_chi = chisq.test(table(df$Group, df$Sex.at.Birth))
hand_chi = chisq.test(table(df$Group, df$Dominant.Hand))
race_chi = chisq.test(table(df$Group, df$Race))
ethni_chi = chisq.test(table(df$Group, df$Ethnicity))

# All but Education are matched across the two groups