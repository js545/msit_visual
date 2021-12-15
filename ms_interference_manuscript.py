import numpy as np
import pandas as pd
import pingouin as pg
import scikit_posthocs as sp
import matplotlib.pyplot as plt
import scipy

########################################################################################################################
# Summary Stats

df = pd.read_csv('E:/Data/MSIT_MIND/merged_demo_behavioral_129.csv')

########################################################################################################################
# Compare RT in a mixed model ANOVA
df = pd.read_csv('E:/Data/MSIT_MIND/group_RT_129_rmanova.csv')

aovrm = pg.mixed_anova(dv='Response_Time', within='Condition', between='Group', subject='MIND_ID', data=df)
print(aovrm.round(3))
pg.print_table(aovrm)

df_con = df[df['Group'] == 'Control']
ms = df_con[df['Condition'] == 'MultiSource']['Response_Time'].values
flanker = df_con[df['Condition'] == 'Flanker']['Response_Time'].values
simon = df_con[df['Condition'] == 'Simon']['Response_Time'].values
control = df_con[df['Condition'] == 'Control']['Response_Time'].values

ms_con = scipy.stats.ttest_rel(ms, control)
ms_flanker = scipy.stats.ttest_rel(ms, flanker)
ms_simon = scipy.stats.ttest_rel(ms, simon)

flanker_con = scipy.stats.ttest_rel(flanker, control)
flanker_simon = scipy.stats.ttest_rel(flanker, simon)

simon_con = scipy.stats.ttest_rel(simon, control)

df_hiv = df[df['Group'] == 'HIV']
ms = df_hiv[df['Condition'] == 'MultiSource']['Response_Time'].values
flanker = df_hiv[df['Condition'] == 'Flanker']['Response_Time'].values
simon = df_hiv[df['Condition'] == 'Simon']['Response_Time'].values
control = df_hiv[df['Condition'] == 'Control']['Response_Time'].values

ms_con = scipy.stats.ttest_rel(ms, control)
ms_flanker = scipy.stats.ttest_rel(ms, flanker)
ms_simon = scipy.stats.ttest_rel(ms, simon)

flanker_con = scipy.stats.ttest_rel(flanker, control)
flanker_simon = scipy.stats.ttest_rel(flanker, simon)

simon_con = scipy.stats.ttest_rel(simon, control)

########################################################################################################################
# Compare Accuracy in a Friedman's test

df = pd.read_csv('E:/Data/MSIT_MIND/accuracy_129_friedman.csv')

df_con = df[df['Group'] == 'Control']
ftest = scipy.stats.friedmanchisquare(df_con[df_con['Condition'] == 'Control']['Accuracy'].values,
                              df_con[df_con['Condition'] == 'Simon']['Accuracy'].values,
                              df_con[df_con['Condition'] == 'Flanker']['Accuracy'].values,
                              df_con[df_con['Condition'] == 'MultiSource']['Accuracy'].values)

posthoc = sp.posthoc_conover(df_con, val_col='Accuracy', group_col='Condition', p_adjust='bonf')

df_hiv = df[df['Group'] == 'HIV']
ftest = scipy.stats.friedmanchisquare(df_hiv[df_hiv['Condition'] == 'Control']['Accuracy'].values,
                              df_hiv[df_hiv['Condition'] == 'Simon']['Accuracy'].values,
                              df_hiv[df_hiv['Condition'] == 'Flanker']['Accuracy'].values,
                              df_hiv[df_hiv['Condition'] == 'MultiSource']['Accuracy'].values)

posthoc = sp.posthoc_conover(df_hiv, val_col='Accuracy', group_col='Condition', p_adjust='bonf')

mann_whit = scipy.stats.mannwhitneyu(df_con['Accuracy'].values, df_hiv['Accuracy'].values)

########################################################################################################################
# Potential superadditive effects in RT

df = pd.read_csv('E:/Data/MSIT_MIND/group_RT_129_modified_interference.csv')
aovrm = pg.mixed_anova(dv='Response_Time', within='Condition', between='Group', subject='MIND_ID', data=df)
pg.print_table(aovrm)

add = df[df['Condition'] == 'Additive']['Response_Time'].values
ms = df[df['Condition'] == 'MultiSource']['Response_Time'].values

p_ttest = scipy.stats.ttest_rel(add, ms)

# Potential Superadditive in Accuracy

df = pd.read_csv('E:/Data/MSIT_MIND/accuracy_129_additive.csv')

df_con = df[df['Group'] == 'Control']
add = df_con[df_con['Condition'] == 'Additive']['Accuracy_Effect'].values
ms = df_con[df_con['Condition'] == 'Multisource']['Accuracy_Effect'].values

wilc = scipy.stats.wilcoxon(add, ms)

df_hiv = df[df['Group'] == 'HIV']
add = df_hiv[df_hiv['Condition'] == 'Additive']['Accuracy_Effect'].values
ms = df_hiv[df_hiv['Condition'] == 'Multisource']['Accuracy_Effect'].values

wilc = scipy.stats.wilcoxon(add, ms)

mann_whit = scipy.stats.mannwhitneyu(df_con['Accuracy_Effect'].values, df_hiv['Accuracy_Effect'].values)


########################################################################################################################
# Alpha Interaction

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/additive_alpha_50_20_-3.csv')

aovrm = pg.mixed_anova(dv='pseudo_tvalue', within='Condition', between='Group', subject='ID', data=df)
print(aovrm.round(3))
pg.print_table(aovrm)

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', between='Group', subject='ID', padjust='bonf', data=df)
print(ttest[['Contrast', 'Condition', 'A', 'B', 'T', 'p-corr']])

########################################################################################################################
# Plot for Correlation Between Neural and Behavioral Interference Effect

import seaborn as sns

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/MS_interference_correlation/alpha/alpha_-7_29_1.csv')

plt.figure(figsize=(10, 10))
p = sns.lmplot(x='pseudo_tvalue', y='MS_Int', hue='Group', data=df, height=5,
                 palette=dict(Control='#F8766D', HIV='#00BFC4'), truncate=False, ci=False)
p.set(xlabel='Pseudo-t value', ylabel='Multisource Interference (ms)')
plt.savefig('E:/Data/MSIT_MIND/Manuscript/v1/Figures/alpha_-7_29_1_correlation.png', dpi=500)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/MS_interference_correlation/alpha/alpha_-14_12_37.csv')

plt.figure(figsize=(10, 10))
p = sns.lmplot(x='pseudo_tvalue', y='MS_Int', hue='Group', data=df, height=5,
                 palette=dict(Control='#F8766D', HIV='#00BFC4'), truncate=False, ci=False)
p.set(xlabel='Pseudo-t value', ylabel='Multisource Interference (ms)')
plt.savefig('E:/Data/MSIT_MIND/Manuscript/v1/Figures/alpha_-14_12_37_correlation.png', dpi=500)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/MS_interference_correlation/gamma/gamma_-44_9_41.csv')

plt.figure(figsize=(10, 10))
p = sns.lmplot(x='pseudo_tvalue', y='MS_Int', hue='Group', data=df, height=5,
                 palette=dict(Control='#F8766D', HIV='#00BFC4'), truncate=False, ci=False)
p.set(xlabel='Pseudo-t value', ylabel='Multisource Interference (ms)')
plt.savefig('E:/Data/MSIT_MIND/Manuscript/v1/Figures/gamma_-44_9_41_correlation.png', dpi=500)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/MS_interference_correlation/gamma/gamma_17_-62_23.csv')

plt.figure(figsize=(10, 10))
p = sns.lmplot(x='pseudo_tvalue', y='MS_Int', hue='Group', data=df, height=5,
                 palette=dict(Control='#F8766D', HIV='#00BFC4'), truncate=False, ci=False)
p.set(xlabel='Pseudo-t value', ylabel='Multisource Interference (ms)')
plt.savefig('E:/Data/MSIT_MIND/Manuscript/v1/Figures/gamma_17_-62_23_correlation.png', dpi=500)


