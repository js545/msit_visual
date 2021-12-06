import numpy as np
import pandas as pd
import pingouin as pg
import scikit_posthocs as sp
import matplotlib.pyplot as plt

pd.options.display.width = 0

# Compare RT in a mixed model ANOVA
df = pd.read_csv('E:/Data/MSIT_MIND/group_RT_129_rmanova.csv')
aovrm = pg.mixed_anova(dv='Response_Time', within='Condition', between='Group', subject='MIND_ID', data=df)
print(aovrm.round(3))
pg.print_table(aovrm)

# Superadditive Test

df = pd.read_csv('E:/Data/MSIT_MIND/group_RT_129_modified_interference.csv')
aovrm = pg.mixed_anova(dv='Response_Time', within='Condition', between='Group', subject='MIND_ID', data=df)
pg.print_table(aovrm)

df = df[df['Group'] == 'Control']
aovrm = pg.rm_anova(dv='Response_Time', within='Condition', subject='MIND_ID', data=df)
pg.print_table(aovrm)

# Post hoc paired t-test
df = pd.read_csv('E:/Data/MSIT_MIND/group_RT_129_rmanova.csv')
df = df[df['Group'] == 'HIV']
control = df[df['Condition'] == 'Control']
simon = df[df['Condition'] == 'Simon']
flanker = df[df['Condition'] == 'Flanker']
ms = df[df['Condition'] == 'MultiSource']

pd.set_option('display.max_columns', None)
pg.ttest(control['Response_Time'].values, simon['Response_Time'].values, paired=True)
pg.ttest(control['Response_Time'].values, flanker['Response_Time'].values, paired=True)
pg.ttest(control['Response_Time'].values, ms['Response_Time'].values, paired=True)
pg.ttest(simon['Response_Time'].values, flanker['Response_Time'].values, paired=True)
pg.ttest(simon['Response_Time'].values, ms['Response_Time'].values, paired=True)
pg.ttest(flanker['Response_Time'].values, ms['Response_Time'].values, paired=True)


########################################################################################################################
# Alpha Condition 2x4

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Alpha/Alpha_Condition_16_-12_61_Composite_age_regressed.csv')
aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df)
pg.print_table(aovrm)
ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Alpha/Alpha_Condition_-20_-35_68_Composite_age_regressed.csv')
aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df)
pg.print_table(aovrm)
ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Alpha/Alpha_Condition_56_-32_41_Composite_age_regressed.csv')
aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df)
pg.print_table(aovrm)
ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Alpha/Alpha_Interaction_-38_-8_-19_Composite_age_regressed.csv')
aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df)
pg.print_table(aovrm)


########################################################################################################################
# Alpha Interaction 2x4

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Alpha/Alpha_Interaction_-38_-8_-19_Composite_age_regressed.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', between='Group', subject='ID', padjust='bonf', data=df)
print(ttest[['Contrast', 'Condition', 'A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Alpha/Alpha_Interaction_50_20_-4_Composite_age_regressed.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', between='Group', subject='ID', padjust='bonf', data=df)
print(ttest[['Contrast', 'Condition', 'A', 'B', 'T', 'p-corr']])

########################################################################################################################
# Gamma Condition 2x4

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Gamma/Gamma_Condition_18_-89_-1_Composite_age_regressed.csv')
aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df)
pg.print_table(aovrm)
ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Gamma/Gamma_Condition_-58_-17_38_Composite_age_regressed.csv')
aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df)
pg.print_table(aovrm)
ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Gamma/Gamma_Condition_-29_-72_-49_Composite_age_regressed.csv')
aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df)
pg.print_table(aovrm)
ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Gamma/Gamma_Condition_-58_-17_38_Composite.csv')
aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df)
pg.print_table(aovrm)
ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'T', 'p-corr']])

########################################################################################################################
# Gamma Interaction 2x4

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Gamma/Gamma_Interaction_-14_-16_-34_Composite_age_regressed.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)


df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Gamma/Gamma_Interaction_-26_-61_-43_age_regressed.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='Pseudot_value', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

# Gamma Interaction Timeseries

filename = 'E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/gamma_-14_-88_5/timeseries/gamma_-14_-88_5_rel_timeseries_concat.csv'
filename = 'E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/gamma_30_-88_9/timeseries/gamma_30_-88_9_rel_timeseries_concat.csv'

df = pd.read_csv(filename)
df = df[(df['condition'] == 'Control') | (df['condition'] == 'Flanker') | (df['condition'] == 'Simon') | (df['condition'] == 'MultiSource')]

aovrm = pg.friedman(dv='mean', within='condition', subject='ID', data=df, method='chisq')
print(aovrm)
posthoc = sp.posthoc_conover(df, val_col='mean', group_col='condition', p_adjust='bonf')
print(posthoc)

########################################################################################################################
# Alpha Interaction 2x3

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x3Interaction/Alpha_2x3Interaction_-38_-8_-19.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', data=df_control)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x3Interaction/Alpha_2x3Interaction_61_-51_-17.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', data=df_control)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x3Interaction/Alpha_2x3Interaction_50_20_-7.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', data=df_control)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x3Interaction/Alpha_2x3Interaction_28_-14_64.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', data=df_hiv)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x3Interaction/Alpha_2x3Interaction_22_-85_23.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', data=df_control)
print(ttest[['A', 'B', 'T', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x3Interaction/Alpha_2x3Interaction_-45_-75_-14.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', data=df_hiv)
print(ttest[['A', 'B', 'T', 'p-corr']])

########################################################################################################################
# Alpha Interaction 2x2 Additive vs. MultiSource

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x2Interaction/Alpha_2x2Interaction_50_20_-4.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', alternative='less', data=df_hiv)
print(ttest[['Contrast', 'A', 'B', 'T', 'p-unc']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x2Interaction/Alpha_2x2Interaction_-18_-26_-43.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'T', 'p-corr']])

########################################################################################################################
# Alpha Interaction 2x3

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Manuscript_Alpha/2x3Interaction/Alpha_2x3Interaction_-38_-8_-19.csv')

df_control = df[df['Group'] == 'Control']
df_hiv = df[df['Group'] == 'HIV']

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_control)
pg.print_table(aovrm)

aovrm = pg.rm_anova(dv='pseudo_tvalue', within='Condition', subject='ID', data=df_hiv)
pg.print_table(aovrm)

ttest = pg.pairwise_ttests(dv='pseudo_tvalue', within='Condition', subject='ID', padjust='bonf', data=df_control)
print(ttest[['A', 'B', 'T', 'p-corr']])

########################################################################################################################
# Grand Average from Pseudo-t Differences

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Grand_Average/alpha_-30_-76_5_subtraction.csv')

aovrm = pg.mixed_anova(dv='pseudo_tvalue', within='Condition', between='Group', subject='ID', data=df)
pg.print_table(aovrm)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Grand_Average/alpha_34_-76_9_subtraction.csv')

aovrm = pg.mixed_anova(dv='pseudo_tvalue', within='Condition', between='Group', subject='ID', data=df)
pg.print_table(aovrm)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Grand_Average/gamma_-14_-88_5.csv')

aovrm = pg.mixed_anova(dv='pseudo_tvalue', within='Condition', between='Group', subject='ID', data=df)
pg.print_table(aovrm)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Grand_Average/gamma_30_-88_9.csv')

aovrm = pg.mixed_anova(dv='pseudo_tvalue', within='Condition', between='Group', subject='ID', data=df)
pg.print_table(aovrm)

# Grand Average Pseudo-t from Timeseries Differences
# SIGNIFICANT FINDINGS

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/alpha_-30_-76_5/timeseries/alpha_-30_-76_5_rel_timeseries_concat.csv')
df = df[(df['condition'] == 'Simon_Control') | (df['condition'] == 'Flanker_Control') | (df['condition'] == 'MS_Control')]

aovrm = pg.mixed_anova(dv='mean', within='condition', between='group', subject='ID', data=df)
pg.print_table(aovrm)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/alpha_34_-76_9/timeseries/alpha_34_-76_9_rel_timeseries_concat.csv')
df = df[(df['condition'] == 'Simon_Control') | (df['condition'] == 'Flanker_Control') | (df['condition'] == 'MS_Control')]

aovrm = pg.mixed_anova(dv='mean', within='condition', between='group', subject='ID', data=df)
pg.print_table(aovrm)


df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/gamma_-14_-88_5/timeseries/gamma_-14_-88_5_rel_timeseries_concat.csv')
df = df[(df['condition'] == 'Simon_Control') | (df['condition'] == 'Flanker_Control') | (df['condition'] == 'MS_Control')]

aovrm = pg.mixed_anova(dv='mean', within='condition', between='group', subject='ID', data=df)
pg.print_table(aovrm)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/gamma_30_-88_9/timeseries/gamma_30_-88_9_rel_timeseries_concat.csv')
df = df[(df['condition'] == 'Simon_Control') | (df['condition'] == 'Flanker_Control') | (df['condition'] == 'MS_Control')]

aovrm = pg.mixed_anova(dv='mean', within='condition', between='group', subject='ID', data=df)
pg.print_table(aovrm)

########################################################################################################################
# Grand Average 2x2 Additive vs. MultiSource from Pseudo-t Differences

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Grand_Average/alpha_-30_-76_5_subtraction_2x2.csv')

aovrm = pg.mixed_anova(dv='pseudo_tvalue', within='Condition', between='Group', subject='ID', data=df)
pg.print_table(aovrm)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Grand_Average/alpha_34_-76_9_subtraction_2x2.csv')

aovrm = pg.mixed_anova(dv='pseudo_tvalue', within='Condition', between='Group', subject='ID', data=df)
pg.print_table(aovrm)

# Grand Average 2x2 Additive vs. MultiSource from Timeseries Differences

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/alpha_-30_-76_5/timeseries/alpha_-30_-76_5_rel_timeseries_concat_2x2.csv')

aovrm = pg.mixed_anova(dv='mean', within='condition', between='group', subject='ID', data=df)
pg.print_table(aovrm)

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/virutal_sensors/alpha_34_-76_9/timeseries/alpha_34_-76_9_rel_timeseries_concat_2x2.csv')

aovrm = pg.mixed_anova(dv='mean', within='condition', between='group', subject='ID', data=df)
pg.print_table(aovrm)




