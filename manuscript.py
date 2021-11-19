import numpy as np
import pandas as pd
import pingouin as pg
import matplotlib.pyplot as plt

# Compare RT in a mixed model ANOVA
df = pd.read_csv('E:/Data/MSIT_MIND/group_RT_129_rmanova.csv')
aovrm = pg.mixed_anova(dv='Response_Time', within='Condition', between='Group', subject='MIND_ID', data=df)
print(aovrm.round(3))
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
print(aovrm.round(3))
pg.print_table(aovrm)
ttest = pg.pairwise_ttests(dv='Pseudot_value', within='Condition', subject='ID', padjust='bonf', data=df)
print(ttest[['A', 'B', 'p-corr']])

df = pd.read_csv('E:/Data/MSIT_MIND/VMPs/Visual_4mm/Extracted_Peaks/Alpha/Alpha_Condition_')


