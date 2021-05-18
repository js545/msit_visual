# Jake Son
# 5/18/2021
# Dicon Lab

# Read in file

filename = 'E:/Data/MSIT_MIND/MEG/010/010_mind_unmc_msit_raw_tsss_mc.evt'

df = read.table(filename, sep='\t', header=TRUE)

correct_trial = vector(length = dim(df)[1])

# Label correct trials
for (i in 1:(dim(df)[1]-1)){
  
  if (df$TriNo[i] == 12309 && df$TriNo[i+1] == 12544 ||
      df$TriNo[i] == 12310 && df$TriNo[i+1] == 12800 ||
      df$TriNo[i] == 12311 && df$TriNo[i+1] == 13312){
    
    correct_trial[i] = 'Control'
    
  }
  
  else if(df$TriNo[i] == 12319 && df$TriNo[i+1] == 12544 ||
          df$TriNo[i] == 12320 && df$TriNo[i+1] == 12800 ||
          df$TriNo[i] == 12321 && df$TriNo[i+1] == 13312){
    
    correct_trial[i] = 'Simon'
    
  }
  
  else if(df$TriNo[i] == 12329 && df$TriNo[i+1] == 12544 ||
          df$TriNo[i] == 12330 && df$TriNo[i+1] == 12800 ||
          df$TriNo[i] == 12331 && df$TriNo[i+1] == 13312){
    
    correct_trial[i] = 'Flanker'
    
  }
  
  else if(df$TriNo[i] == 12349 && df$TriNo[i+1] == 12544 ||
          df$TriNo[i] == 12350 && df$TriNo[i+1] == 12800 ||
          df$TriNo[i] == 12351 && df$TriNo[i+1] == 13312){
    
    correct_trial[i] = 'MultiSource'
    
  }
  
  else{
    
    correct_trial[i] = 0
    
  }
  
}

df$correct_trial = correct_trial

# Find lowest correct trial count across conditions
# Does NOT account for artifact scan; thus, if lowest number of trials changes after artifact
# scan, this will not be reflected here.

control_correct = length(grep('Control', df$correct_trial))
flanker_correct = length(grep('Flanker', df$correct_trial))
simon_correct = length(grep('Simon', df$correct_trial))
ms_correct = length(grep('MultiSource', df$correct_trial))

max_count = min(control_correct, flanker_correct, simon_correct, ms_correct)

control_rand = floor(runif((control_correct - max_count), min=1, max=control_correct))
flanker_rand = floor(runif((flanker_correct - max_count), min=1, max=flanker_correct))
simon_rand = floor(runif((simon_correct - max_count), min=1, max=simon_correct))

control_count = 0
flanker_count = 0
simon_count = 0

for (i in 1:length(df$correct_trial)){
  
  if (df$correct_trial[i] == 'Control'){
    
    control_count = control_count+1
    
    if (is.element(control_count, control_rand)){
      
      df$TriNo[i+1] = '1'
      
    }
    
  }
    
  else if (df$correct_trial[i] == 'Flanker'){
    
    flanker_count = flanker_count+1
    
    if (is.element(flanker_count, flanker_rand)){
      
      
      df$TriNo[i+1] = '1'
      
    }
    
  }
    
  else if (df$correct_trial[i] == 'Simon'){
    
    simon_count = simon_count+1
    
    if (is.element(simon_count, simon_rand)){
      
      df$TriNo[i+1] = '1'
      
    }
    
  }
  
}

df = subset(df, select = -correct_trial)

savename = "E:/Data/MSIT_MIND/MEG/010/010_mind_unmc_msit_raw_tsss_mc_replaced.evt"

write.table(df, savename, sep='\t', row.names=FALSE, quote=FALSE)

