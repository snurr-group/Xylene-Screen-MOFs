# plot the Pc versus diff
import pathlib

import numpy as np
import pandas as pd
import os
file_dir = os.getcwd()
Sim_dir = file_dir + '/'
Dataset = pd.read_table(Sim_dir + '/' + "mixCycle_vs_loading.txt", sep = ' ')
Dataset['Total'] = Dataset['pX'] + Dataset['oX'] + Dataset['mX']
Vars=['pX', 'Total']
f = open('GelmanRubin.txt', 'w')
for v in range(0, len(Vars)):
  Variable = Vars[v]
  # only production
  Dataset = Dataset[Dataset['Cycle'] > 11500]
  total_rows = len(Dataset['Cycle']) + 1
  # Look at this website for help
  # https://bookdown.org/rdpeng/advstatcomp/monitoring-convergence.html
  # divide the cycles into 5 blocks
  nblock=5
  L = np.int(total_rows/nblock)
  # pX convergence
  list_of_blocks = []
  list_averages = []
  list_block_var = []

  for a in range(0,nblock):
    list_of_blocks.append(Dataset[a*L: (a+1)*L])
    list_averages.append(np.mean(list_of_blocks[a][Variable]))
    list_block_var.append(np.var(list_of_blocks[a][Variable]))
  # then calculate the r score
  B = np.var(list_averages)*L
  W = np.mean(list_block_var)
  R = ((L-1)/L*W + 1/L*B)/W
  if(R > 1.15):
      string = Variable + " GR is " + str(R) + ", need restart\n"
  else:
      string = Variable + " GR is " + str(R) + ", all good\n"
  f.write(string)
f.close()
