#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Yikun Zhang
Last Editing: May 20, 2023

Description: Synthesize the slurm array results (Monte Carlo simulations for 
derivative estimation).
"""

import numpy as np
import pandas as pd

deri_first3 = pd.DataFrame()
deri_first4 = pd.DataFrame()
B = 100
for b in range(1, B+1):
    first3 = pd.read_csv('./deri_res/deri_first_order_'+str(b)+'_sim3_new.csv')
    deri_first3 = pd.concat([deri_first3, first3])
    first4 = pd.read_csv('./deri_res/deri_first_order_'+str(b)+'_sim4_new.csv')
    deri_first4 = pd.concat([deri_first4, first4])
deri_first3.to_csv('deri_first_sim3_new.csv', index=False)
deri_first4.to_csv('deri_first_sim4_new.csv', index=False)


deri_first1 = pd.DataFrame()
B = 100
for b in range(1, B+1):
    first1 = pd.read_csv('./deri_res/deri_first_order_'+str(b)+'_sim1_new.csv')
    deri_first1 = pd.concat([deri_first1, first1])
deri_first1.to_csv('deri_first_sim1_new.csv', index=False)


deri_second6 = pd.DataFrame()
B = 100
for b in range(1, B+1):
    second6 = pd.read_csv('./deri_res/deri_second_order_'+str(b)+'_sim6_new.csv')
    deri_second6 = pd.concat([deri_second6, second6])
deri_second6.to_csv('deri_second_sim6_new.csv', index=False)


deri_second4 = pd.DataFrame()
B = 100
for b in range(1, B+1):
    second4 = pd.read_csv('./deri_res/deri_second_order_'+str(b)+'_sim4_new.csv')
    deri_second4 = pd.concat([deri_second4, second4])
deri_second4.to_csv('deri_second_sim4_new.csv', index=False)

