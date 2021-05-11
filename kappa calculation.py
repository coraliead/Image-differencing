#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 11:10:48 2021

@author: coralie
"""

# kappa coefficient calculation!



import numpy as np
import sys

var_ref2 = sys.argv[1]
def_flag = sys.argv[2]
sd_flag = sys.argv[3]
def_thresh = sys.argv[4]
for_thresh = sys.argv[5]
data_quality = sys.argv[6]
edge_pix = sys.argv[7]
patch_size = sys.argv[8]
patch_thresh = sys.argv[9]
yrInQ = sys.argv[10]
yrEnd = sys.argv[11]
if def_flag == '0':
    def_ref = 'def-forest'
    def_ = 'Deforested-forested'
if def_flag == '1':
    def_ref = 'def'
    def_ = 'Deforested'
    
if float(var_ref2) == 0:
    var_ref = 'EVI'
    y1_0, y2_0 = 0.15, 0.9
    y1_1, y2_1 = -0.4, 0.2
    p1, p2 = -.75, .26
elif float(var_ref2) == 1:
    var_ref = 'NDVI'
    y1_0, y2_0 = 0.35, 1.1
    y1_1, y2_1 = -0.4, 0.2
    p1, p2 = -.75, .26

if patch_size == '1':
    patch_size_ref = '250m'
elif patch_size == '2':
    patch_size_ref = '500m'
elif patch_size == '3':
    patch_size_ref = '750m'
elif patch_size == '4':
    patch_size_ref = '1km'
    
yrInQ = int(yrInQ)
yrEnd = int(yrEnd)

fname_def = (var_ref + '_' + def_ref + '_' + sd_flag + '_sd_' + def_thresh + 'd_' + for_thresh +'f_' + data_quality + 'q_' + edge_pix +
    '_' + str(yrInQ) + '-' + str(yrEnd) + '_' + patch_size_ref + '_patches_' + patch_thresh.replace('.', '') )
fname_for = (var_ref + '_' + 'forest' + '_' + sd_flag + '_sd_' + def_thresh + 'd_' + for_thresh +'f_' + data_quality + 'q_' + edge_pix +
    '_' + str(yrInQ) + '-' + str(yrEnd) + '_' + patch_size_ref + '_patches_' + patch_thresh.replace('.', '') )

fpath = def_ + '_' + sd_flag + 'sd_' + def_thresh + 'd_' + for_thresh +'f_' + data_quality + 'q_' + edge_pix + '_' + patch_size_ref + '_' + str(yrInQ) + '-' + str(yrEnd)
filepathImgDiff = '/home/coralie/Documents/Project_work/Remote_sensing/Deforestation_image_differencing/lat2.0-3.6_lon20.8_lon23.8/' + var_ref +'/' + fpath + '/'

def_data = np.loadtxt(filepathImgDiff + 'Data/Deforested/var_all_data_' + fname_def)
for_data = np.loadtxt(filepathImgDiff + 'Data/Forested/var_all_data_' + fname_for)
filler = np.zeros([11,1])
kappa_all = np.zeros([11,1])
kappa_all[kappa_all == 0] = -9999
filler[filler == 0] = -9999
for k in range(0,4):
    
    if k == 0:
        d = 3
        f = 4
    else:
        d = d + 2
        f = f + 3

    kappa = np.zeros([11,3])
    kappa[kappa == 0] = -9999

    # add 3 each time
    kappa[0,0] = for_data[f]
    kappa[2,0] = for_data[14]
    
    # add 4 each time
    kappa[1,1] = def_data[d]
    kappa[2,1] = def_data[10]
    
    kappa[0,1] = kappa[2,1] - kappa[1,1]
    kappa[1,0] = kappa[2,0] - kappa[0, 0]
    kappa[0,2] = kappa[0,0] + kappa[0,1]
    kappa[1,2] = kappa[1,0] + kappa[1,1]
    kappa[2,2] = kappa[0,2] + kappa[1,2]
    
    kappa[3,0] = kappa[0,0]
    kappa[3,1] = kappa[1,1]
    kappa[3,2] = kappa[3,0] + kappa[3,1]
    kappa[4,0] = kappa[2,0] / kappa[2,2]
    kappa[4,1] = kappa[2,1] / kappa[2,2]
    kappa[5,0] = kappa[0,2] / kappa[2,2]
    kappa[5,1] = kappa[1,2] / kappa[2,2]
    kappa[6,0] = kappa[4,0] * kappa[5,0]
    kappa[6,1] = kappa[4,1] * kappa[5,1]
    
    kappa[7,0] = kappa[6,0] + kappa[6,1]
    kappa[8,0] = (kappa[3,2] / kappa[7,0]) / (kappa[2,2] / kappa[7,0])
    kappa[9, 0] = kappa[0,0] / kappa[2,0]
    kappa[10, 0] = kappa[1,1] / kappa[2,1]
    
    y = np.append(kappa, filler, axis = 1)
    kappa_all = np.append(kappa_all, y, axis = 1)
    
np.savetxt(filepathImgDiff + 'Kappa/Kappa_' + fname_def, kappa_all,  delimiter=',')
