#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 15:33:13 2021

@author: coralie
"""
# i want to do two date image processing using dates when there's EVI data > 80%

import warnings
import xarray as xr
import numpy as np
import iris
import iris.coord_categorisation
import matplotlib.pyplot as plt
import cftime
from datetime import datetime
import rpy2.robjects as ro
from dateutil.relativedelta import relativedelta
from matplotlib import colors
import sys


# var_ref2 = sys.argv[1]
# def_flag = sys.argv[2]
# sd_flag = sys.argv[3]
# def_thresh = sys.argv[4]
# for_thresh = sys.argv[5]
# data_quality = sys.argv[6]
# edge_pix = sys.argv[7]
# patch_size = sys.argv[8]
# patch_thresh = sys.argv[9]
# yrInQ = sys.argv[10]
# yrEnd = sys.argv[11]

var_ref2 = '0'
def_flag = '1'
sd_flag = '1'
def_thresh = '60'
for_thresh = '5'
data_quality = '80'
edge_pix = 'n'
patch_size = '1'
patch_thresh = '0.75'
yrInQ = '2012'
yrEnd = '2013'

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

def_thresh = int(def_thresh)
for_thresh = int(for_thresh)
data_quality = int(data_quality)
patch_size = int(patch_size)
patch_thresh = float(patch_thresh)


filepathEVI = '/home/coralie/Documents/Project_work/Remote_sensing/Data/MODIS_EVI/'  
filepathHansen = '/home/coralie/Documents/Project_work/Remote_sensing/Data/Hansen/'  
filepath = '/home/coralie/Documents/Project_work/Remote_sensing/Deforestation_image_differencing/'

latmin, latmax, lonmin, lonmax = 2.0, 3.6, 20.8, 23.4
# basankusu 
# filepathEVIFig = '/home/coralie/Documents/Project_work/Remote_sensing/Deforestations_impact_on_EVI_BEAST/Basankusu_Case_Study/Figs/'
# latmin, latmax, lonmin, lonmax =  1.05, 1.25, 19.5, 19.9
StandardNomenclature = str(lonmin) + '-' + str(lonmax) + '_lat' + str(latmin) + '-' + str(latmax)

var =  xr.open_dataset(filepathEVI + '/Processed/' +'Proc NDVI Land Mask Applied Aqua lon' + StandardNomenclature +'.nc')
# EVI =  xr.open_dataset(filepathEVI + '/Processed/' +'Proc NDVI Land Mask Applied Aqua lon' + StandardNomenclature +'.nc')
var = var.to_array()
var = var[0,:,:,:]

CumulativeArray = np.load(filepathHansen + 'Processed/TreeCoverCumulativeArray_' + str(StandardNomenclature) + '.npy')
PercentArray = np.load(filepathHansen + 'Processed/TreeCoverPercentArray_' + str(StandardNomenclature) + '.npy')

TreeCoverLinearDeforested = np.load(filepathHansen + 'Processed/TreeCoverLinearDeforested_' + str(StandardNomenclature) + '.npy')
TreeCoverLinearForested = np.load(filepathHansen + 'Processed/TreeCoverLinearForested2000_' + str(StandardNomenclature) + '.npy')

# this section of the program is running through the years and outputting the cell references of cells deforested that 
# year

countHere = 0
totalLat, totalLon = var.lat.values, var.lon.values
patchLat, patchLon = patch_size, patch_size
patchLeeway = (patchLat * patchLon) - ((patchLat * patchLon) * patch_thresh)
if patch_size > 1:
    patch_size = patch_size - 1

ForestMaskAll = np.zeros([yrEnd-yrInQ, np.shape(var)[1], np.shape(var)[2]])

for year in range(yrInQ,yrEnd):   
    print(year)
    ForestMask = np.zeros_like(PercentArray[year-2000,:,:])
    ForestMask[np.isnan(PercentArray[year-2000,:,:]) == True] = np.nan 
   
    CodeCount = 0

    refStore = []
    latRange = range(len(var.lat))
    lonRange = range(len(var.lon))
    pixCount = 0
    # iterating through whole dataset looking for 1km squares of deforestation
    for latCount in latRange[:-int(patch_size):]:
        for lonCount in lonRange[:-int(patch_size):]:     
            # this is to check whether the square i'm looking at has already been categorised as forested or deforested
            # it looks through all pixels in the square, if any have been categorised before then it skips to next 
            # iteration
            if 1 in ForestMask[latCount:latCount + patchLat, lonCount:lonCount + patchLon]: continue
            # extracting the forest loss for looped year for the 1km2 square
            Square = PercentArray[year-2000, latCount:latCount + patchLat, lonCount:lonCount + patchLon]
            warnings.simplefilter('ignore', RuntimeWarning)
                        
            # if enough of a square has higher deforestation than the threshold then its counted as deforested
            if np.size(Square[Square>=def_thresh]) >= (np.size(Square) - patchLeeway):
             #  print("hey")
                # saving the coords so i know exactly which pixels pass the threshold
                
                SquareCoords = np.where(Square>=def_thresh)                   
                ArrCount = 0
                PixRef = []         
                for SqRef in range(len(SquareCoords[0])):
                    # Co1 and Co2 convert the coords on the deforested pixel to the coords of the EVI array
                    Co1 = SquareCoords[0][SqRef] + latCount
                    Co2 = SquareCoords[1][SqRef] + lonCount 
                    ForestMask[Co1, Co2] = Square[SquareCoords[0][SqRef], SquareCoords[1][SqRef]]
                                    
                pixCount = pixCount + 1
                       
            lonCount = lonCount + 1               
        latCount = latCount + 1
    print(countHere)
    TotalForestPix = ForestMask[ForestMask > 0]

    print(str(year)+ ' ' + str(np.shape(TotalForestPix)[0]))

    ForestMaskAll[countHere,:,:] = ForestMask
    countHere = countHere + 1 
    
    
#%%
# now go through the years and the forest mask for each year. compare year before, year of def, and years after. only for EVI > 80%

var =  xr.open_dataset(filepathEVI + '/Processed/' +'Proc ' + var_ref + ' Land Mask Applied Aqua lon' + StandardNomenclature +'.nc')
var = var.to_array()
var = var[0,:,:,:]
yr_arr = np.arange(yrInQ, yrEnd)
bin_arr_all = np.zeros([1,11])
x = np.arange(-2,5)

def_minus_forest_store = np.zeros([7,1])
def_minus_forest_diff_store = np.zeros([6,1])
fig, axs = plt.subplots(2)
fig.subplots_adjust(hspace = 0)
var_store = np.zeros([7,1])
ChangeDetectionMask = np.zeros_like(CumulativeArray[2010-2000,:,:])
sd_arr = np.zeros([3,1])
clump_mask = np.load(filepath + 'clump_mask_' + str(def_thresh) + 'd_lon' + StandardNomenclature + '.npy')
frag_mask = np.load(filepath + 'frag_mask_' + str(def_thresh) + 'd_lon' + StandardNomenclature + '.npy')

pix_size = 250*250
for yr in range(yrInQ,yrEnd):
    forest_mask_yr = ForestMaskAll[int(np.where(yr_arr == yr)[0]),:,:]
    def_loc =np.where(forest_mask_yr > 0)
    
    ForestPixMask = np.zeros_like(CumulativeArray[yr-2000,:,:])
    ForestPixMask[CumulativeArray[yr-2000,:,:] <= for_thresh] = 1
    
    for k in range(0,23):
        print(str(k) + '-' + str(yr) + ', size = ' + str(np.shape(def_loc)[1]) + ' ' + str(var[np.arange(12 + k, 403, 23),:,:].time[0].values)[5:10])
        for m in range(len(def_loc[0])):

            lat = int(def_loc[0][m])
            lon = int(def_loc[1][m])
            varpoint = var[:,lat,lon].values
            clump_no = clump_mask[lat, lon]
            clump_tot = (clump_mask == clump_no).sum()
            clump_size = clump_tot * pix_size / 1000000
            frag = frag_mask[lat, lon]
            PL = int((10/0.25)/2)
            latminPL = lat-PL
            latmaxPL = lat+PL
            lonminPL = lon-PL
            lonmaxPL = lon+PL
            if latminPL < 0: latminPL = 0
            if lonminPL < 0: lonminPL = 0
            # extracting EVI for surrounding 10km and extracting forest mask for this area 
            var10km = var[:, latminPL:latmaxPL, lonminPL:lonmaxPL]
            Forest10kmMask = ForestPixMask[latminPL:latmaxPL, lonminPL:lonmaxPL]
            var10km_timerange = var10km[np.arange(12 + k, 403, 23), :, :]
            
            # looping through each date and selecting the associated var10km variable, using the  forest mask to average all forested 
            # values within the variable and saving this into an array
            var_time_avg = []
            for i in range(len(var10km[:,0,0])):
                var_time = var10km[i,:,:].values
                var_time_avg = np.append(var_time_avg,np.nanmean(var_time[Forest10kmMask == 1]))
           
            def_minus_forest = varpoint - var_time_avg
           
            # recurring var extracts the variable for every year on the same date            
            recurringvar = var[np.arange(12 + k, 403, 23),:,:]
            larger_than, larger_than2 = np.where(recurringvar.time.dt.year >= yr - 2), np.where(recurringvar.time.dt.year >= yr - 1)
            smaller_than, smaller_than2 = np.where(recurringvar.time.dt.year <= yr + 4), np.where(recurringvar.time.dt.year <= yr + 1)
            to_extract = np.intersect1d(larger_than, smaller_than)
            to_extract2 = np.intersect1d(larger_than2, smaller_than2)
           
            nan_flag = 0
            # the year before, the year of and the year after deforestation need to have the surrounding forestry with data availability 
            # over 80%
            for j in to_extract2:
                var_time = var10km_timerange[j,:,:]
                var_time = var_time.values[Forest10kmMask == 1]
                non_nan = np.isnan(var_time)[np.isnan(var_time) == False].size / var_time.size * 100
                if non_nan < data_quality:
                    nan_flag = nan_flag + 1
            #if nan_flag < 1:
            recurringvarpoint1 = varpoint[np.arange(12 + k, 403, 23)]
            recurringvarpoint2 = varpoint[np.arange(11 + k, 403, 23)]
            recurringvarpoint3 = varpoint[np.arange(13 + k, 403, 23)]
            varpoint_plot = recurringvarpoint1[to_extract]
               
            recurring_def_minus_forest_point1 = def_minus_forest[np.arange(12 + k, 403, 23)]
            recurring_def_minus_forest_point2 = def_minus_forest[np.arange(11 + k, 403, 23)]
            recurring_def_minus_forest_point3 = def_minus_forest[np.arange(13 + k, 403, 23)]
            def_minus_forest_point_plot = recurring_def_minus_forest_point1[to_extract]

            if def_flag == '0':
                if sd_flag == '3':
                    if np.shape(recurring_def_minus_forest_point1) != np.shape(recurring_def_minus_forest_point3):
                        recurring_def_minus_forest_point3 = np.append(recurring_def_minus_forest_point3, np.nan)
                    if np.shape(recurring_def_minus_forest_point1) != np.shape(recurring_def_minus_forest_point2):
                        recurring_def_minus_forest_point2 = recurring_def_minus_forest_point2[0:np.shape(recurring_def_minus_forest_point1)[0]]
                    recurring_def_minus_forest_point1 = np.expand_dims(recurring_def_minus_forest_point1, axis = 0)
                    recurring_def_minus_forest_point2 = np.expand_dims(recurring_def_minus_forest_point2, axis = 0)
                    recurring_def_minus_forest_point3 = np.expand_dims(recurring_def_minus_forest_point3, axis = 0)
                    recurring_point = np.append(recurring_def_minus_forest_point1, recurring_def_minus_forest_point2, axis = 0)
                    recurring_point= np.append(recurring_point, recurring_def_minus_forest_point3, axis = 0)
                elif sd_flag == '1':
                    recurring_point = recurring_def_minus_forest_point1
                def_arr = def_minus_forest_point_plot
                
            elif def_flag == '1':
                if sd_flag == '3':
                    if np.shape(recurringvarpoint1) != np.shape(recurringvarpoint3):
                        recurringvarpoint3 = np.append(recurringvarpoint3, np.nan)
                    if np.shape(recurringvarpoint1) != np.shape(recurringvarpoint2):
                        recurringvarpoint2 = recurringvarpoint2[0:np.shape(recurringvarpoint1)[0]]
                    recurringvarpoint1 = np.expand_dims(recurringvarpoint1, axis = 0)
                    recurringvarpoint2 = np.expand_dims(recurringvarpoint2, axis = 0)
                    recurringvarpoint3 = np.expand_dims(recurringvarpoint3, axis = 0)
                    recurring_point = np.append(recurringvarpoint1, recurringvarpoint2, axis = 0)
                    recurring_point = np.append(recurring_point, recurringvarpoint3, axis = 0)
                elif sd_flag == '1':
                    recurring_point = recurringvarpoint1
                def_arr = varpoint_plot
            # this code means that if looking at def-forest then the surrounding 10km forest must have higher quality data, if looking at just
            # deforest then the surrounding 10km forest is irrelevent so does not need to have high quality data
            
            if def_flag == '0' and nan_flag > 0: 
                continue
  
            def_minus_forest_diff = np.zeros(6)
            for g in range(0,6):
                 def_minus_forest_diff[g] = def_arr[g+1] - def_arr[g]
            
            # adding this is so that if either of the differences are a nan then its not included 
            if  np.isnan([def_minus_forest_diff[1],  def_minus_forest_diff[2]]).any() == False:
                 bin_arr = np.zeros([11])
                 if clump_size == 0.0625:
                    bin_arr[8] = 1
                 elif 0.0625 < clump_size <= 0.25:
                    bin_arr[8] = 2
                 elif 0.25 < clump_size <= 0.5625:
                    bin_arr[8] = 3
                 elif clump_size > 0.5625:
                    bin_arr[8] = 4
                 else: print(clump_size)
                
                 if frag < 20:
                    bin_arr[9] = 1
                 elif frag < 40:
                    bin_arr[9] = 2
                 elif frag < 60:
                    bin_arr[9] = 3
                 else: bin_arr[9] = 4
                 
                 if def_thresh < 70:
                    bin_arr[10] = 1
                 elif def_thresh < 80:
                    bin_arr[10] = 2
                 elif def_thresh < 90:
                    bin_arr[10] = 3
                 else: bin_arr[10] = 4
                 if def_minus_forest_diff[1] < def_minus_forest_diff[2]:
                     if def_minus_forest_diff[1] < 0:
                         bin_arr[1] = def_minus_forest_diff[1]
                         
                         if sd_flag == '3':
                             sd_1 = np.nanstd(recurring_point[:, 0:to_extract[1]])
                             bin_arr[5] = np.count_nonzero(~np.isnan(recurring_point[:, 0:to_extract[1]]))
                             bin_arr[7] = sd_1 / np.sqrt(np.count_nonzero(~np.isnan(recurring_point[:, 0:to_extract[1]])))
                             sd_arr = np.append(sd_arr,recurring_point[:, 0:to_extract[1]], axis=1)
                         if sd_flag == '1':
                             sd_1 = np.nanstd(recurring_point[0:to_extract[1]])
                             bin_arr[5] = np.count_nonzero(~np.isnan(recurring_point[0:to_extract[1]]))
                             bin_arr[7] = sd_1 / np.sqrt(np.count_nonzero(~np.isnan(recurring_point[0:to_extract[1]])))
                      #       sd_arr = np.append(sd_arr,recurring_point[0:to_extract[1]], axis=1)
                         
                         sd_half = sd_1 * 0.5
                         sd_2 = sd_1 * 2
                         sd_3 = sd_1 * 3
                         bin_arr[0] = 1
                         bin_arr[6] = sd_1
                         
                         if def_minus_forest_diff[1] <=  -sd_3:
                             bin_arr[2] = 3
                             ChangeDetectionMask[lat,lon] = 2
                         elif def_minus_forest_diff[1] <= -sd_2:
                             bin_arr[2] = 2
                             ChangeDetectionMask[lat,lon] = 2
                         elif def_minus_forest_diff[1] <= -sd_1:
                             bin_arr[2] = 1
                         elif def_minus_forest_diff[1] <= -sd_half:
                             bin_arr[2] = 0.5
                 if def_minus_forest_diff[2] < def_minus_forest_diff[1]:
                     if def_minus_forest_diff[2] < 0:
                         bin_arr[1] = def_minus_forest_diff[2]
                         
                         if sd_flag == '3':
                             sd_1 = np.nanstd(recurring_point[:, 0:to_extract[2]])
                             bin_arr[5] = np.count_nonzero(~np.isnan(recurring_point[:, 0:to_extract[2]]))
                             bin_arr[7] = sd_1 / np.sqrt(np.count_nonzero(~np.isnan(recurring_point[:, 0:to_extract[2]])))
                             sd_arr = np.append(sd_arr,recurring_point[:, 0:to_extract[2]], axis=1)
                         if sd_flag == '1':
                             sd_1 = np.nanstd(recurring_point[0:to_extract[2]])
                             bin_arr[5] = np.count_nonzero(~np.isnan(recurring_point[0:to_extract[2]]))
                             bin_arr[7] = sd_1 / np.sqrt(np.count_nonzero(~np.isnan(recurring_point[0:to_extract[2]])))
                             
                         sd_half = sd_1 * 0.5
                         sd_2 = sd_1 * 2
                         sd_3 = sd_1 * 3
                         bin_arr[0] = 2
                         bin_arr[6] = sd_1
                         
                         if def_minus_forest_diff[2] <= -sd_3:
                             bin_arr[2] = 3
                             ChangeDetectionMask[lat,lon] = 2
                         elif def_minus_forest_diff[2] <= -sd_2:
                             bin_arr[2] = 2
                             ChangeDetectionMask[lat,lon] = 2
                         elif def_minus_forest_diff[2] <= -sd_1:
                             bin_arr[2] = 1
                         elif def_minus_forest_diff[2] <= -sd_half:
                             bin_arr[2] = 0.5
                 if ChangeDetectionMask[lat,lon] != 2:
                    ChangeDetectionMask[lat,lon] = 1
                 bin_arr[3] = k
                 bin_arr[4] = yr
                 bin_arr = np.expand_dims(bin_arr, axis = 0)
                 bin_arr_all = np.append(bin_arr_all, bin_arr, axis = 0 )
#%%
# need to repeat the above but for all forested 
NoChangeDetectionMask = np.zeros_like(CumulativeArray[2010-2000,:,:])
NoChangeDetectionMaskFalse = np.zeros_like(CumulativeArray[2010-2000,:,:])
for_bin_arr_all = np.zeros([6720430,8])
for_count = 0
yr_arr = np.arange(yrInQ, yrEnd)
for yr in range(yrInQ, yrEnd):

    forest_mask_yr = ForestMaskAll[int(np.where(yr_arr == yr)[0]),:,:]
    def_loc =np.where(forest_mask_yr > 0)
    
    ForestPixMask = np.zeros_like(CumulativeArray[yr-2000,:,:])
    ForestPixMask[CumulativeArray[yr-2000,:,:] <= for_thresh] = 1
    
    for k in range(0,23):
        print(str(k) + '-' + str(yr) + ', size = ' + str(np.shape(def_loc)[1]) + ' ' + str(var[np.arange(12 + k, 403, 23),:,:].time[0].values)[5:10])
        # will be used to ensure some forested pixels aren't double counted
        ForestCountMask = np.zeros_like(CumulativeArray[yr-2000,:,:])
        # looping through the def_locs and then using the surrounding 10km forest for experiment
        for m in range(len(def_loc[0])):
            lat = int(def_loc[0][m])
            lon = int(def_loc[1][m])
            varpoint = var[:,lat,lon].values
            
            PL = int((10/0.25)/2)
            latminPL = lat-PL
            latmaxPL = lat+PL
            lonminPL = lon-PL
            lonmaxPL = lon+PL
            if latminPL < 0: latminPL = 0
            if lonminPL < 0: lonminPL = 0

            # extracting EVI for surrounding 10km and extracting forest mask for this area 
            var10km = var[:, latminPL:latmaxPL, lonminPL:lonmaxPL]
            Forest10kmMask = ForestPixMask[latminPL:latmaxPL, lonminPL:lonmaxPL]
            ForestCount10kmMask = ForestCountMask[lat-PL:lat+PL, lon-PL:lon+PL]
            var10km_timerange1 = var10km[np.arange(12 + k, 403, 23), :, :]
            var10km_timerange2 = var10km[np.arange(11 + k, 403, 23), :, :]
            var10km_timerange3 = var10km[np.arange(13 + k, 403, 23), :, :]
            
            lat_ref = np.arange(latminPL, latmaxPL)
            lon_ref = np.arange(lonminPL, lonmaxPL)
            
            larger_than = np.where(var10km_timerange1.time.dt.year >= yr - 1)
            smaller_than = np.where(var10km_timerange1.time.dt.year <= yr + 1)
            to_extract = np.intersect1d(larger_than, smaller_than)
            # extracting year before, year of and yr after deforestation
            var10km_point_plot = var10km_timerange1[to_extract]

            if np.shape(var10km_timerange1) != np.shape(var10km_timerange3):
                y = np.zeros([1, np.shape(var10km_timerange3)[1], np.shape(var10km_timerange3)[2]])
                y[y == 0] = np.nan
                var10km_timerange3 = np.append(var10km_timerange3, y, axis = 0)
            if np.shape(var10km_timerange1) != np.shape(var10km_timerange2):
                var10km_timerange2 = var10km_timerange2[0:np.shape(var10km_timerange1)[0]]
                
            var10km_timerange1 = np.expand_dims(var10km_timerange1, axis = 0)
            var10km_timerange2 = np.expand_dims(var10km_timerange2, axis = 0)
            var10km_timerange3 = np.expand_dims(var10km_timerange3, axis = 0)
            
            if sd_flag == '3':
                recurring_point = np.zeros([3, np.shape(var10km_timerange1)[1], np.shape(var10km_timerange1)[2], np.shape(var10km_timerange1)[3]])
                recurring_point[0,:,:,:] = var10km_timerange1
                recurring_point[1,:,:,:] = var10km_timerange2
                recurring_point[2,:,:,:] = var10km_timerange3
            if sd_flag == '1':
                recurring_point = var10km_timerange1
            
            Forest10km = var10km_point_plot
            ForestNanMask = np.zeros([np.shape(Forest10km)[1], np.shape(Forest10km)[2]])
            # making a forest nan mask which saves the location of every nan across the three time steps
            for d in range(0,3):
                ForestNanMask[np.isnan(Forest10km[d,:,:].values) == True] = 1
            # applying the nan mask to ensure that all three time steps have the same data and nanning out all locations which aren't
            # classed 
            for g in range(0,3):
                Forest10km[g,:,:].values[Forest10kmMask != 1] = np.nan
                Forest10km[g,:,:].values[ForestNanMask == 1] = np.nan
                Forest10km[g,:,:].values[ForestCount10kmMask == 1] = np.nan

            ForestCountMask[lat-PL:lat+PL, lon-PL:lon+PL] = 1
            forest_diff = np.zeros([2, np.shape(Forest10km)[1], np.shape(Forest10km)[2]])
            for g in range(0,2):
                forest_diff[g,:,:] = Forest10km[g+1,:,:] - Forest10km[g,:,:]
            # looping through each point in the forested10km
            for lat_10km in range(0,np.shape(Forest10km)[1]):
                for lon_10km in range(0,np.shape(Forest10km)[2]):
                    point = forest_diff[:, lat_10km, lon_10km]
                    if np.isnan(point).any() == False:
                        forest_bin_arr = np.zeros([8])
                        NoChangeTrue = 0
            # if the drop between the year before and year of deforestation is the largest drop and if its below 0 then:
                        if point[0] < point[1]:
                            if point[0] < 0:
                                forest_bin_arr[1] = point[0]
                                sd_1 = np.nanstd(recurring_point[:,0:to_extract[0], lat_10km, lon_10km])
                                forest_bin_arr[5] = np.count_nonzero(~np.isnan(recurring_point[:,0:to_extract[0], lat_10km, lon_10km]))
                                forest_bin_arr[7] = sd_1 / np.sqrt(np.count_nonzero(~np.isnan(recurring_point[:,0:to_extract[0], lat_10km, lon_10km])))
                                sd_half = sd_1 * 0.5
                                sd_2 = sd_1 * 2
                                sd_3 = sd_1 * 3
                                forest_bin_arr[0] = 1
                                forest_bin_arr[6] = sd_1
                                if point[0] <= -sd_3:
                                    forest_bin_arr[2] = 3
                                    NoChangeDetectionMask[lat_ref[lat_10km],lon_ref[lon_10km]] = NoChangeDetectionMask[lat_ref[lat_10km],lon_ref[lon_10km]] + 1
                                    NoChangeTrue = 1
                                elif point[0] <= -sd_2:
                                    forest_bin_arr[2] = 2
                                    NoChangeDetectionMask[lat_ref[lat_10km],lon_ref[lon_10km]] = NoChangeDetectionMask[lat_ref[lat_10km],lon_ref[lon_10km]] + 1
                                    NoChangeTrue = 1

                                elif point[0] <= -sd_1:
                                    forest_bin_arr[2] = 1
                                elif point[0] <= -sd_half:
                                    forest_bin_arr[2] = 0.5
                        # if the drop between the year of and year after deforestation is the largest drop and if its below 0 then:
                        if point[1] < point[0]:
                            if point[1] < 0:
                                forest_bin_arr[1] = point[1]
                                sd_1 = np.nanstd(recurring_point[:,0:to_extract[1], lat_10km, lon_10km])
                                forest_bin_arr[5] = np.count_nonzero(~np.isnan(recurring_point[:,0:to_extract[1], lat_10km, lon_10km]))
                                forest_bin_arr[7] = sd_1 / np.sqrt(np.count_nonzero(~np.isnan(recurring_point[:,0:to_extract[1], lat_10km, lon_10km])))
                                sd_half = sd_1 * 0.5
                                sd_2 = sd_1 * 2
                                sd_3 = sd_1 * 3
                                forest_bin_arr[0] = 2
                                forest_bin_arr[6] = sd_1
                                # sort this
                                forest_bin_arr[5] = np.shape(recurring_point[:,0:to_extract[1], lat_10km, lon_10km])[0] * np.shape(recurring_point[:,0:to_extract[1], lat_10km, lon_10km])[1]
                                if point[1] <= -sd_3:
                                    forest_bin_arr[2] = 3
                                    NoChangeDetectionMask[lat_ref[lat_10km],lon_ref[lon_10km]] = NoChangeDetectionMask[lat_ref[lat_10km],lon_ref[lon_10km]] + 1
                                    NoChangeTrue = 1

                                elif point[1] <= -sd_2:
                                    forest_bin_arr[2] = 2
                                    NoChangeDetectionMask[lat_ref[lat_10km],lon_ref[lon_10km]] = NoChangeDetectionMask[lat_ref[lat_10km],lon_ref[lon_10km]] + 1
                                    NoChangeTrue = 1

                                elif point[1] <= -sd_1:
                                    forest_bin_arr[2] = 1 
                                elif point[1] <= -sd_half:
                                    forest_bin_arr[2] = 0.5
                    # if either of the points is a nan then it doesn't save into a new line
                        for_count = for_count + 1
                        if NoChangeTrue != 1:
                            NoChangeDetectionMaskFalse[lat_ref[lat_10km],lon_ref[lon_10km]]  = NoChangeDetectionMaskFalse[lat_ref[lat_10km],lon_ref[lon_10km]] + 1
                        forest_bin_arr[3] = k
                        forest_bin_arr[4] = yr
                        forest_bin_arr = np.expand_dims(forest_bin_arr, axis = 0)
                        for_bin_arr_all[for_count,:] = forest_bin_arr

for k in range(0,2):
    # first loop is for forested. second is for deforested
    if k == 0:
        y = np.where(for_bin_arr_all[:,4] == 0)
        # check this bit!!!
        for_bin_arr_all2 = for_bin_arr_all[1:y[0][1], : ]
        data_all = for_bin_arr_all2
        fname = fname_for
        f = 'Forested'
    elif k == 1:
        bin_arr_all1 = bin_arr_all[1:len(bin_arr_all), :]
        data_all = bin_arr_all1
        fname = fname_def
        f = 'Deforested'
        
    minus_diff = np.where(data_all[:,1] < 0)
    sd_half = np.where(data_all[:,2] >= 0.5)
    sd_1 = np.where(data_all[:,2] >= 1)
    sd_2 = np.where(data_all[:,2] >= 2)
    sd_3 = np.where(data_all[:,2] == 3)
    diff_arr = data_all[minus_diff, 1]
    
    total_minus_diff = np.shape(minus_diff)[1] / np.shape(data_all)[0] * 100
    
    total_sdhalf_diff = np.shape(sd_half)[1] / np.shape(data_all)[0] * 100
    total_sd1_diff = np.shape(sd_1)[1] / np.shape(data_all)[0] * 100
    total_sd2_diff = np.shape(sd_2)[1] / np.shape(data_all)[0] * 100
    total_sd3_diff = np.shape(sd_3)[1] / np.shape(data_all)[0] * 100
    sd_val = np.where(data_all[:,6 > 0])
    sd_val_data = data_all[sd_val[0], 6]
    se_val_data = data_all[sd_val[0], 7]
    sample_size = data_all[sd_val[0], 5]
    
    if k == 1:
        var_all_data = [[np.round(total_minus_diff, 2), np.shape(minus_diff)[1], 
                     np.round(total_sdhalf_diff, 2), np.shape(sd_half)[1], 
                     np.round(total_sd1_diff, 2),  np.shape(sd_1)[1], 
                     np.round(total_sd2_diff, 2), np.shape(sd_2)[1] ,
                     np.round(total_sd3_diff, 2),  np.shape(sd_3)[1],
             np.shape(data_all)[0], np.round(np.median(diff_arr), 2), np.round(np.mean(diff_arr), 2), np.round(np.std(diff_arr), 3), 
             np.round(np.std(diff_arr) / np.sqrt(np.shape(data_all)[0]), 4), np.nanmean(sample_size), np.nanmedian(sample_size), np.nanmean(sd_val_data), np.nanmean(se_val_data)]]
    if k == 0:
        var_all_data = [[np.round(total_minus_diff, 2), np.shape(minus_diff)[1], 
                 np.round(total_sdhalf_diff, 2), np.shape(sd_half)[1], np.shape(data_all)[0] - np.shape(sd_half)[1],
                 np.round(total_sd1_diff, 2),  np.shape(sd_1)[1], np.shape(data_all)[0] - np.shape(sd_1)[1],
                 np.round(total_sd2_diff, 2), np.shape(sd_2)[1] ,np.shape(data_all)[0] - np.shape(sd_2)[1],
                 np.round(total_sd3_diff, 2),  np.shape(sd_3)[1],np.shape(data_all)[0] - np.shape(sd_3)[1],
                 np.shape(data_all)[0], np.round(np.median(diff_arr), 2), np.round(np.mean(diff_arr), 2), np.round(np.std(diff_arr), 3), 
                 np.round(np.std(diff_arr) / np.sqrt(np.shape(data_all)[0]), 4), np.nanmean(sample_size), np.nanmedian(sample_size), np.nanmean(sd_val_data), np.nanmean(se_val_data)]]
    month = {}
    months = np.unique(data_all[:,3])
    for m in months:
        m = int(m)
        month[m] = str(var[12 + m, 0, 0].time.values)[5:10]
        
    month_table = np.zeros([23,19])
    for k in months:
        k = int(k)
        data_ref = np.where(data_all[:,3] == k)
        data = data_all[data_ref, :][0]
        minus_diff = np.where(data[:,1] < 0)
        sd_half = np.where(data[:,2] >= 0.5)
        sd_1 = np.where(data[:,2] >= 1)
        sd_2 = np.where(data[:,2] >= 2)
        sd_3 = np.where(data[:,2] == 3)
        diff_arr = data[minus_diff, 1]
        
        sd_val = np.where(data[:,6 > 0])
        sd_val_data = data[sd_val[0], 6]
        se_val_data = data[sd_val[0], 7]
        sample_size = data[sd_val[0], 5]
        
        total_minus_diff = np.shape(minus_diff)[1] / np.shape(data)[0] * 100
        total_sdhalf_diff = np.shape(sd_half)[1] / np.shape(data)[0] * 100
        total_sd1_diff = np.shape(sd_1)[1] / np.shape(data)[0] * 100
        total_sd2_diff = np.shape(sd_2)[1] / np.shape(data)[0] * 100
        total_sd3_diff = np.shape(sd_3)[1] / np.shape(data)[0] * 100
        
        month_table[k,0] =  (np.round(total_minus_diff, 2))
        month_table[k,1] =  np.shape(minus_diff)[1]
        month_table[k,2] =  (np.round(total_sdhalf_diff, 2))
        month_table[k,3] =  np.shape(sd_half)[1]
        month_table[k,4] =  (np.round(total_sd1_diff, 2))
        month_table[k,5] =  np.shape(sd_1)[1]
        month_table[k,6] =  (np.round(total_sd2_diff, 2)) 
        month_table[k,7] =  np.shape(sd_2)[1] 
        month_table[k,8] =  (np.round(total_sd3_diff, 2)) 
        month_table[k,9] =  np.shape(sd_3)[1]
        month_table[k,10] =  (np.shape(data)[0])
        month_table[k,11] =  (np.round(np.nanmedian(diff_arr), 2))
        month_table[k,12] =  (np.round(np.nanmean(diff_arr), 2))
        month_table[k,13] =  (np.round(np.nanstd(diff_arr), 3))
        month_table[k,14] =  (np.round(np.nanstd(diff_arr) / np.sqrt(np.shape(data)[0]), 4))
        month_table[k,15] =  np.nanmean(sample_size)
        month_table[k,16] =  np.nanmedian(sample_size)
        month_table[k,17] =  np.nanmean(sd_val_data)
        month_table[k,18] =  np.nanmean(se_val_data)
        
    yr_table = np.zeros([7,19])
    yr_count = 0
    for k in range(yrInQ, yrEnd):
        data_ref = np.where(data_all[:,4] == k)
        data = data_all[data_ref, :][0]
        minus_diff = np.where(data[:,1] < 0)
        
        if np.shape(minus_diff)[1] != 0:
            sd_half = np.where(data[:,2] >= 0.5)
            sd_1 = np.where(data[:,2] >= 1)
            sd_2 = np.where(data[:,2] >= 2)
            sd_3 = np.where(data[:,2] == 3)
            diff_arr = data[minus_diff, 1]
            
            sd_val = np.where(data[:,6 > 0])
            sd_val_data = data[sd_val[0], 6]
            se_val_data = data[sd_val[0], 7]
            sample_size = data[sd_val[0], 5]
            
            total_minus_diff = np.shape(minus_diff)[1] / np.shape(data)[0] * 100
            total_sdhalf_diff = np.shape(sd_half)[1] / np.shape(data)[0] * 100
            total_sd1_diff = np.shape(sd_1)[1] / np.shape(data)[0] * 100
            total_sd2_diff = np.shape(sd_2)[1] / np.shape(data)[0] * 100
            total_sd3_diff = np.shape(sd_3)[1] / np.shape(data)[0] * 100
            
            yr_table[yr_count,0] = (np.round(total_minus_diff, 2))
            yr_table[yr_count,1] = np.shape(minus_diff)[1]
            yr_table[yr_count,2] = (np.round(total_sdhalf_diff, 2)) 
            yr_table[yr_count,3] = np.shape(sd_half)[1] 
            yr_table[yr_count,4] = (np.round(total_sd1_diff, 2)) 
            yr_table[yr_count,5] = np.shape(sd_1)[1]
            yr_table[yr_count,6] = (np.round(total_sd2_diff, 2))
            yr_table[yr_count,7] = np.shape(sd_2)[1] 
            yr_table[yr_count,8] = (np.round(total_sd3_diff, 2))
            yr_table[yr_count,9] = np.shape(sd_3)[1] 
            yr_table[yr_count,10] = (np.shape(data)[0])
            yr_table[yr_count,11] = (np.round(np.median(diff_arr), 2))
            yr_table[yr_count,12] = (np.round(np.mean(diff_arr), 2))
            yr_table[yr_count,13] = (np.round(np.std(diff_arr), 3))
            yr_table[yr_count,14] = (np.round(np.std(diff_arr) / np.sqrt(np.shape(data)[0]), 4))
            yr_table[yr_count,15] =  np.nanmean(sample_size)
            yr_table[yr_count,16] =  np.nanmedian(sample_size)
            yr_table[yr_count,17] =  np.nanmean(sd_val_data)
            yr_table[yr_count,18] =  np.nanmean(se_val_data)
        yr_count = yr_count + 1
    
    filepathImgDiff = '/home/coralie/Documents/Project_work/Remote_sensing/Deforestation_image_differencing/lat2.0-3.6_lon20.8_lon23.8/' + var_ref +'/' + fpath + '/'
    
    np.savetxt(filepathImgDiff + 'Data/' + f + '/var_all_data_' + fname , np.transpose(var_all_data[0]), delimiter=',')
         #      header='Total pixels with drop (%), Total pixels with a drop (no.), 0.5 sd (%), 0.5 sd (no.), 1 sd (%), 1 sd (no.), 2 sd (%), 2 sd (no.), 3 sd (%), 3 sd (no.), Total, Median drop, Mean drop, Sd drop, SE drop, Mean sample size, Med sample size, Mean SD,Mean SE',
          #     comments = '')
    np.savetxt(filepathImgDiff + 'Data/' + f + '/yr_table_' + fname , yr_table, delimiter=',')
            #   header='Total pixels with drop (%), Total pixels with a drop (no.), 0.5 sd (%), 0.5 sd (no.), 1 sd (%), 1 sd (no.), 2 sd (%), 2 sd (no.), 3 sd (%), 3 sd (no.), Total, Median drop, Mean drop, Sd drop, SE drop, Mean sample size, Med sample size, Mean SD,Mean SE',
            #   comments = '')
    np.savetxt(filepathImgDiff + 'Data/' + f + '/month_table_' + fname , month_table, delimiter=',')
          #     header='Total pixels with drop (%), Total pixels with a drop (no.), 0.5 sd (%), 0.5 sd (no.), 1 sd (%), 1 sd (no.), 2 sd (%), 2 sd (no.), 3 sd (%), 3 sd (no.), Total, Median drop, Mean drop, Sd drop, SE drop, Mean sample size, Med sample size, Mean SD,Mean SE',
           #    comments = '')
#%%

# need to look at detection according to bins. primarily interetsed in 2sd threshold. for month and year the clump total sizes need to be
# calculated again 
clumps = np.unique(clump_mask)
clumps = clumps[1:np.shape(clumps)[0]]

clumpSizeStorage = np.zeros([np.shape(clumps)[0]], dtype = 'float64')

# also should be noted that clump_mask contains all clumps from 2001-19 and therefore needs to be  narrowed down 
bin_arr_all1 = bin_arr_all[1:len(bin_arr_all), :]
data_all = bin_arr_all1
sd_2 = np.where(data_all[:,2] >= 2)


for clump in clumps:
    clumpCells = np.size(clump_mask[clump_mask == clump])
    clump_size = clumpCells * pix_size   
    clumpSizeStorage[clump] = clump_size
    
    
    if clump_size == 0.0625:
        bin_arr[8] = 1
    elif 0.0625 < clump_size <= 0.25:
        bin_arr[8] = 2
    elif 0.25 < clump_size <= 0.5625:
        bin_arr[8] = 3
    elif clump_size > 0.5625:
        bin_arr[8] = 4
    
    


# area bins: needs to be a percentage, so the number of changes detected vs the total number of changes detection for that area