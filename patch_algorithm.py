#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 09:54:00 2021

@author: coralie
"""

import numpy as np
import xarray as xr

filepath = '/home/coralie/Documents/Project_work/Remote_sensing/Deforestation_image_differencing/'
filepathHansen = '/home/coralie/Documents/Project_work/Remote_sensing/Data/Hansen/'  
filepathEVI = '/home/coralie/Documents/Project_work/Remote_sensing/Data/MODIS_EVI/'  

latmin, latmax, lonmin, lonmax = 2.0, 3.6, 20.8, 23.4
StandardNomenclature = str(lonmin) + '-' + str(lonmax) + '_lat' + str(latmin) + '-' + str(latmax)

ForestLossXr = xr.open_dataarray(filepathHansen + 'Processed/ForestLoss2000-19 lon20.8-23.4_lat2.0-3.6.nc')

PercentArray = np.load(filepathHansen + 'Processed/TreeCoverPercentArray_' + str(StandardNomenclature) + '.npy')
#start iterating at 0,0

numOfClump = 0
cellsToLookAt = []
def_thresh = 60

mask = np.zeros_like(PercentArray[0, :, :], dtype = 'uint64')

for k in range(1,20):
    data = PercentArray[k, :, :]
    numOfRows = np.shape(data)[0]
    numOfCols = np.shape(data)[1]
    #create a matching 2d array of -1s to track where we have checked
    row = 0
    col = 0
    mask2 = np.zeros_like(PercentArray[0, :, :], dtype = 'uint64')
  #  print(numOfClump)
    numOfClump = k * 10000
    print(numOfClump)
    while row < numOfRows and col < numOfCols:
      mask2[row][col] = 1
      #this is tests if that cell has been deforested and whether it has already been accounted for
   #   if data[row][col] > 0:
      #    print(data[row][col] )
      if data[row][col] >= def_thresh and mask[row][col] == 0:
      #  print("h")
        numOfClump += 1
        
        #create an empty queue of items to look at
        cellsToLookAt = [[row, col]]        
        #these below 10 lines set the limiters for possible adjacencies - 
        #normally 9, but 4 if a corner and 6 if an edge
        minX, maxX = -1, 2
        minY, maxY = -1, 2
        
        if row == 0:
          minX = 0
        elif row == numOfRows - 1:
          maxX = 1
        
        if col == 0:
          minY = 0
        elif col == numOfCols - 1:
          maxY = 1
        
        #first pass will have the original cell, will continue until empty
        while cellsToLookAt:
          thisItem = cellsToLookAt.pop(0)
          cItemX = thisItem[0]
          cItemY = thisItem[1]
          if cItemY == numOfCols - 1: continue
          if cItemX == numOfRows - 1: continue
          for x in range(minX, maxX):
            for y in range(minY, maxY):
        #      print(cellsToLookAt)
              #test the adjacent 8 - if we enter then one of the surrounding is deforested
              if data[cItemX+x][cItemY+y] >= def_thresh and mask[cItemX+x][cItemY+y] == 0:
                cellsToLookAt.append([cItemX+x, cItemY+y])
                
                mask[cItemX+x][cItemY+y] = numOfClump
    
    #  print("H")       
    #  if col == numOfCols-50: print(row, col)
      if col <= numOfCols - 2:
        col += 1
    
      else:
        row += 1
        col = 0
        
np.save(filepath + 'clump_mask_' + str(def_thresh) + 'd_lon' + StandardNomenclature, mask)

#%%
mask = np.load(filepath + 'clump_mask_' + str(def_thresh) + 'd_lon' + StandardNomenclature + '.npy')

clumps = np.unique(mask)
# removing the 0 from this
clumps = clumps[1:np.shape(clumps)[0]]
print(len(clumps))
count=0
frag_mask = np.zeros([np.shape(mask)[0], np.shape(mask)[1]])
pix_size = 250*250
for clump in clumps:
    edge_count = 0
    clump_ref = np.where(mask == clump)
    print(count)
    count = count + 1
    count_mask = np.zeros([np.shape(mask)[0], np.shape(mask)[1]])
    for k in range(np.shape(clump_ref)[1]):
        lat = clump_ref[0][k]
        lon = clump_ref[1][k]
        if lat != (np.shape(frag_mask)[0] - 1):
            if mask[lat+1, lon] != clump and count_mask[lat+1, lon] == 0:
                edge_count = edge_count+1
                count_mask[lat+1, lon] = 1
        if lon != (np.shape(frag_mask)[1] - 1):
            if mask[lat, lon+1] != clump and count_mask[lat, lon+1] == 0:
                edge_count = edge_count+1        
                count_mask[lat, lon+1] = 1
        if lat != 0:
            if mask[lat-1, lon] != clump and count_mask[lat-1, lon] == 0:
                edge_count = edge_count+1
                count_mask[lat-1, lon] = 1
        if lon != 0:
            if mask[lat, lon-1] != clump and count_mask[lat, lon-1] == 0:
                edge_count = edge_count+1    
                count_mask[lat, lon-1] = 1
    clump_tot = np.shape(clump_ref)[1]
    clump_size = clump_tot * pix_size / 1000000  
    frag_mask[clump_ref] = edge_count / clump_size
np.save(filepath + 'frag_mask_' + str(def_thresh) + 'd_lon' + StandardNomenclature, frag_mask)
                        
                            
                         
#%%
# this section is calculating the area of the clumps
PixelSize = 30*30
clumpSizeStorage = np.zeros([numOfClump], dtype = 'float64')

# clump no = 13,5535 so 1355 times it will print the clump number

for clump in range(1,numOfClump):
    clumpCells = np.size(mask[mask == clump])
    clumpSizeStorage[clump] = clumpCells * PixelSize
    
    if clump % 100 == 0:
        print(clump)
        
# ClumpArr is an array of each clump's size in m2

# setting 0's to nan to not skew the data. I think the zeros are present bc they are cells which had no adjacency
# and so should actually be 900m. however, need to check that this is correct before changing the algorithm

clumpSizeStorage[clumpSizeStorage == 0] = np.nan
# convert m2 to km2
ClumpArr_km = clumpSizeStorage / 1000000
# this step is only neccessary if i haven't redone the clumpSizeStorage bit again (the loop used to start at 0 and
# so it counted zero as a clump)
ClumpArr_km1 = ClumpArr_km[1:]

SizeMean = np.nanmean(ClumpArr_km1)
SizeMed = np.nanmedian(ClumpArr_km1)

SizeSD = np.nanstd(ClumpArr_km1)
Size = np.size(ClumpArr_km1)

SizeSE = SizeSD / np.sqrt(Size)
z_score = 1.96  # this is 95%
CI = SizeSE * z_score
CI1 = SizeMean - CI
CI2 = SizeMean + CI

