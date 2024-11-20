# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 12:24:37 2023

@author: Tom
"""

from General_Code import create_filter_dictionaries
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.dpi'] = 500
plt.style.use('default')

# Call the function from the imported module
filters, filter_errors, hf, ID = create_filter_dictionaries()

# Set what filters you want
filter_JWST_num = 150
filter_HST_num = 160
filter_JWST = filters[filter_JWST_num]
filter_err_JWST = filter_errors[filter_JWST_num]
filter_HST = filters[filter_HST_num]
filter_err_HST = filter_errors[filter_HST_num]

#Calculate the ratio between the 2 filters
data_ratio = np.log10((filter_JWST)/(filter_HST))

#Plot the flux of the 2 fliters against each other
plt.figure()
plt.title('HST vs JWST')
plt.scatter(np.log10(filter_JWST), np.log10(filter_HST), s=5)
plt.xlabel('JWST 1.50 $\mu$m ($\log_{10}(nJy)$)')
plt.ylabel('HST 1.60 $\mu$m ($\log_{10}(nJy)$)')
plt.show()

#Plot the errors from both filters
plt.figure()
plt.title('Error on the flux for JWST and HST')
plt.scatter(range(len(filter_err_JWST)), np.log10(filter_err_JWST), s=5, label='JWST')
plt.scatter(range(len(filter_err_HST)), np.log10(filter_err_HST), s=5, label='HST')
plt.ylabel('Flux Error ($\log_{10}(nJy)$)')
plt.xlabel('Object Number')
plt.legend()
plt.show()

#Set all the errors which are above 10^6 to 0 for cleaning
filter_err_JWST[filter_err_JWST > 10**6] = 0
filter_err_HST[filter_err_HST > 10**6] = 0

#Plot the new cleaned data
plt.figure()
plt.title('Cleaned error on the flux for JWST and HST')
plt.scatter(range(len(filter_err_JWST)), np.log10(filter_err_JWST), s=5, label='JWST')
plt.scatter(range(len(filter_err_HST)), np.log10(filter_err_HST), s=5, label='HST')
plt.ylabel('Flux Error ($\log_{10}(nJy)$)')
plt.xlabel('Object Number')
plt.legend()
plt.show()

# Create a mask for values that are below the error or 'nan'
mask_150 = np.logical_or(filter_JWST <= filter_err_JWST, np.isnan(filter_JWST))
mask_160 = np.logical_or(filter_HST <= filter_err_HST, np.isnan(filter_HST))

# Replace the corresponding values with the errors
filter_JWST[mask_150] = filter_err_JWST[mask_150]
filter_HST[mask_160] = filter_err_HST[mask_160]

#Plot the new distribution of flux between the two filters
plt.figure()
plt.title('HST (F160) vs JWST (F150)')
plt.scatter(np.log10(filter_JWST), np.log10(filter_HST), s=5)
plt.xlabel('JWST 1.50 $\mu$m ($\log_{10}(nJy)$)')
plt.ylabel('HST 1.60 $\mu$m ($\log_{10}(nJy)$)')
plt.show()

#Apply mask to data
plt.figure()
plt.title('HST (F160) vs JWST (F150)')
plt.scatter(np.log10(filter_JWST), np.log10(filter_HST), s=5, label='All Data')
plt.scatter(np.log10(filter_JWST[mask_160]), np.log10(filter_HST[mask_160]), s=5, label='Data that is smaller than the noise in HST')
plt.xlabel('JWST 1.50 $\mu$m ($\log_{10}(nJy)$)')
plt.ylabel('HST 1.60 $\mu$m ($\log_{10}(nJy)$)')
plt.legend()
plt.show()

#Zoom into the part of the plot where the objects of interest are
plt.figure()
plt.title('HST (F160) vs JWST (F150)')
plt.scatter(np.log10(filter_JWST), np.log10(filter_HST), s=5, label='All Data')
plt.scatter(np.log10(filter_JWST[mask_160]), np.log10(filter_HST[mask_160]), s=5, label='Data that is smaller than the noise in HST')
plt.xlabel('JWST 1.50 $\mu$m ($\log_{10}(nJy)$)')
plt.ylabel('HST 1.60 $\mu$m ($\log_{10}(nJy)$)')
plt.legend()
plt.xlim(0, 4)
plt.ylim(-1, 3)
plt.show()

#Zoom in further
plt.figure()
plt.title('HST (F160) vs JWST (F150)')
plt.scatter(np.log10(filter_JWST), np.log10(filter_HST), s=5, label='All Data')
plt.scatter(np.log10(filter_JWST[mask_160]), np.log10(filter_HST[mask_160]), s=5, label='Data that is smaller than the noise in HST')
plt.xlabel('JWST 1.50 $\mu$m ($\log_{10}(nJy)$)')
plt.ylabel('HST 1.60 $\mu$m ($\log_{10}(nJy)$)')
plt.legend()
plt.xlim(1.5, 4)
plt.ylim(0.5, 3)
plt.show()

#Filter by objects which are bright in JWST
mask_upper = (np.log10(filter_JWST) >= 1.5)
data_ratio_filt = data_ratio[mask_upper]

# Filter out 'inf' values from 'data_ratio'
inf_mask = ~np.isinf(data_ratio_filt)
data_ratio_inf = data_ratio_filt[inf_mask]

#Plot histogram of the ratio of the 2 filters
plt.figure()
plt.hist(data_ratio_inf, bins=1000, label = f'Number of objects = {len(data_ratio)}')
plt.title('Histogram of JWST vs HST')
plt.xlabel('Ratio of JWST:HST in log scale')
plt.ylabel('Frequency')
plt.legend()
plt.show()

# Create a mask to filter out 'inf' values and values outside the range [1, 5]
data_ratio_mask = ~np.isinf(data_ratio) & (data_ratio >= 1) & (data_ratio <= 5) & (np.log10(filter_JWST) >= 1.5)

# Apply the mask to filter the 'data_ratio' array
data_ratio_filt = data_ratio[data_ratio_mask]

#Plot a histogram with the filtered data
plt.figure()
plt.hist(data_ratio_filt, bins=len(data_ratio_filt), label = f'Number of objects = {len(data_ratio_filt)}')
plt.title('Histogram of JWST vs HST')
plt.xlabel('Ratio of JWST:HST in log scale')
plt.ylabel('Frequency')
plt.legend()
plt.show()

Selected_ID = ID[data_ratio_mask]

print(Selected_ID)