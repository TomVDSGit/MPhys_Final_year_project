# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 13:31:40 2023

@author: Tom
"""

import h5py

def create_filter_dictionaries():
    hf = h5py.File('CEERS_NIRCam1_v0.2.h5', 'r')
    photom_group = hf['photom']
    ID = hf['photom/ID'][()]

    # Dictionary to store filter datasets
    filters = {}
    filter_errors = {}

    # Filter all datasets based on SN condition
    for dataset_name in photom_group:
        if dataset_name.startswith('F') and len(dataset_name) == 4:
            # Extract the filter number (e.g., 'F200' -> '200')
            filter_number = int(dataset_name[1:])
            # Get the values from the dataset
            dataset_values = photom_group[dataset_name][:]
            
            # Add to the filters dictionary
            filters[filter_number] = dataset_values
            
        elif dataset_name.startswith('DF') and len(dataset_name) == 5:
            # Extract the filter error number (e.g., 'DF200' -> '200')
            filter_error_number = int(dataset_name[2:])
            # Get the values from the dataset
            dataset_values = photom_group[dataset_name][:]
            
            # Add to the filter_errors dictionary
            filter_errors[filter_error_number] = dataset_values

    return filters, filter_errors, hf, ID
