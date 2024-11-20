# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 16:06:05 2024

@author: Tom
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits
from General_Code import create_filter_dictionaries

# Call the function from the imported module
filters, filter_errors, hf, ID = create_filter_dictionaries()

plt.style.use('default')

### Stephens Thumbnail Code (edited by me) ###

def make_cutout(data, x, y, width):

    """extracts a cut out from arbitrary data"""

    cutout = np.zeros((width, width))
    
    x = int(np.round(x, 0))
    y = int(np.round(y, 0))
    
    xmin = x - width // 2
    xmax = x + width // 2
    ymin = y - width // 2
    ymax = y + width // 2

    xstart = 0
    ystart = 0
    xend = width
    yend = width

    if xmin < 0:
        xstart = -xmin
        xmin = 0
    if ymin < 0:
        ystart = -ymin
        ymin = 0
    if xmax > data.shape[0]:
        xend -= xmax - data.shape[0]
        xmax = data.shape[0]
    if ymax > data.shape[1]:
        yend -= ymax - data.shape[1]
        ymax = data.shape[1]
    
    if (width % 2) != 0:
        xmax += 1
        ymax += 1
        
    data = np.array(data)
    cutout[xstart:xend,ystart:yend] = data[xmin:xmax,ymin:ymax]

    return cutout

def plot_cutout(ax, im, vmin=None, vmax=None, scaling=False, cmap=cm.magma):
    
    """plots an image cut out"""
    
    if vmin is None:
        vmin = np.min(im)
    if vmax is None:
        vmax = np.max(im)

    if scaling:
        im = scaling(im)

    ax.axis('off')
    ax.imshow(im, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')


if __name__ == "__main__":

    width = 100
    catalogue_file = r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\Python\CEERS_NIRCam1_v0.2.h5"
    
    id_list = [301, 597, 1417, 3917, 3957, 4325, 4797, 6972, 7977, 8728]
    
    #4591, 44, 1057  is edge of imaging
    #1798 possible diffraction spike
    #476, 1663 big diffraction spike
    #2109 very bright object
    #5891 bright in 150
    
    #id_list = [4325, 7977]
    
    # id_list =  [44, 92, 127, 288, 301, 449, 476, 494, 505, 517, 597, 637, 672, 764, 855, 
    #             864, 975, 1057, 1067, 1073, 1417, 1510, 1615, 1645, 1651, 1654, 1663, 1798, 
    #             2042, 2043, 2044, 2093, 2109, 2294, 2471, 3917, 3957, 4325, 4422, 4591, 4797, 
    #             5891, 6253, 6972, 7336, 7977, 8549, 8605, 8728]
    
    image_files = {
        'Hubble.F606W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\Hubble Images\egs_all_acs_wfc_f606w_030mas_v1.9_nircam1_drz.fits",
        'Hubble.F814W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\Hubble Images\egs_all_acs_wfc_f814w_030mas_v1.9_nircam1_drz.fits",    
        'Hubble.F125W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\Hubble Images\egs_all_wfc3_ir_f125w_030mas_v1.9_nircam1_drz.fits",
        'Hubble.F160W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\Hubble Images\egs_all_wfc3_ir_f160w_030mas_v1.9_nircam1_drz.fits",
        'JWST/NIRCam.F115W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\NIRcam Images\ceers_nircam1_f115w_sci_bkgsub_match.fits",
        'JWST/NIRCam.F150W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\NIRcam Images\ceers_nircam1_f150w_sci_bkgsub_match.fits",
        'JWST/NIRCam.F200W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\NIRcam Images\ceers_nircam1_f200w_sci_bkgsub_match.fits",
        'JWST/NIRCam.F277W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\NIRcam Images\ceers_nircam1_f277w_sci_bkgsub_match.fits",
        'JWST/NIRCam.F356W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\NIRcam Images\ceers_nircam1_f356w_sci_bkgsub_match.fits",
        'JWST/NIRCam.F444W': r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\NIRcam Images\ceers_nircam1_f444w_sci_bkgsub.fits",
    
    }

    images = {}
    for image_id, image_file in image_files.items():
        hdu = fits.open(image_file)
        images[image_id] = hdu[0].data

        print(hdu[0].data.shape)

    n_images = len(list(images.keys()))
    print(n_images)

    with h5py.File(catalogue_file, 'r') as hf:

        for id in id_list:

            # determine the x and y coordinate

            s = hf['photom/ID'][()] == id

            x = hf['photom/X'][s][0]
            y = hf['photom/Y'][s][0]

            x,y = y,x


            print(hf['photom/ID'][s][0], x, y)
            
            # Position values for each subplot
            # positions = [(0.606, 0.8), (0.814, 0.8), (1.25, 0.8), (1.60, 0.8), 
            #               (1.15, 0), (1.50, 0), (2.00, 0), (2.77, 0), (3.56, 0), (4.44, 0)]
            #positions = [(x*2, y) for x, y in positions]
            positions = [(0, 0.8), (0.5, 0.8), (1.2, 0.8), (1.9, 0.8),
                          (1, 0), (1.7, 0), (2.4, 0), (2.9, 0), (3.4, 0), (3.9, 0)]

            # initialise figure
            fig, axes = plt.subplots(n_images)
            
            # loop over axes and image_ids
            for image_id, ax, pos in zip(images.keys(), axes, positions):
                print(image_id)
                cutout = make_cutout(images[image_id], x, y, width)
                
                vmin = 0
                vmax = 0.02
            
                print(vmin, vmax)
            
                # Adjust the position of each subplot
                ax.set_position([pos[0], pos[1], 0.7, 0.7])  
                ax.axis('off')
                ax.imshow(cutout, vmin=vmin, vmax=vmax, origin='lower')
                ax.set_title(image_id.split('.')[-1], fontsize=20)
            
            #plt.figure(figsize=(15,15))
            fig.suptitle(f'Object {id}', x=3.5, y=1.4, fontsize=50)

            plt.show()
            
### My SED code ###

            # Define the specific set of filters and galaxies
            selected_filters = ['606', '814', '115', '125', '150', '160', '200', '277', '356', '444']
            
            pivot_wavelength = {'606': 0.606, '814': 0.814, '115': 1.15, '125': 1.25, '150': 1.50,
                                '160': 1.60, '200': 2.00, '277': 2.77, '356': 3.56, '444': 4.44} # pivot wavelength in um
            
            photom_group = hf['photom']  
                        
            # Initialize empty lists to store filter numbers and normalized values
            filter_numbers = []
            normalized_values = []
            normalized_error_values = []
        
            # Filter all datasets based on SN condition and selected galaxy
            for dataset_name in photom_group:
                if dataset_name.startswith('F') and len(dataset_name) == 4:
                    # Extract the filter number (e.g., 'F200' -> '200')
                    filter_number = dataset_name[1:]
                    
                    # Convert filter numbers to wavelengths (microns)
                    if filter_number in pivot_wavelength:
                        wavelength = pivot_wavelength[filter_number]
        
                    if filter_number not in selected_filters:
                        continue  # Skip filters not in the selected set
        
                    # Get the values from the dataset
                    dataset_values = photom_group[dataset_name][:]
        
                    # Filter the dataset values based on selected galaxy
                    filtered_dataset_values = dataset_values[id-1]  # Indexing is 0-based
        
                    # Normalize the filtered dataset values
                    normalized_dataset_values = filtered_dataset_values / filters[277][id-1]  # Indexing is 0-based
        
                    # Append the filter number and normalized values to the respective lists
                    filter_numbers.extend([wavelength])
                    normalized_values.extend([normalized_dataset_values])
                
                if dataset_name.startswith('DF') and len(dataset_name) == 5 and dataset_name != 'DF105' and dataset_name != 'DF140':
                    # Extract the filter number (e.g., 'F200' -> '200')
                    filter_number = dataset_name[2:]
                    
                    # Convert filter numbers to wavelengths (microns)
                    if filter_number in pivot_wavelength:
                        wavelength = pivot_wavelength[filter_number]
        
                    if filter_number not in selected_filters:
                        continue  # Skip filters not in the selected set            
                    
                    # Extract the error dataset corresponding to each filter
                    error_dataset = photom_group[dataset_name][:]
                    
                    # Filter the dataset values based on selected galaxy
                    filtered_error_dataset_values = error_dataset[id-1]  # Indexing is 0-based
        
                    # Normalize the filtered dataset values
                    normalized_error_dataset_values = filtered_error_dataset_values / filters[277][id-1]  # Indexing is 0-based
        
                    # Append the filter number and normalized values to the respective lists
                    normalized_error_values.extend([normalized_error_dataset_values])
            
            # Desired numbers to be filtered
            JWST_filters = [1.15, 1.5, 2.0, 2.77, 3.56, 4.44]
            
            # Initialize lists to store filtered values
            JWST_numbers = []
            JWST_data = []
            JWST_errors = []
            HST_numbers = []
            HST_data = []
            HST_errors = []
            
            # Iterate over filter_numbers and data simultaneously
            for num, value, errors in zip(filter_numbers, normalized_values, normalized_error_values):
                if num in JWST_filters:
                    JWST_numbers.append(num)
                    JWST_data.append(value)
                    JWST_errors.append(errors)
                else:
                    HST_numbers.append(num)
                    HST_data.append(value)
                    HST_errors.append(errors)    
            
            #JWST_yerr = [(np.log10(np.array(JWST_data) + np.array(JWST_errors)) - np.log10(JWST_data)), (np.log10(JWST_data) - np.log10(np.array(JWST_data) - np.array(JWST_errors)))]
            #HST_yerr = [(np.log10(np.array(HST_data) + np.array(HST_errors)) - np.log10(HST_data)), (np.log10(HST_data) - np.log10(np.array(HST_data) - np.array(HST_errors)))]
            
            # Create a scatter plot to display points for the normalized values for each galaxy
            plt.figure(figsize=(6, 5), dpi=300)
            # plt.errorbar(JWST_numbers, np.log10(JWST_data), yerr=JWST_yerr, capsize=5, fmt='o', label='JWST')
            # plt.errorbar(HST_numbers, np.log10(HST_data), yerr=HST_yerr, capsize=5, fmt='o', label='HST')
            plt.errorbar(JWST_numbers, (JWST_data), yerr=(JWST_errors), capsize=5, fmt='o', label='JWST')
            plt.errorbar(HST_numbers, (HST_data), yerr=(HST_errors), capsize=5, fmt='o', label='HST')
            #plt.axhline(0, linestyle='--', linewidth=1)
            
            # Set labels for the axes and a title
            plt.xlabel('Wavelength ($\mu$m)')
            plt.ylabel('Flux (277 Normalized)')
            plt.title(f'Normalized SED for Object {id}')
            plt.ylim(-1,2)
            #plt.yscale('log')#, nonposy='clip')
            plt.grid(alpha=0.5, axis='x')    
            plt.xticks(filter_numbers, filter_numbers, rotation=75)
            plt.legend()
            
            plt.show()