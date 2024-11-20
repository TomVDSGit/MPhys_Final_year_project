# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 09:38:29 2024

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

    width = 150
    catalogue_file = r"C:\Users\Tom\OneDrive\Documents\A.Uni\Year 4\Final Year Project\Python\CEERS_NIRCam1_v0.2.h5"

    id_list = [3892, 813, 2200, 2199, 4100, 1681, 1682, 1683, 1685, 1684]
    
    #id_list = [1684]

    image_files = {
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

            # initialise figure
            fig, axes = plt.subplots(1, n_images, figsize=(n_images*2, 4))
            
            # loop over axes and image_ids
            for image_id, ax in zip(images.keys(), axes):
                print(image_id)
                cutout = make_cutout(images[image_id], x, y, width)
                
                
                # # Assuming 'images' is your array of values and 'image_id' is the index of the image
                # values = images[image_id]
                
                # # Calculate mean and standard deviation
                # mean = np.mean(values)
                # std_dev = np.std(values)
                
                # # Define the lower and upper bounds for filtering outliers
                # lower_bound = mean - 0.5 * std_dev
                # upper_bound = mean + 0.5 * std_dev
                
                # # Filter values within the bounds
                # filtered_values = values[(values >= lower_bound) & (values <= upper_bound)]
                
                # # Calculate min and max from filtered values
                # vmin = np.min(filtered_values)
                # vmax = np.max(filtered_values)
                
                vmin = 0
                vmax = 0.01

                print(vmin, vmax)
                    
                ax.axis('off')
                ax.imshow(cutout, vmin=vmin, vmax=vmax, origin='lower')#, interpolation='none')
                
                ax.set_title(image_id.split('.')[-1])
            
            
            fig.suptitle(f'Object {id}', fontsize=30)
        
            plt.tight_layout()
            plt.show()

#%%
### My SED code ###
            
            # Define the specific set of filters and galaxies
            selected_filters = ['115', '150', '200', '277', '356', '410', '444']
            
            pivot_wavelength = {'606': 0.606, '814': 0.814, '105': 1.05, '115': 1.15, '125': 1.25, '140': 1.40, '150': 1.50,
                                '160': 1.60, '200': 2.00, '277': 2.77, '356': 3.56, '410': 4.10, '444': 4.44} # pivot wavelength in um
            
            photom_group = hf['photom']  
                        
            # Initialize empty lists to store filter numbers and normalized values
            filter_numbers = []
            normalized_values = []
            normalized_error_values = []
        
            # Filter all datasets based on SN condition and selected galaxy
            for dataset_name in photom_group:
                if dataset_name.startswith('F') and len(dataset_name) == 4 and dataset_name != 'F105' and dataset_name != 'F140':
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
                    filtered_dataset_values = dataset_values[id - 1]  # Indexing is 0-based
        
                    # Normalize the filtered dataset values
                    normalized_dataset_values = filtered_dataset_values / filters[277][id - 1]  # Indexing is 0-based
        
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
                    filtered_error_dataset_values = error_dataset[id - 1]  # Indexing is 0-based
        
                    # Normalize the filtered dataset values
                    normalized_error_dataset_values = filtered_error_dataset_values / filters[277][id - 1]  # Indexing is 0-based
        
                    # Append the filter number and normalized values to the respective lists
                    normalized_error_values.extend([normalized_error_dataset_values])
           
            average_values = [5.883691272899731, 12.005726309264817, 12.771619042735326, 1.0,
                              17.82430630001262, 19.748579964051892, 12.407297323167095]

            average_error_values = [1.425064206996843, 1.6353227900326095, 1.2772761116351894, 1.0,
                                    0.9221523692095028, 1.7359660459648385, 1.356512752436968]


            # Create a scatter plot to display points for the normalized values for each galaxy
            plt.figure(figsize=(6, 5), dpi=300)
            plt.errorbar(filter_numbers, normalized_values, yerr=normalized_error_values, capsize=5, fmt='o', label='Values')
            plt.errorbar(filter_numbers, average_values, yerr=average_error_values, label='Average Values', capsize=5, fmt='o', c='C1', alpha=0.5)
            
            # Set labels for the axes and a title
            plt.xlabel('Wavelength ($\mu$m)')
            plt.ylabel('Log Flux (277 Normalized)')
            plt.title(f'Normalized SED for Object {id}')
            #plt.ylim(-15,15)
            plt.yscale('log')
            plt.grid(alpha=0.5, axis='x')    
            plt.xticks(filter_numbers, filter_numbers, rotation=75)
            plt.legend()
            
            plt.show()