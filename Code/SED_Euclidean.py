# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 14:30:58 2023

@author: Tom
"""

from General_Code import create_filter_dictionaries
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.dpi'] = 500
plt.style.use('default')

# Call the function from the imported module
filters, filter_errors, hf, ID = create_filter_dictionaries()

# Define the specific set of filters and galaxies
selected_filters = ['115', '150', '200', '277', '356', '410', '444']
selected_galaxies = [3892, 813, 2200, 2199, 4100, 1681, 1682, 1683, 1685, 1684]
selected_galaxies = [1683, 1681, 813]
#8646 most similar 
#7272 most disimilar (just the 277 band)

pivot_wavelength = {'606': 0.606, '814': 0.814, '105': 1.05, '115': 1.15, '125': 1.25, '140': 1.40, '150': 1.50,
                    '160': 1.60, '200': 2.00, '277': 2.77, '356': 3.56, '410': 4.10, '444': 4.44} # pivot wavelength in um

photom_group = hf['photom']

average_values = [0.7425718368850297, 0.8683878101368937, 0.9423520053500087, 1.0,
                  0.9362228270542196, 0.8439348130759242, 0.8089566333471095]

for galaxy_num in selected_galaxies:
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
            filtered_dataset_values = dataset_values[galaxy_num - 1]  # Indexing is 0-based

            # Normalize the filtered dataset values
            normalized_dataset_values = filtered_dataset_values / filters[277][galaxy_num - 1]  # Indexing is 0-based

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
            filtered_error_dataset_values = error_dataset[galaxy_num - 1]  # Indexing is 0-based

            # Normalize the filtered dataset values
            normalized_error__dataset_values = filtered_error_dataset_values / filters[277][galaxy_num - 1]  # Indexing is 0-based

            # Append the filter number and normalized values to the respective lists
            normalized_error_values.extend([normalized_error__dataset_values])
            
    filter_numbers = np.sort(filter_numbers)  

    # Create a scatter plot to display points for the normalized values for each galaxy
    plt.figure(figsize=(5, 9))
    plt.subplot(2, 1, 1)  # Creating subplot for the first plot
    #plt.scatter(filter_numbers, np.log(normalized_values), label='Values', s=15)
    plt.errorbar(filter_numbers, normalized_values, yerr=normalized_error_values, label='Values', capsize=5, fmt='o')
    plt.scatter(filter_numbers, average_values, label='Average Values', c='C1')

    # Set labels for the axes and a title
    plt.ylabel('Flux (Normalized)')
    plt.title(f'Normalized SED for Galaxy {galaxy_num}')
    plt.ylim(0.5, 1.5)
    #plt.yscale('log')#, nonposy='mask')
    plt.grid(alpha=0.5, axis='x')    
    plt.xticks(filter_numbers, filter_numbers, rotation=75)
    

    plt.legend()

    normalized_values = np.array(normalized_values)
    average_values = np.array(average_values)

    diff = normalized_values / average_values
    diff_error = diff * np.sqrt((normalized_error_values/normalized_values)**2)
    # Adding another subplot for the second plot
    plt.subplot(2, 1, 2)
    plt.errorbar(filter_numbers, diff, yerr=diff_error, capsize=5, fmt='o', c='C3')
    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('Flux / Average')
    plt.ylim(0.5,2)
    #plt.yscale('log')
    plt.grid(alpha=0.5, axis='x')
    plt.xticks(filter_numbers, filter_numbers, rotation=75)

    # Show the plot
    #plt.tight_layout()
    plt.show()


filter_num = [1.15, 1.5, 2.0, 2.77, 3.56, 4.10, 4.44]

#Average SED Plot
plt.figure()
plt.scatter(filter_num, average_values, c='C1')
plt.ylabel('Flux (277 Normalized)')
plt.title('Normalized Average SED')
plt.xlabel('Wavelength ($\mu$m)')
#plt.yscale('log')#, nonposy='mask')
plt.ylim(0,2)
plt.grid(alpha=0.5, axis='x')    
plt.xticks(filter_numbers, filter_numbers, rotation=75)
plt.show()

