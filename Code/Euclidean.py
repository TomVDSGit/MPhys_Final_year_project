# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 18:52:19 2024

@author: Tom
"""

from General_Code import create_filter_dictionaries
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 500
plt.style.use('default')

# Call the function from the imported module
filters, filter_errors, hf, ID = create_filter_dictionaries()

# Initialize empty lists to store filter numbers and normalized values
filter_numbers = []
average_values = []
average_error_values = []

photom_group = hf['photom'] 

pivot_wavelength = {'115': 1.15, '150': 1.50, '200': 2.00, '277': 2.77, '356': 3.56, '410': 4.10, '444': 4.44}

# Filter all datasets based on SN condition
for dataset_name in photom_group:
    if dataset_name.startswith('F') and len(dataset_name) == 4: # Extract the filter number (e.g., 'F200' -> '200')
        filter_number = dataset_name[1:]
        # Get the values from the dataset
        dataset_values = photom_group[dataset_name][:]
        
        # Convert filter numbers to wavelengths (microns)
        if filter_number in pivot_wavelength:
            wavelength = pivot_wavelength[filter_number]
        
        if filter_number not in pivot_wavelength:
                continue  # Skip filters not in the selected set

        # Normalize the filtered dataset values
        normalized_dataset_values = np.array(filters[int(filter_number)]) / np.array(filters[277])
        
        #normalized_dataset_values = normalized_dataset_values[(normalized_dataset_values < 1e2) & 
                                                           #  (normalized_dataset_values > -1e2)]
        #print( filters[int(filter_number)].shape)
        
        # Find the average value
        normalized_average_values = np.nanmedian(normalized_dataset_values)
        average_values.append(normalized_average_values)
        
        #print(filter_number, normalized_average_values)

# List of desired filters
JWST_filters = [115, 150, 200,  277, 356, 410, 444]

# Create a list to store arrays of values
values_array = []

# Iterate through the keys of the dictionary and collect values from arrays with desired numbers
for key in sorted(filters.keys()):  # Sort the keys to ensure order
    if key in JWST_filters:
        values = filters[key] / filters[277]
        for i, value in enumerate(values):
            if len(values_array) <= i:
                values_array.append([])
            values_array[i].append(value)
            
values_array = np.array(values_array)

# Create a list to store arrays of values
errors_array = []

# Iterate through the keys of the dictionary and collect values from arrays with desired numbers
for key in sorted(filters.keys()):  # Sort the keys to ensure order
    if key in JWST_filters:
        errors = filter_errors[key] / filters[277]
        for i, value in enumerate(errors):
            if len(errors_array) <= i:
                errors_array.append([])
            errors_array[i].append(value)     
            
errors_array = np.array(errors_array)            
            
# Perform element-wise division
SN = values_array / errors_array

# Create a mask based on the condition
SN_mask = SN > 2

# Use the mask to filter values_array
filtered_values = values_array[SN_mask]

# Apply log to all values in values_array
log_flux_values = []
for values in filtered_values:
    log_flux_values.append(np.log10(values))
    
log_average_values = np.log10(average_values)

# Calculate the Euclidean distance for each array in values_array
euclidean_distances = []

for ID_val, values in zip(ID, log_flux_values):
    euclidean_distance = np.linalg.norm(log_average_values - np.array(values))
    if not np.isnan(euclidean_distance):  # Filter out nan distances
        euclidean_distances.append((ID_val, euclidean_distance))

# Sort the list of tuples based on the Euclidean distances
euclidean_distances.sort(key=lambda x: x[1])

# Extract the first 100 ID values
top_ID = [item[0] for item in euclidean_distances[-10:]]
#top_100_ID.sort()
print(top_ID)


# Extracting the values
first_values = [pair[0] for pair in euclidean_distances]
second_values = [pair[1] for pair in euclidean_distances]

# Filtering the array by the second value being greater than 4
threshold = 4
filtered_distances = np.array(second_values)[np.array(second_values) > threshold]

# Plotting
plt.scatter(first_values, second_values, marker='o', s=5)
plt.axhline(y=threshold, color='r', linestyle='--', label=f'Threshold:\n{len(filtered_distances)} objects > Threshold')
plt.xlabel('ID')
plt.ylabel('Euclidean distance')
plt.title('Euclidean distance for each ID')
plt.legend()
plt.grid(alpha=0.5)
plt.show()

flux_values = []

# Iterate through the ID array and fetch the relevant flux value
for obj_id in first_values:
    flux_value = filters[277][obj_id - 1]  # Assuming IDs start from 1
    flux_values.append(flux_value)

flux_values = np.array(flux_values)
second_values = np.array(second_values)

mask = (flux_values < 1e5) & (flux_values > -1e5)
flux_values = flux_values[mask]
second_values_filt = second_values[mask]

# Plotting
plt.scatter(flux_values, second_values_filt, marker='o', s=5)
plt.xlabel('Flux')
plt.ylabel('Euclidean distance')
plt.title('Euclidean distance for each ID')
plt.grid(alpha=0.5)
#plt.xscale('log')
plt.show()