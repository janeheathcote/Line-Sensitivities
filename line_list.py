import re
import numpy as np
import matplotlib.pyplot as plt
import os

        

# file & directory setup
directory = r'C:/Users/janee/Documents/Astrophotonics/Methane_Mapping/Input/'
filename = r'CO2_list_2.out'
molecule_name = os.path.basename(filename).split('_')[0]
path = directory + filename
with open(path, 'r') as file:
    lines = file.readlines()


# var 1: molecule ID (6 = CH4)
# var 3: isotopologue ID
# var 4: wavenumber (nu)
# var 5: intensity
# var 6: einstein A coefficient

molec_id, iso_id, nu, intensity, A, sensitivity = [], [], [], [], [], []

for line in lines:
    if line.strip(): 
        
        columns = line.split()
        molec_id.append(int(columns[0]))
        iso_id.append(int(columns[2]))
        nu.append(float(columns[3]))
        intensity.append(float(columns[4]))
        A.append(float(columns[5]))

# conversion from wavenumber (cm^-1) to wavelength (nm)     
wavelength = [(1/val)*1e7 for val in nu]
wavelength = wavelength[::-1]

# flip orders for all other variables
molec_id = molec_id[::-1]
iso_id = iso_id[::-1]
nu = nu[::-1]
intensity = intensity[::-1]
A = A[::-1]


# specify the wavelength range you want to analyze
start_wavelength = 2200
end_wavelength = 2400
bin_size = 0.5

# filter data based on the specified wavelength range
filtered_data = [(wl, inten) for wl, inten in zip(wavelength, intensity) if start_wavelength <= wl <= end_wavelength]
filtered_wavelength, filtered_intensity = zip(*filtered_data) if filtered_data else ([], [])

# define bins
bins = np.arange(start_wavelength, end_wavelength + bin_size, bin_size)  # bins from start_wavelength to end_wavelength (inclusive)
binned_sensitivity = np.zeros(len(bins) - 1)
binned_line_counts = np.zeros(len(bins) - 1, dtype=int)  # array to count number of lines per bin

# accumulate intensities into the bins
for wl, inten in zip(filtered_wavelength, filtered_intensity):
    bin_index = np.digitize(wl, bins) - 1
    if 0 <= bin_index < len(binned_sensitivity):
        binned_sensitivity[bin_index] += inten
        binned_line_counts[bin_index] += 1  # increment line count for this bin

# calculate midpoints of the bins
midpoints = (bins[:-1] + bins[1:]) / 2



def find_sensitivity_peaks(num):
    """
    Prints the bins with the highest sensitivites
    
    Parameters:
        num (int): the number of bins you want to investigate 
    """
    top_bins = sorted(zip(midpoints, binned_sensitivity), key=lambda x: x[1], reverse=True)[:num]

    print("\nTop bins with highest sensitivities:")
    for wavelength, sens in top_bins:
        print(f"Wavelength (nm): {wavelength:.2f}, Sensitivity: {sens:.2e}")
    print('\n')


def find_line_number(start, end):
    """
    Prints the number of lines included in a given bin
    
    Parameters:
        start (float): the starting wavelength of the bin (nm)
        end (float): the ending wavelength of the bin (nm)
    """
    bin_indices = np.where((midpoints >= start) & (midpoints < end))
    if bin_indices[0].size > 0:
        specific_bin_count = binned_line_counts[bin_indices]
        print(f"Number of lines in the bin {start}-{end} nm: {specific_bin_count[0]}")
    else:
        print(f"No data in the bin range {start}-{end} nm.")
        



# plot sensitivity vs. wavelength
plot_sensitivity = False

if plot_sensitivity:
    plt.figure(figsize=(10, 6))
    plt.plot(midpoints, binned_sensitivity)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Sensitivity')
    plt.title('{} Sensitivity vs. Wavelength ({}-{} nm)'.format(molecule_name, start_wavelength, end_wavelength))
    plt.grid(True)
    plt.show()
    
    
plt.figure(figsize=(10, 6))
plt.plot(wavelength, intensity)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity')
plt.title('{} Intensity vs. Wavelength')
plt.grid(True)
plt.show()