import numpy as np
import matplotlib.pyplot as plt
import os

def read_hitran_data(directory, filename):
    """
    Reads a HITRAN data file and extracts relevant variables. Filename should start
    with "[molecule_name]_", where molecule_name is extracted before the first underscore.

    This function assumes:
      - Wavenumber is in column 4 (HITRAN default).
      - Intensity is in column 5 (HITRAN default).
    
    Extracted values:
      - Wavelength (nm) (converted from wavenumber)
      - Intensity of spectral lines
      - Molecule name (from filename)
    """
    molecule_name = os.path.basename(filename).split('_')[0]
    path = os.path.join(directory, filename)
    
    with open(path, 'r') as file:
        lines = file.readlines()
    
    nu, intensity = [], []
    
    for line in lines:
        if line.strip():
            columns = line.split()
            nu.append(float(columns[3]))
            intensity.append(float(columns[4]))
    
    # conversion from wavenumber (cm^-1) to wavelength (nm)     
    wavelength = [(1/val)*1e7 for val in nu][::-1]
    
    # flip intensity to match
    intensity = intensity[::-1]

        
    return molecule_name, wavelength, intensity


def process_hitran_data(molecule_name, wavelength, intensity, start_wavelength, end_wavelength, bin_size=0.5, Plot=True):
    """
    Filters data within a specified wavelength range, bins the intensity values, and 
    optionally plots this sensitivity vs. wavelength.
    """
    
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

    # optionally plot sensitivity vs. wavelength
    if Plot:
        plt.figure(figsize=(10, 6))
        plt.plot(midpoints, binned_sensitivity)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Sensitivity')
        plt.title('{} Sensitivity vs. Wavelength ({}-{} nm)'.format(molecule_name, start_wavelength, end_wavelength))
        plt.grid(True)
        plt.show()
    
    return midpoints, binned_sensitivity, binned_line_counts
        

def find_sensitivity_peaks(midpoints, binned_sensitivity, num=5):
    """
    Prints bins with the highest sensitivites. Default = 5 highest bins.
    """
    
    top_bins = sorted(zip(midpoints, binned_sensitivity), key=lambda x: x[1], reverse=True)[:num]

    print("\nTop bins with highest sensitivities:")
    for wavelength, sens in top_bins:
        print(f"Wavelength (nm): {wavelength:.2f}, Sensitivity: {sens:.2e}")
    print('\n')


def find_line_number(start, end, midpoints, binned_line_counts):
    """
    Prints the number of lines included in a given bin, given the start and end
    wavelength of the bin (in nm).
    """
    
    bin_indices = np.where((midpoints >= start) & (midpoints < end))
    
    if bin_indices[0].size > 0:
        specific_bin_count = binned_line_counts[bin_indices]
        print(f"Number of lines in the bin {start}-{end} nm: {specific_bin_count[0]}")
    else:
        print(f"No data in the bin range {start}-{end} nm.")
        


# Filename & input directory setup
input_dir = r'C:/Users/janee/Documents/Astrophotonics/Methane_Mapping/Input/'
filename = r'CO2_list.out'

# Step 1: Read in HITRAN data file
molecule, wavelength, intensity = read_hitran_data(input_dir, filename)

# Step 2: Bin data and (optionally) plot
start_wavelength = 2000   # [nm]
end_wavelength = 2400
bin_size = 0.05           # [nm]

midpoints, binned_sensitivity, binned_line_counts = process_hitran_data(molecule, wavelength, intensity, start_wavelength, end_wavelength, bin_size, Plot=True)

# Step 3: Print top 5 bins with highest sensitivities
find_sensitivity_peaks(midpoints, binned_sensitivity)

