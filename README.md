# HITRAN Data Processing and Sensitivity Analysis

This repository contains code for processing HITRAN spectral line data and performing sensitivity analysis on various molecules. It includes functions for reading, filtering, and binning HITRAN data, and identifying regions with the highest sensitivity for the specified wavelength range.

Sample input data is provided in the Input/ directory for CH4, CO, CO2, H2O, N2O, and NH3. There are two sub-folders:
- 1200-1700/: Contains spectral data for molecules in the wavelength range of 1200-1700 nm.
- 2000-2500/: Contains spectral data for molecules in the wavelength range of 2000-2500 nm.

Any data file from the [HITRANonline line-by-line search](https://hitran.org/lbl/) should work as an input.

### Example:

The following example reads data for CO2, processes it in the wavelength range of 2200-2400 nm, plots the sensitivity, and prints data on the top 5 bins with the highest sensitivity:

```python
input_dir = r'path/to/your/repository/Input/2000-2500/'
filename = r'CO2_list.out'

# Read and process data
molecule, wavelength, intensity = read_hitran_data(input_dir, filename)
start_wavelength = 2200.0  # nm
end_wavelength = 2400.0    # nm
bin_size = 0.05            # nm

midpoints, binned_sensitivity, binned_line_counts = process_hitran_data(wavelength, intensity, start_wavelength, end_wavelength, bin_size)

# Print top sensitivity peaks
find_sensitivity_peaks(midpoints, binned_sensitivity)

# Plot the sensitivity
plot_sensitivity(molecule, midpoints, binned_sensitivity)
