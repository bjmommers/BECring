#!/usr/bin/python3

# Python script for calculating Fisher information from the output of three
#   BECring simulations. Designed to be used as part of fisherrun.sh

import sys
import pandas as pd
import numpy as np

print("Calculating Fisher information...")

# Use the first argument as the root directory
root_dir = sys.argv[1]

# Pull the value of delta from file in the root directory
delta_file = open(root_dir + "/fisher-delta.dat",'r')
delta_omega = float(delta_file.read())

# First calculate the classical Fisher info (CFI), which only requires the density
# densitytime-labelled.dat has time in the first column, and x-values in the first row
densitytime_minus = pd.read_csv(root_dir + "/fisher-output-minus/" + "densitytime.dat",sep='\t',header=None)
densitytime_omega = pd.read_csv(root_dir + "/fisher-output-omega/" + "densitytime.dat",sep='\t',header=None)
densitytime_plus = pd.read_csv(root_dir + "/fisher-output-plus/" + "densitytime.dat",sep='\t',header=None)
densitytime_labelled_omega = pd.read_csv(root_dir + "/fisher-output-omega/" + "densitytime-labelled.dat",sep='\t',header=0)

# Isolate time values into a dataframe column
time = densitytime_labelled_omega["time"]

# Test that the dataframes are all the same size (they should be!)
if (densitytime_minus.size != densitytime_plus.size):
    print("Warning: densitytime dataframes are not the same size!")
    print("This could be due to incorrect simulation input,")
    print("or corrupted data. Results may not be viable!")

# Extract the spatial step size from the labelled dataframe
dx = float(densitytime_labelled_omega.columns[2]) - float(densitytime_labelled_omega.columns[1])

# To calculate CFI, use the central finite difference to approximate the derivative
cfi_finite_diff = densitytime_plus.sub(densitytime_minus)/(2.0*delta_omega)
# Construct the integrand
cfi_integrand = ((cfi_finite_diff ** 2.0)/densitytime_omega)*dx
# Integrate by summing over row contents (i.e. over space)
cfi = cfi_integrand.sum(axis=1)
# Combine CFI and time into a neat data file for plotting with gnuplot
cfi_time = pd.concat([time,cfi],axis=1)
cfi_time.to_csv(root_dir + "/classical_fisher_info.dat",sep='\t')


# Check we want to calculate QFI, otherwise exit gracefully
if int(sys.argv[2]) != 1:
    print("Skipping QFI calculation...")
    sys.exit()



# Now calculate the quantum Fisher info

# Start by declaring a dataframe
qfi_time = pd.DataFrame({'time':[], 'QFI':[]})

# Iterate through files for the complex wavefunction
for t in range(0,len(cfi_time.index)):

    #load files
    frame_data_minus = pd.read_csv(root_dir + "/fisher-output-minus/frames/frame{:05d}.dat".format(t),sep='\t')
    frame_data_omega = pd.read_csv(root_dir + "/fisher-output-omega/frames/frame{:05d}.dat".format(t),sep='\t')
    frame_data_plus = pd.read_csv(root_dir + "/fisher-output-plus/frames/frame{:05d}.dat".format(t),sep='\t')

    #create dataframes that contain the complex wavefunction
    wavefunction_minus = frame_data_minus["Re(psi1)"] + frame_data_minus["Im(psi1)"] * 1j
    wavefunction_omega = frame_data_omega["Re(psi1)"] + frame_data_omega["Im(psi1)"] * 1j
    wavefunction_plus = frame_data_plus["Re(psi1)"] + frame_data_plus["Im(psi1)"] * 1j

    #calculate central finite difference of wavefunction 
    wf_finite_diff = (wavefunction_plus - wavefunction_minus)/(2.0*delta_omega)

    #compute the first term of the QFI
    first_term = (np.conjugate(wf_finite_diff)*wf_finite_diff*dx).sum(axis=0)

    #compute the second term of the QFI
    second_term = pow(abs((np.conjugate(wavefunction_omega)*wf_finite_diff*dx).sum(axis=0)),2.0)

    #compute the final QFI for this timestep
    qfi = 4.0*(first_term - second_term)

    #append time and QFI to final dataframe
    #abs used to cast back to float for plotting, test data showed no imaginary component
    qfi_time = qfi_time.append({'time':abs(cfi_time["time"][t]),'QFI':abs(qfi)},ignore_index=True)

#Once the QFI dataframe has been filled, save it to file
qfi_time.to_csv(root_dir + "/quantum_fisher_info.dat",sep='\t')
