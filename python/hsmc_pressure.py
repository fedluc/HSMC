import os, sys
import numpy as np
import math
from scipy import optimize as opt
from scipy import special as spec
import glob

# ------ Compute pressure from virial route ------
def pressure_virial(data_dir,samples_block=1):

    # Lines for the header of each sample (assumed fixed)
    lines_header = 7

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,'press_virial.dat'))
    out_dir = data_dir
    n_files = len(file_names)
    if n_files == 0:
        sys.exit('hsmc_pressure.pressv: No pressure file was found')
    
    # Samples for the radial distribution function at contact
    [coeff, vol, n_part] = read_pressure_output(file_names,samples_block,lines_header,"virial")
    rdf_contact = coeff[:,0] + coeff[:,1]

    # Average and standard deviation of the pressure
    dens = n_part/vol
    n_samples = len(coeff)
    press = 1 + (2.*math.pi/3.) * dens * rdf_contact
    press_ave = np.average(press)
    press_var = np.var(press)/np.sqrt(n_samples)
    print("Pressure from virial calculations")
    print("Samples, average (reduced units), variance (reduced units), average (HS units), variance (HS units))")
    print("%d %.8f %.8f %.8f %.8f" % (n_samples, press_ave, press_var, press_ave*dens, press_var*dens))
    

# ------ Compute pressure from virial route ------
def pressure_thermo(data_dir=os.getcwd(),samples_block=1,npt=False):

    # Lines for the header of each sample (assumed fixed)
    lines_header = 7

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,'press_thermo.dat'))
    out_dir = data_dir
    if len(file_names) == 0:
        sys.exit('hsmc_pressure.pressv: No pressure file was found')
    
    # Samples for the excess pressure
    [coeff, vol, n_part] = read_pressure_output(file_names,samples_block,lines_header,"thermo")
    p_ex = coeff[:,1]/n_part

    # Samples for the density
    if npt:
        lines_header = 3
        file_names = glob.glob(os.path.join(data_dir,'density.dat'))
        if len(file_names) == 0:
            sys.exit('hsmc_pressure.pressv: No density file was found')
        dens = read_density_output(file_names,samples_block,lines_header)
    else:
        dens = n_part/vol

    # Average and standard deviation of the pressure
    if npt:
        n_press_samples = len(coeff)
        n_dens_samples = len(dens)
        press = dens * (1 + p_ex)
        press_ave = np.average(press)
        press_var = np.var(press)/np.sqrt(n_press_samples)
        dens_ave = np.average(dens)
        dens_var = np.var(dens)/np.sqrt(n_dens_samples)
        print("Pressure from thermodynamic calculations")
        print("Samples, average (HS units), variance (HS units))")
        print("%d %.8f %.8f" % (n_press_samples, press_ave, press_var))
        print("Density from thermodynamic calculations")
        print("Samples, average (HS units), variance (HS units))")
        print("%d %.8f %.8f" % (n_dens_samples, dens_ave, dens_var))
    else:
        n_samples = len(coeff)
        press = 1 + p_ex
        press_ave = np.average(press)
        press_var = np.var(press)/np.sqrt(n_samples)
        print("Pressure from thermodynamic calculations")
        print("Samples, average (reduced units), variance (reduced units), average (HS units), variance (HS units))")
        print("%d %.8f %.8f %.8f %.8f" % (n_samples, press_ave, press_var, press_ave*dens, press_var*dens))

    

# ------ Read  HSMC output for the pressure -----
def read_pressure_output(file_names,samples_block,lines_header,press_flag):

    # Global variables declaration
    global lines_press

    # Read files listed in file_names
    init_hist_flag = True
    init_file_flag = True
    for name in file_names:
            
        # Open file, read all lines and close
        ff = open(name)
        lines_press = ff.readlines()
        ff.close()
        lines_file = len(lines_press)
                
        # Extract number of bins, volume and number of particles
        [n_bins, vol, n_part]  = [float(kk) for kk in lines_press[3].split()]
        n_bins = int(n_bins)

        # Initialize histograms and pressure
        hist = np.zeros((n_bins,2))

        # Number of lines per sample and per block
        lines_sample = n_bins + lines_header
        lines_block = samples_block*lines_sample

        # Number of usable samples in file
        n_samples_file = lines_file//lines_block
        coeff_file = np.zeros((n_samples_file,2))
        
        # Read file
        n_samples = 0
        n_blocks = 0;
        lines_read = 0;        
        while lines_read <= lines_file - lines_block:

            # Collect samples for average
            while n_samples < samples_block:
            
                # Read one sample
                histTmp = read_pressure_sample(lines_read+lines_header, n_bins)

                # Create vector with abscissa for fit 
                if init_hist_flag:
                    hist[:,0] = histTmp[:,0]
                    init_hist_flag = False

                # Update running average
                hist[:,1] += histTmp[:,1]/samples_block

                # Update number of samples
                n_samples += 1

                # Update the number of read lines
                lines_read += lines_sample

                
            # Prepare data for fitting
            if (press_flag == "thermo"):

                # Remove indeterminate forms
                del_idx = np.argwhere(hist[:,1]==0)
                hist = np.delete(hist,del_idx,axis=0) 
                
                # Create function to fit                
                hist[:,1] = -np.log(hist[:,1])

            # Fit data
            coeff_file[n_blocks,:] = data_fit(hist)

            # Update number of blocks 
            n_blocks+=1
                     
            # Reset histogram
            hist = np.zeros((n_bins,2))
            init_hist_flag = True

            # Reset number of samples
            n_samples = 0

            # Remove nan values 
            del_idx = np.argwhere(np.isnan(coeff_file))
            coeff_file = np.delete(coeff_file,del_idx,axis=0)
        
        if init_file_flag:
            coeff = coeff_file
            init_file_flag = False;
        else:
            coeff = np.append(coeff,coeff_tmp,axis=0)

    # Output
    return [coeff, vol, n_part]



def read_pressure_sample(line_start, n_bins):

    hist = np.zeros((n_bins,2))

    for jj in range(n_bins):
        histTmp = np.array([float(kk) for kk in lines_press[jj+line_start].split()])
        hist[jj,0] = histTmp[0]
        hist[jj,1] = histTmp[1] 
    
    return hist

def linear_fit(x, b0, b1):
    return b0 + b1*x

def data_fit(hist):

    if (np.shape(hist)[0] >= 2):
        [fit_par, fit_par_cov] = opt.curve_fit(linear_fit,
                                               hist[:,0], hist[:,1],
                                               p0=[10.0,1.0])
    else:

        fit_par = [np.nan, np.nan]

    
    return fit_par



# ------ Read  HSMC output for the density -----

def read_density_output(file_names,samples_block,lines_header):

    # Read files listed in file_names
    init_file_flag = True
    for name in file_names:
            
        # Open file, read all lines and close
        ff = open(name)
        lines_dens = ff.readlines()
        ff.close()
        lines_file = len(lines_dens)
                
        # Number of usable samples in file
        lines_sample = 1
        lines_block = samples_block*lines_sample
        n_samples_file = (lines_file-1)//lines_block
        dens_file = np.zeros((n_samples_file,1))
        
        # Read file
        n_samples = 0
        n_blocks = 0;
        lines_read = 0;        
        while lines_read <= lines_file - lines_block:

            # Collect samples for average
            while n_samples < samples_block:
            
                # Update running average
                dens_file[n_blocks] += float(lines_dens[lines_read+lines_header])/samples_block

                # Update number of samples
                n_samples += 1

                # Update the number of read lines
                lines_read += lines_sample

            # Update number of blocks 
            n_blocks+=1
                       
            # Reset number of samples
            n_samples = 0
        
        if init_file_flag:
            dens = dens_file
            init_file_flag = False
        else:
            dens = np.append(dens,dens_file,axis=0)


    # Output
    return dens

