import os, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import glob
import hsmc_stat

# ---------------------------------------------------------
# Note: these functions analyze the cavity output produced
# when hsmc is run with the cavity option
# ---------------------------------------------------------

# ------ Radial distribution function ------

def rdf(data_dir,file_id='cavity_distance.dat',file_comments='#',
        dr=0.001,rmax=1.4,rmin=0.001,nPart=125, normalize=True, plot=True):

    # Grid for the rdf extraction
    rr = np.arange(rmin-dr/2,rmax+dr,dr)
    bin_vol = (4.*math.pi/3) * (rr[1:]**3 - rr[:-1]**3)
    nn_rr = len(rr) - 1

    # Read data 
    data = read_hsmc_output(data_dir, file_id, file_comments)
    
    # Histograms of positions
    hist = np.histogram(data, bins=rr)
    
    # Radial distribution function
    rdf = np.zeros((nn_rr,2))
    rdf[:,0] = (rr[1:] + rr[:-1])/2.
    rdf[:,1] = hist[0]/(bin_vol*4./nPart)
    if normalize:
        rdf[:,1] /= bin_vol*4./nPart

    # Plot
    if plot:
        rdfPlot = rdf[rdf[:,1]>0,:];
        plt.plot(rdfPlot[:,0],np.log(rdfPlot[:,1]),'b')
        plt.xlabel('Distance [sigma]')
        plt.ylabel('ln[g(r)]')
        plt.show()
    
    # Output
    return rdf


# ------ Interaction potential ------

def psi(data_dir,file_id='cavity_psi.dat',file_comments='#',
        plot=True):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments)
    
    # Plot
    if plot:
        plt.plot(data[:,0],data[:,1],'b')
        plt.ylabel('Cavity interaction potential [k_BT]')
        plt.xlabel('Distance [sigma]')
        plt.show()
    
    # Output
    return data

# ------ Logarithm of the cavity function ------

def lny(data_dir,file_id_dist='cavity_distance.dat',
        file_id_psi='cavity_psi.dat',file_comments='#',
        nPart=125, normalize=True, plot=True):

    # Interaction potential between the cavities
    pot = psi(data_dir,file_id_psi,file_comments,False)
    pot = pot[1:,:]
    rmax = pot[pot.shape[0]-1,0]
    rmin = pot[0,0]
    dr = pot[1,0] - pot[0,0]

    # Radial distribution function
    gg = rdf(data_dir,file_id_dist,file_comments,
             dr,rmax,rmin,nPart,True,False)

    # Remove indeterminate forms
    mask = gg[:,1]>0
    gg = gg[mask,:]
    pot = pot[mask,:]
    
    # Logarithm of the cavity
    lny = pot
    lny[:,1] += np.log(gg[:,1])

    # Plot
    if plot:
        plt.plot(lny[:,0],lny[:,1],'b')
        plt.xlabel('Distance [sigma]')
        plt.ylabel('ln[y(r)]')
        plt.show()
    
    # Output
    return lny


# ------ Read cavity output of hsmc ------

def read_hsmc_output(data_dir,file_id,file_comments):

    # Get names of files in data directory
    file_names = glob.glob(os.path.join(data_dir,file_id))
    out_dir = data_dir
    if len(file_names) == 0:
        sys.exit('hsmc_cavity.read_hsmc_output: No data file was found')

    # Read files listed in file_names
    init_file_flag = True
    for name in file_names:

        # Open file, read all lines and close
        data_file = np.loadtxt(name, comments=file_comments)

        # Update output
        if init_file_flag:
            data = data_file
            init_file_flag = False
        else:
            data = np.append(data,data_file,axis=0)

    # Output
    return data











