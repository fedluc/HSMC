import os, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import glob
import scipy.interpolate as scint
import hsmc_stat

# ---------------------------------------------------------
# Note: these functions analyze the cavity output produced
# when hsmc is run with the cavity option
# ---------------------------------------------------------

# ------ Logarithm of the cavity function (physical system) ------

def lny(intervals,data_dir_lr,data_dir_sr,out_dir=None,
        file_rdf_lr=None, file_id_dist=None, file_id_psi=None,
        file_comments='#',
        nPart=1000, plot=False):

    # Assign default input
    if out_dir == None:
        out_dir = os.getcwd()
    if file_rdf_lr == None:
        file_rdf_lr = 'rdf_average*config.dat'
    if file_id_dist == None:
        file_id_dist = 'cavity_distance.dat'
    if file_id_psi == None:
        file_id_psi = 'cavity_psi.dat'
    
    # Check that the input is consistent
    n_intervals = intervals.shape[0]
    n_data_dir_sr = len(data_dir_sr)
    if n_intervals != n_data_dir_sr:
        sys.exit('The number of intervals (%d) is not equal to the number of short range directories (%d)' % (n_intervals,n_data_dir_sr))

    # Order the intervals 
    order = intervals[:,0].argsort()[::-1]
    intervals = intervals[order]
    data_dir_sr = np.array(data_dir_sr)[order]

    # Load long-range rdf
    file_name = sorted(glob.glob(os.path.join(data_dir_lr,file_rdf_lr)))
    n_files = len(file_name)
    if n_files > 1:
        sys.exit('Error: more than one file matches the name given for the file with the long range radial distribution function')
    elif n_files == 0:
        sys.exit('Error: file with long range radial distribution function not found')
    file_name = file_name[0]
    rdf_lr = np.loadtxt(file_name)
    
    # Grid for the matching
    rr = rdf_lr[:,0]
    dr = rr[1] - rr[0]
    rr = np.arange(dr/2.0,rr[-1]+dr/2,dr)
    nr = len(rr)

    # Initialize cavity
    lny = np.zeros(nr)
    lny_tmp = np.zeros(nr)

    # Long-range cavity
    mask = rr>=1
    lny[mask] = np.log(rdf_lr[:,1])

    # Short range-cavity
    for ii in range(n_intervals):

        # Confining positions 
        mask_confine = (rr>=intervals[ii,0]) & (rr<=intervals[ii,1]) 
    
        # Cavity from cavity simulations
        lny_tmp[mask_confine] = lny_sim(data_dir_sr[ii],rr[mask_confine],file_id_dist,file_id_psi,file_comments,nPart,False)[:,1]

        # Matching positions
        r_max_match = intervals[ii,1]
        if ii==0:
            r_min_match = 1.0
        else:
            r_min_match = intervals[ii-1,0]
        mask_match = (rr>=r_min_match) & (rr<=r_max_match)

        # Matching constant
        cc = lny_tmp[mask_match] - lny[mask_match]

        # Update cavity 
        lny[mask_confine] = lny_tmp[mask_confine] - np.mean(cc)
        
        # Plot
        plt.plot(rr[mask_match],cc,'bo')
        plt.ylabel('Matching constant')
        plt.xlabel('x = r/sigma')
        plt.show()

    mask = rr < 2.0
    plt.plot(rr[mask],lny[mask],'b')
    plt.ylabel('ln[y(r)]')
    plt.xlabel('x = r/sigma')
    plt.show()

        
    

# ------ Logarithm of the cavity function (simulated system) ------

def lny_sim(data_dir,grid,file_id_dist='cavity_distance.dat',
            file_id_psi='cavity_psi.dat',file_comments='#',
            nPart=1000, plot=False):

    # Interaction potential between the cavities
    pot = psi(data_dir,grid,file_id_psi,file_comments,False)

    # Radial distribution function
    gg = rdf(data_dir,grid,file_id_dist,file_comments,
             nPart,True,False)

    # Remove indeterminate forms
    mask = gg[:,1]<=0
    gg[mask,1] = float("nan")

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

# ------ Interaction potential ------

def psi(data_dir,grid,file_id='cavity_psi.dat',file_comments='#',
        plot=False):

    # Read data
    data = read_hsmc_output(data_dir, file_id, file_comments)

    # Interpolate to the desired grid
    data_interp = scint.interp1d(data[:,0], data[:,1], 
                                 kind='cubic', fill_value="extrapolate")

    # Plot
    if plot:
        plt.plot(grid,data_interp(grid),'b')
        plt.ylabel('Cavity interaction potential [k_BT]')
        plt.xlabel('Distance [sigma]')
        plt.show()
    
    # Output
    out = np.zeros((len(grid),2))
    out[:,0] = grid
    out[:,1] = data_interp(grid)
    return out



# ------ Radial distribution function ------

def rdf(data_dir,grid,file_id='cavity_distance.dat',file_comments='#',
        nPart=1000, normalize=True, plot=False):

    # Grid for the rdf extraction
    dr = grid[1] - grid[0]
    grid_shift = np.arange(grid[0]-dr/2,grid[-1]+dr,dr)
    bin_vol = (4.*math.pi/3) * (grid_shift[1:]**3 - grid_shift[:-1]**3)
    nn = len(grid)

    # Read data 
    data = read_hsmc_output(data_dir, file_id, file_comments)
    
    # Histograms of positions
    hist = np.histogram(data, bins=grid_shift)
    
    # Radial distribution function
    rdf = np.zeros((nn,2))
    rdf[:,0] = grid
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











