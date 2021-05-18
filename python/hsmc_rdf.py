import os, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import glob
import gzip
import hsmc_stat

# ------ Plot average radial distribution function ------

def plot(rdf_file):

    # Load rdf data
    rdf = np.loadtxt(rdf_file)

    # Plot
    plt.plot(rdf[:,0],rdf[:,1],'b')
    plt.ylabel('g(x)')
    plt.xlabel('x = r/sigma')
    plt.show()


# ------ Compute average radial distribution function ------

def average(file_names=None,data_dir=None, out_dir=None,lines_header=7):
    
    # Assign default input
    if file_names == None:
        file_names = 'rdf_*.dat.gz'
    if data_dir == None:
        data_dir = os.getcwd()
    if out_dir == None:
        out_dir = os.getcwd()

    # Get names of files in data directory
    file_names = sorted(glob.glob(os.path.join(data_dir,file_names)))
    n_files = len(file_names)
    if n_files == 0:
        sys.exit('hsmc_rdf.average: No rdf file was found')
    
    # Number of samples 
    n_samples = 0

    # Compute average
    for name in file_names:

        print(name)

        # Read file
        rdf_tmp = read_rdf_file(name,lines_header)
     
        # Update average
        if (n_samples == 0):
            rdf = np.zeros((rdf_tmp.shape[0],2))
        rdf[:,1] += np.sum(rdf_tmp[:,1:],axis=1)

        # Update sample counter
        n_samples += rdf_tmp.shape[1] - 1
        
                        
    # Normalize the rdf
    rdf[:,1] /= n_samples

    # Define the set of interparticle distances
    rdf[:,0] = rdf_tmp[:,0]
            
    # Output
    out_name = os.path.join(out_dir, 'rdf_average_'+str(n_samples)+'config.dat')
    np.savetxt(out_name, rdf, fmt='%.16e')


# ------ Perform blocking analysis via Flyvberg method ------

def blocking(target_pos, file_names=None, data_dir=None, out_dir=None,lines_header=7,plt=False):
    
    # Assign default input
    if file_names == None:
        file_names = 'rdf_*.dat.gz'
    if data_dir == None:
        data_dir = os.getcwd()
    if out_dir == None:
        out_dir = os.getcwd()

    # Get names of files in data directory
    file_names = sorted(glob.glob(os.path.join(data_dir,file_names)))
    n_files = len(file_names)
    if n_files == 0:
        sys.exit('hsmc_rdf.average: No rdf file was found')

        
    # Number of samples 
    n_samples = 0

    # Index for target positions 
    n_target = len(target_pos)
    target_pos_idx = np.full((n_target,1),-1,dtype=int)

    # Get data from output files
    for name in file_names:

        print(name)

        # Read file
        rdf_tmp = read_rdf_file(name,lines_header)
     
        # Find index for target positions
        if (n_samples == 0):
            for ii in range(n_target):
                for jj in range(len(rdf_tmp[:,0])):
                    if rdf_tmp[jj,0] >= target_pos[ii]:
                        target_pos_idx[ii] = jj
                        break
                if target_pos_idx[ii] == -1:
                    sys.exit('Target position %.8f could not be found in the rdf files' % target_pos[ii])
                
                target_pos[ii] = rdf_tmp[target_pos_idx[ii],0]

        # Update rdf data
        if (n_samples == 0):
            rdf = np.squeeze(rdf_tmp[target_pos_idx,1:])
        else:
            rdf = np.concatenate((rdf,np.squeeze(rdf_tmp[target_pos_idx,1:])),axis=1)

        # Update sample counter
        n_samples += rdf_tmp.shape[1] - 1

    # Initialize output file for blocking analysis
    ff = open('rdf_blocking_analysis.dat','w')
    ff.write('Blocking analysis of the radial distribution function data in ' + data_dir + '\n')
    ff.write('Total number of samples: %d\n' % n_samples)
    ff.write('Target positions:\n')
    for ii in range(n_target):
        ff.write('%.8f\n' % target_pos[ii])

    # Blocking analysis for the target distances
    for ii in range(n_target):
        
        # Print on screen
        print('###############################################')
        print('Blocking analysis for target distance %.8f' % target_pos[ii] )
        sigma = hsmc_stat.blocking_std(rdf[ii,:],plt_flag=plt)
    
        # Print to file
        ff.write('###############################################\n')
        ff.write('Blocking analysis for target distance %.8f\n' % target_pos[ii] )
        ff.write('Transformation, variance, variance uncertainty\n')
        nn_trans = len(sigma[:,0])
        for jj in range(nn_trans):
            ff.write('%d %.8f %.8f\n' % (jj, sigma[jj,0], sigma[jj,1]))
       
    ff.close()
        

    
# ------ Read  HSMC output for the radial distribution function ------

def read_rdf_file(file_name,lines_header=7):

    # Open file, read all lines and close
    with gzip.open(file_name) as ff:
        lines_data = ff.readlines()
        lines_file = len(lines_data)
        
    # Extract number of bins, volume and number of particles
    [sweep, n_bins, vol, n_part]  = [float(kk) for kk in lines_data[3].split()]
    n_bins = int(n_bins)
        
    # Number of lines per sample and per block
    lines_sample = n_bins + lines_header

    # Number of samples
    n_samples = lines_file/lines_sample

    # Initialize rdf
    rdf = np.zeros((n_bins,n_samples+1))

    # Read file
    lines_read = 0        
    n_samples_read = 0
    while lines_read <= lines_file - lines_sample:
            
        # Read one sample and update output
        for ii in range(n_bins):
            rdf_tmp = np.array([float(kk) for kk in 
                                lines_data[ii+lines_read+lines_header].split()])
            if n_samples_read == 0:
                rdf[ii,n_samples_read] = rdf_tmp[0]
            rdf[ii,n_samples_read+1] = rdf_tmp[1] 
                
        # Update sample counter
        n_samples_read += 1
        if n_samples_read > n_samples:
            sys.exit('hsmc_rdf.read_rdf_file: Bad file structure, more samples then expected')
    
        # Update the number of read lines
        lines_read += lines_sample

    # Check that the correct number of samples was read
    if n_samples_read != n_samples:
        sys.exit('hsmc_rdf.read_rdf_file: Bad file structure, less samples then expected')

    # Output
    return rdf

