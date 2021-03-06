import os, sys
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.fftpack as fft
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


# ------ Compute direct correlation function ------

def dcf(rdf_file,rho,nPart=13500,out_dir=None,plot=False):

    # Assign default input
    if out_dir == None:
        out_dir = os.getcwd()

    # Load rdf data
    rdf = np.loadtxt(rdf_file)

    # Finite size correction
    rdf[:,1] = finite_size_correction(rdf[:,1],rho,nPart)

    # Grid for Fourier transform
    rr = rdf[:,0]
    dr = rr[1] - rr[0]
    rr = np.arange(dr/2.0,rr[-1]+dr/2,dr)
    nn = len(rr)
    dq = math.pi/(nn*dr)
    qq = np.arange(nn)*dq + dq

    # Total correlation function
    hh = np.zeros(nn)
    mask = rr>=1
    hh[mask] = rdf[:,1]
    hh -= 1.0

    # Fourier transform of the total correlation function
    hhft = 4*math.pi*dr/qq * fft.dst(rr*hh,type=2)/2.0

    # Static structure factor (to check)
    ssf = 1.0 + rho*hhft

    # Fourier transform of the indirect correlation function
    gammaft = rho*hhft**2.0/(1.0 + rho*hhft)

    # Indirect correlation function
    gamma = dq/(2*math.pi**2*rr) * fft.dst(qq*gammaft,type=3)/2.0

    # Direct correlation function
    cc = np.zeros((nn,2))
    cc[:,0] = rr
    cc[:,1] = hh - gamma

    # Plot
    if plot:
        plt.plot(cc[:,0],cc[:,1],'b')
        plt.ylabel('c(x)')
        plt.xlabel('x = r/sigma')
        plt.show()

    # Output
    out_name = os.path.join(out_dir, 'dcf.dat')
    np.savetxt(out_name, cc, fmt='%.16e')

    return rdf,cc,gamma


# ------ Compute bridge function ------

def bf(rdf_file,rho,nPart=13500,lny_file=None,
       out_dir=None,plot=False):

    # Assign default input
    sr_flag = True
    if lny_file == None:
        sr_flag = False
    if out_dir == None:
        out_dir = os.getcwd()

    # Correlation functions
    rdf, cc, gamma = dcf(rdf_file,rho,nPart,out_dir,False)

    # Grid
    rr = cc[:,0]
    nn = len(rr)

    # Cavity function
    if sr_flag:
        lny = np.loadtxt(lny_file)
    else:
        lny = np.zeros((nn,2))
        lny[:,0] = rr
        lny[:,1] = float("nan")
        mask = rr >= 1.0
        lny[mask,1] = np.log(rdf[:,1])

    # Bridge function
    bb = np.zeros((len(rr),2))
    bb[:,0] = rr
    bb[:,1] = lny[:,1] - gamma

    if plot:
        mask = rr<3.0
        plt.plot(bb[mask,0],-bb[mask,1],'b')
        plt.ylabel('-B(x)')
        plt.xlabel('x = r/sigma')
        plt.show()

        mask = (rr > 1.2) &  (rr<5.0)
        plt.plot(bb[mask,0],-bb[mask,1],'b')
        plt.ylabel('-B(x)')
        plt.xlabel('x = r/sigma')
        plt.show()

    # Output
    out_name = os.path.join(out_dir, 'bf.dat')
    np.savetxt(out_name, bb, fmt='%.16e')


def finite_size_correction(rdf,rho,nPart):
    
    # Inverse isothermal compressibility
    mu_T = iic_ew(rho*math.pi/6.0)/rho

    # Finite size corrected rdf
    rdf *= 1.0 + 1.0/(nPart*mu_T)
    
    return rdf


def iic_ew(eta):

    # Inverse isothermal compressibility from Erpenbeck-Wood EOS
    
    # Compressibility factor
    zz1 = 0.890851
    zz2 = 0.8924486
    zz3 = 0.3430298
    zz4 = -2.277287
    zz5 = 1.3262418
    zz = eta * ( (4.0 +zz1*eta +zz2*eta**2.0 + zz3*eta**3.0)
                /(1.0 + zz4*eta + zz5*eta**2) )
    
    zz_der_1 = ((eta * (zz1 + 2.0*zz2*eta + 3.0*zz3*eta**2.0))/
              (1 + zz4*eta + zz5*eta**2.0))  
    zz_der_2 = -((eta * (zz4 + 2.0*zz5*eta) * 
                  (4.0 + zz1*eta + zz2*eta**2.0 + zz3*eta**3.0))
                 /(1 + zz4*eta + zz5*eta**2)**2.0) 
    zz_der_3 = ((4.0 + zz1*eta + zz2*eta**2.0 + zz3*eta**3.0)/
                (1.0 + zz4*eta + zz5*eta**2.0))
    zz_der = zz_der_1 + zz_der_2 + zz_der_3

    return 1.0 + eta*zz_der + zz

     
