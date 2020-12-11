import numpy as np

# Get names of files in data directory
file_id = 'samples*'
file_names = glob.glob(os.path.join(data_dir,file_id))
out_dir = data_dir
if len(file_names) == 0:
    sys.exit('hsmc_pressure.pressv: No pressure file was found')
