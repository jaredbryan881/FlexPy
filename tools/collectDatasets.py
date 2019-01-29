import numpy as np
import h5py

def collectDatasets(basePath, outFile, datasets, numFiles, dims):
    """Combine text files containing stresses, deflections, etc. and output hdf5.

    Args:
        :param basePath: list of str
            list of base paths to files to collect
        :param outFile: str
            path to output hdf5
        :param datasets: list of str
            name of datasets to write
        :param numFiles: int
            number of files to iterate over
        :param dims: tuple
            dimension of output data
    """
    hf = h5py.File(outFile, 'w')
    for i, file in enumerate(basePath):
        data = np.empty(dims)
        for f in range(numFiles):
            data[:, f, :] = np.genfromtxt(file.format(f+1), delimiter=',')
        hf.create_dataset(datasets[i], data=data)
    hf.close()

collectDatasets(['/Users/jared/Desktop/URCO/Matlab_codes/normStress_{}_0_20.txt', '/Users/jared/Desktop/URCO/Matlab_codes/shearStress_{}_0_20.txt'],
                'stressStateConf_20.h5',
                ['sigmayy', 'sigmaxy'],
                17,
                (201, 17, 101))
