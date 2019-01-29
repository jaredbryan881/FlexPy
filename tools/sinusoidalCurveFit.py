import matplotlib.pyplot as plt
from scipy import optimize
import numpy as np
import argparse
import h5py
import sys

sys.path.append('../../')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-pd', '--plotData', type=int, default=None,
                        help='Plot vertical normal stresses every n time steps.')
    parser.add_argument('-pw', '--plotWavelengths', type=bool, default=False,
                        help='Plot load wavelengths through time (True) or not (False)')
    parser.add_argument('-i', '--inFile', type=str,
                        default='/Users/jared/Desktop/URCO/Matlab_Codes/vert_norm_stress_50_60.txt',
                        help='Path to input vertical normal stresses, comma delimited csv')
    parser.add_argument('-s', '--saveData', type=bool, default=False,
                        help='Save data to hdf5 (True) or not (False)')
    args = parser.parse_args()

    vertNormStress = np.genfromtxt(args.inFile, delimiter=',')
    # number of data points
    N = vertNormStress[:].shape[1]
    xlin = np.linspace(0, 50, N)

    # fit sinusoid to data for each time step and store sinusoid params
    wavelength = []
    time = []
    params = np.empty((vertNormStress.shape[0], 4))
    zeros = np.empty((vertNormStress.shape[0], 2))

    for t in range(vertNormStress.shape[0]):
        # fit sinusoid
        params[t], paramsCov = optimize.curve_fit(testFunc, xlin, vertNormStress[t], p0=[-2000000, 0.1, 0, 0])

        # get two zero crossings
        # assumes single cycle of sinusoid
        zeros[t, 0] = (np.arcsin(-params[t, 3]/params[t, 0]) - params[t, 2]) / params[t, 1]
        midpoint = ((np.pi/2) - params[t, 2]) / params[t, 1]
        zeros[t, 1] = midpoint + (midpoint - zeros[t, 0])

        print('Timestep {}'.format(t))
        print('Wavelength: {} km'.format((2*np.pi/params[t, 1]) * 30))
        print('Amplitude: {}'.format(params[t, 0]))
        print('\n')

        wavelength.append((2*np.pi/params[t, 1]) * 30000)
        time.append(t)

        if args.plotData is not None and t % args.plotData == 0:
            # plot best fit sinusoid over data
            plt.plot([zeros[t, 0], zeros[t, 1]], [0, 0], 'ro')
            plt.plot(vertNormStress[t], '.', label='Data')
            plt.plot(xlin, testFunc(xlin, params[t, 0], params[t, 1], params[t, 2], params[t, 3]), label='Fitted Curve {}'.format(t))
            plt.title('Vertical Normal Stresses')
            plt.legend()
            plt.show()

    if args.saveData:
        hf = h5py.File('curveFitParams.h5', 'w')
        grp = hf.create_group('params')
        grp.create_dataset('params', data=params)
        grp.create_dataset('zeros', data=zeros)
        grp.create_dataset('wavelength', data=wavelength)
        grp.create_dataset('time', data=time)
        grp.create_dataset('stress', data=vertNormStress)
        hf.close()

    if args.plotWavelengths:
        plt.plot(time, wavelength)
        plt.title('Load Wavelength Over Time')
        plt.xlabel('Time Step')
        plt.ylabel('Wavelength (km)')
        plt.show()


def testFunc(x, a, b, c, d):
    """Return a sinusoid with a given amplitude, frequency, phase shift, and vertical shift.

    Params:
        @param x: domain of function
        @param a: amplitude
        @param b: frequency
        @param c: phase shift
        @param d: vertical shift
    Returns:
        @return: sinusoid
    """

    return a * np.sin(b * x + c) + d


if __name__ == '__main__':
    main()
