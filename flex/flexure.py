import matplotlib.pyplot as plt
from scipy import optimize
from scipy import signal
import numpy as np
import argparse
import h5py
import sys

sys.path.append('../../')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-il', '--inLoad', type=str, required=False, default=None,
                        help='Path to txt or csv containing vertical normal stresses')
    parser.add_argument('-ip', '--inParams', type=str, required=False, default=None,
                        help='Path to hdf5 containing curve fit parameters for vertical norm stress.')
    args = parser.parse_args()
    args.inLoad = '/Users/jared/Desktop/URCO/Matlab_Codes/vertNormStress.txt'

    # create instance of flexure class
    flex = Flexure()
    # get vertical normal stresses
    if args.inLoad is not None:
        # read vertical normal stresses
        lithLoad = flex.readStress(args.inLoad)
    elif args.inParams is not None:
        params, zeros, wavelength, time, lithLoad = flex.readParams(args.inParams)
    else:
        raise ValueError("Input file is None. You must specify an input file of stresses or curve fit parameters.")

    if lithLoad.ndim == 1:
        lithLoad = np.reshape(lithLoad[1:], (1, lithLoad.shape[0]-1))
    stepmax = lithLoad.shape[0]

    # compute flexure calculate for each time step
    for tStep in range(stepmax):
        greenDef = flex.deflectionGreensFunction(flex.x)
        deflection = flex.convGreenLoad(greenDef, lithLoad[tStep])
        greenStress = flex.bendingStressGreensFunction(flex.x)

        sigxx = np.empty((flex.ynum, 2050))
        for d in range(flex.ynum):
            sigxx[d] = flex.convGreenLoad(greenStress[d], lithLoad[tStep])

        plt.subplot(2, 1, 1)
        plt.plot(deflection)
        plt.title('Deflection (m)')
        plt.subplot(2, 1, 2)
        plt.imshow(sigxx, cmap='seismic')
        plt.title('Bending Stress')
        #plt.colorbar()
        plt.clim(-5e+6, 5e+6)
        plt.show()


class Flexure:
    def __init__(self):
        """Initialize class for computing deflection and stresses of a sinusoidal load."""
        # simulation parameters
        self.xstp = 30000  # m -- x spacing
        self.ystp = 20000  # m -- y spacing
        self.xnum = 50  # number of nodes in x
        self.lithH = 100000  # m -- lithospheric thickness
        self.rhom = 3300  # kg/m**3 -- mantle density
        self.rhoin = 1000  # kg/m**3 -- infill density
        self.g = 9.81  # m/s**2 -- gravitational acceleration

        # material properties
        self.E = 6.5e+10  # Young's modulus (Pa)
        self.v = 0.25  # Poisson's ratio (dimensionless)
        self.I = (self.xstp * self.xnum) * (self.lithH ** 3 / 3)  # m**4 -- second moment of area
        self.D = (self.E * self.lithH ** 3) / (12 * (1 - self.v ** 2))  # Nm -- flexural rigidity
        self.alpha = ((4 * self.D) / ((self.rhom - self.rhoin) * self.g))**0.25  # flexural parameter
        self.lamb = 1 / self.alpha  # effective wavelength of the load

        # calculation parameters
        self.lowBound = 0  # m -- lower bound to where we calculate flexure and stress
        self.upBound = 1000000  # m -- upper bound to where we calculate flexure and stress
        self.spacing = 1000  # m -- spacing between sampled points
        self.numPoints = int((self.upBound - self.lowBound) / self.spacing)  # number of points in x
        self.x = np.linspace(self.lowBound, self.upBound, self.numPoints+1)  # range of points in x
        self.ynum = 201  # number of points in y

    def convGreenLoad(self, green, load):
        """Compute the deflection of an infinite plate under a distributed load.

        Params:
            :param green: np.array
                numpy array containing computed Green's function for a line load on an infinite plate
            :param load: np.array
                numpy array containing vertical normal stresses at the base of the lithosphere

        Returns:
            :return conv: np.array
                numpy array containing the convolution of a Green's function and a discretized load
        """
        w = np.zeros(2 * green.shape[0] - 1)
        nw = green.shape[0]
        w[nw:2*nw-1] = green[1:nw]
        w[:nw] = green[::-1]

        # convolve greens function with applied load
        conv = np.convolve(w, load)
        # make this a deflection per unit width
        conv *= self.xstp

        return conv

    def simplySupported(self, params, zeros, wavelength, load):
        """Compute deflection, bending stress, and membrane stress for a simply supported sinusoidal load.
        Args:
            :param params: np.array
                numpy array containing curve fit parameters
            :param zeros: np.array
                numpy array containing coordinates of zero crossings
            :param wavelength: list
                list containing wavelength of sinusoid at each time step
            :param load: np.array
                numpy array containing vertical normal stresses at the base of the lithosphere

        Returns:
            :return ss_deflection:
            :return ss_sigmaBend:
            :return ss_sigmaMembrane:
        """
        # get number of nodes in the simply supported beam
        ss_nodeNum = int((zeros[1] - zeros[0]) / self.xstp)
        # get length of beam
        ss_len = zeros[1] - zeros[0]

        # maximum stress on simply supported beam
        x = np.linspace(zeros[0], zeros[1], 1001)
        sin = -params[0] * np.sin(params[1] * x + params[2]) + params[3]

        zero0_loc = zeros[0] / self.xstp
        zero1_loc = zeros[1] / self.xstp
        s_mag = -np.sign(sin[int((zero1_loc - zero1_loc) / 2)]) * np.max(abs(sin[int(zero0_loc):int(zero1_loc)]))

        # compute deflection along the simply supported beam
        ss_deflection = -(s_mag * ss_len ** 4 * np.sin((np.pi * x) / ss_len)) / (np.pi ** 4 * self.E * self.I)

        plt.plot(ss_deflection)
        plt.show()

        # compute shear moment (V) and bending moment (M)
        M = ((s_mag * ss_len ** 2) / np.pi ** 2) * np.sin((np.pi * x) / ss_len)
        V = ((s_mag * ss_len) / np.pi) * np.cos((np.pi * x) / ss_len)

        # compute bending stresses in the simply supported beam
        neutralAxis = self.lithH / 2
        y = np.linspace(0, self.lithH, 1001)  # distance from neutral axis
        y = neutralAxis - y
        ss_sigmaBend = np.empty((y.shape[0], M.shape[0]))

        for i in range(ss_sigmaBend.shape[0]):
            ss_sigmaBend[i] = (M * y[i]) / self.I

        # compute membrane stresses in the simply supported beam
        ss_sigmaMembrane = 0

        return ss_deflection, ss_sigmaBend, ss_sigmaMembrane

    def cantilevered(self, params, zeros, wavelength, load):
        """Compute deflection, bending stress, and membrane stress for a cantilevered sinusoidal load.
        Args:
            :param params: np.array
                numpy array containing curve fit parameters
            :param zeros: np.array
                numpy array containing coordinates of zero crossings
            :param wavelength: list
                list containing wavelength of sinusoid at each time step
            :param load: np.array
                numpy array containing vertical normal stresses at the base of the lithosphere

        Returns:
            :return cl_deflection:
            :return cl_sigmaBend:
            :return cl_sigmaMembrane:
        """
        cl_mag1 = 0  # maximum stress on left cantilever
        cl_mag2 = 0  # maximum stress on right cantilever

        cl_len1 = self.xstp * zeros[0]  # length of left cantilever
        cl_len2 = self.xstp * (self.xnum - zeros[1])  # length of right cantilever

        cl_deflection = 0
        cl_sigmaBend = 0
        cl_sigmaMembrane = 0
        return cl_deflection, cl_sigmaBend, cl_sigmaMembrane

    def curveFit(self, inFile, saveCurveFit=False):
        """Fit a sinusoid to vertical normal stresses.

        Args:
            :param inFile: str
                path to vertical normal stress data in csv or txt
            :param saveCurveFit: bool
                save curve fit params to hdf5 or not

        Returns:
            :return params: np.array
                numpy array containing curve fit parameters
            :return zeros: np.array
                numpy array containing coordinates of zero crossings
            :return wavelength: list
                list containing wavelength of sinusoid at each time step
            :return time: list
                list containing time steps
            :return vertNormStress: np.array
                numpy array containing vertical normal stresses at the base of the lithosphere
        """
        vertNormStress = np.genfromtxt(inFile, delimiter=',')

        # number of data points
        N = vertNormStress.shape[1]
        xlin = np.linspace(0, self.xstp * self.xnum, self.xnum)

        wavelength = []
        time = []
        params = np.empty((vertNormStress.shape[0], 4))
        zeros = np.empty((vertNormStress.shape[0], 2))
        for t in range(vertNormStress.shape[0]):
            # fit sinusoid
            params[t], paramsCov = optimize.curve_fit(self.sinusoid, xlin, vertNormStress[t],
                                                      p0=[-2000000, 0.000001, 0, 0])

            # plt.plot(xlin, params[t, 0] * np.sin(params[t, 1]*xlin + params[t, 2]) + params[t, 3])
            # plt.plot(xlin, vertNormStress[t], '.')
            # plt.show()

            # get two zero crossings
            # assumes single cycle of sinusoid
            zeros[t, 0] = (np.arcsin(-params[t, 3] / params[t, 0]) - params[t, 2]) / params[t, 1]
            midpoint = ((np.pi / 2) - params[t, 2]) / params[t, 1]
            zeros[t, 1] = midpoint + (midpoint - zeros[t, 0])

            print('Timestep {}'.format(t))
            print('Wavelength: {} km'.format((2 * np.pi / params[t, 1]) * 30))
            print('Amplitude: {}'.format(params[t, 0]))
            print('\n')

            # record wavelength and time step
            wavelength.append((2 * np.pi / params[t, 1]) * 30000)
            time.append(t)

        if saveCurveFit:
            hf = h5py.File('curveFitParams.h5', 'w')
            grp = hf.create_group('params')
            grp.create_dataset('params', data=params)
            grp.create_dataset('zeros', data=zeros)
            grp.create_dataset('wavelength', data=wavelength)
            grp.create_dataset('time', data=time)
            grp.create_dataset('stress', data=vertNormStress)
            hf.close()

        return params, zeros, wavelength, time, vertNormStress

    def readParams(self, inFile):
        """Read an hdf5 containing curve fit parameters.

        Args:
            :param inFile: str
                hdf5 containing curve fit parameters

        Returns:
            :return params: np.array
                numpy array containing curve fit parameters
            :return zeros: np.array
                numpy array containing coordinates of zero crossings
            :return wavelength: list
                list containing wavelength of sinusoid at each time step
            :return time: list
                list containing time steps
            :return data: np.array
                numpy array containing vertical normal stresses at the base of the lithosphere
        """
        hf = h5py.File(inFile, 'r')
        hfp = hf['params']
        # read datasets
        params = hfp['params'][:]
        zeros = hfp['zeros'][:]
        time = hfp['time'][:]
        wavelength = hfp['wavelength'][:]
        data = hfp['stress'][:]
        hf.close()

        return params, zeros, wavelength, time, data

    def readStress(self, inFile):
        """Read vertical normal stresses from a file.

        Params:
            :param inFile: str
                string to file containing vertical normal stresses at each timestep

        Returns:
            :return vertNormStress: np.array
                numpy array containing vertical normal stresses at each timestep
        """
        vertNormStress = np.genfromtxt(inFile, delimiter=',')

        return vertNormStress

    def sinusoid(self, x, a, b, c, d):
        """Return a sinusoid with a given amplitude, frequency, phase shift, and vertical shift.

        Params:
            :param x: np.array
                domain of the function
            :param a: float
                amplitude
            :param b: float
                frequency
            :param c: float
                phase shift
            :param d: float
                vertical shift

        Returns:
            :return sinusoid: np.array
                array containing description of a sinusoid
        """

        return a * np.sin(b * x + c) + d

    def deflectionGreensFunction(self, x):
        """Compute the analytical Green's function for deflection due to a line load on an infinite plate.

        Params:
            :param x: np.array
                numpy array containing x positions of the load

        Returns:
            :return wgreen: np.array
                numpy array containing computed Green's function for deflection
        """
        wgreen = np.zeros(len(x))
        for i in range(len(x)):
            xda = x[i] / self.alpha
            wgreen[i] = ((self.alpha**3) / (8 * self.D)) * np.exp(-xda) * (np.cos(xda) + np.sin(xda))
        return wgreen

    def bendingStressGreensFunction(self, x):
        """Compute the analytical Green's function for bending stress due to a line load on an infinite plate

        Params:
            :param x: np.array
                numpy array containing x positions of the load

        Returns:
            :param sgreen: np.array
                numpy array containing computed Green's function for bending stress
        """
        sgreen = np.zeros((self.ynum, len(x)))
        y = np.linspace(0, self.lithH, self.ynum) - (self.lithH / 2)

        for i in range(len(y)):
            sgreen[i] = (-self.E * self.lamb**3) / ((1 - self.v**2) * (self.rhom - self.rhoin) * self.g) * \
                       np.exp(-self.lamb * x[:]) * (np.cos(self.lamb * x[:]) - np.sin(self.lamb * x[:])) * y[i]
        return sgreen


if __name__ == '__main__':
    main()
