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
    parser.add_argument('-sx', '--surface', type=str, required=False, nargs=2, default=None,
                        help='Path to text file containing x-coordinates of the surface from the flow model.')
    args = parser.parse_args()

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

    # expand dimension if we have only a single timestep
    if lithLoad.ndim == 1:
        lithLoad = np.reshape(lithLoad[1:], (1, lithLoad.shape[0]-1))
    stepmax = lithLoad.shape[0]

    if args.surface is not None:
        surf = flex.readSurface(args.surface)

    dArr = np.linspace(1e+24, 1e+25, 10)
    error = np.empty((dArr.shape[0], stepmax))
    for i, d in enumerate(dArr):
        flex.D = d
        print('Step {} --- {}'.format(i, d))
        # compute flexure and bending stress for each time step
        sigxx = np.empty((stepmax, int(flex.lithH / flex.ystp), flex.xnum))
        deflection = np.empty((stepmax, flex.xnum))
        surface = np.empty((stepmax, flex.xnum))

        for tStep in range(stepmax):
            # define green's functions for deflection and bending stress
            greenDef = flex.deflectionGreensFunction(flex.x)
            greenStress = flex.bendingStressGreensFunction(flex.x)

            # compute deflection
            compDeflection = flex.convGreenLoad(greenDef, lithLoad[tStep])
            deflection[tStep] = signal.resample(compDeflection, num=flex.xnum)

            # compute bending stress
            compSigxx = np.empty((flex.ynum, 301))
            for d in range(flex.ynum):
                compSigxx[d] = flex.convGreenLoad(greenStress[d], lithLoad[tStep])

            # resample x dimension of computed stress
            sigxx_dsampx = np.empty((flex.ynum, flex.xnum))
            sigxx_dsampx = signal.resample(compSigxx, num=flex.xnum, axis=1)
            # resample y dimension of computed stress
            sigxx[tStep] = signal.resample(sigxx_dsampx, num=int(flex.lithH/flex.ystp), axis=0)

            # resample surface deflection from flow model
            surface[tStep] = signal.resample(surf[tStep, 1], flex.xnum)

            error[i, tStep] = flex.rmsError(deflection[tStep], surface[tStep])

    error = np.swapaxes(error, 0, 1)
    plt.imshow(error, cmap='coolwarm')
    plt.colorbar()
    plt.ylabel('Time Step')
    plt.show()

    tStepBest = np.empty(error.shape[0])
    tStepBestError = np.empty(error.shape[0])
    for i in range(tStepBest.shape[0]):
        tStepBest[i] = np.where(error[i] == np.amin(error[i]))[0]
        tStepBestError[i] = np.amin(error[i])

    plt.subplot(1, 2, 1)
    plt.plot(tStepBest)
    plt.subplot(1, 2, 2)
    plt.plot(tStepBestError)
    plt.show()



class Flexure:
    def __init__(self):
        """Initialize class for computing deflection and stresses of a discretized or sinusoidal load."""
        # simulation parameters
        self.xstp = 20000  # m -- x spacing
        self.ystp = 20000  # m -- y spacing
        self.xnum = 101  # number of nodes in x
        self.lithH = 100000  # m -- lithospheric thickness
        self.rhom = 3300  # kg/m**3 -- mantle density
        self.rhoin = 1000  # kg/m**3 -- infill density
        self.g = 9.81  # m/s**2 -- gravitational acceleration

        # material properties
        self.E = 6.5e+10  # Pa -- Young's modulus
        self.v = 0.25  # Dimensionless -- Poisson's ratio
        self.I = (self.xstp * self.xnum) * (self.lithH ** 3 / 3)  # m**4 -- second moment of area
        self.D = (self.E * (self.lithH) ** 3) / (12 * (1 - self.v ** 2))  # Nm -- flexural rigidity
        self.alpha = ((4 * self.D) / ((self.rhom - self.rhoin) * self.g))**0.25  # flexural parameter
        self.lamb = 1 / self.alpha  # effective wavelength of the load

        # calculation parameters
        self.lowBound = 0  # m -- lower bound in x for where we calculate flexure and stress
        self.upBound = 2000000  # m -- upper bound in x for where we calculate flexure and stress
        self.spacing = 20000  # m -- spacing between sampled points
        self.numPoints = int((self.upBound - self.lowBound) / self.spacing) + 1  # number of points in x
        self.x = np.linspace(self.lowBound, self.upBound, self.numPoints)  # range of points in x
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

        return ss_deflection, ss_sigmaBend

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

            # get two zero crossings
            # assumes single cycle of sinusoid
            zeros[t, 0] = (np.arcsin(-params[t, 3] / params[t, 0]) - params[t, 2]) / params[t, 1]
            midpoint = ((np.pi / 2) - params[t, 2]) / params[t, 1]
            zeros[t, 1] = midpoint + (midpoint - zeros[t, 0])

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

    def readSurface(self, surface):
        """Read the x and y coordinates of the surface of the lithosphere from the flow model.

        Parameters:
            :param surface: list
                path to the files containing x and y coordinates of the surface of the lithosphere.

        Returns:
            :return surf: np.array
                numpy array containing x and zeroed y coordinates of the surface of the lithosphere.
        """
        surfX = np.genfromtxt(surface[0], delimiter=',')
        surfY = np.genfromtxt(surface[1], delimiter=',')

        surf = np.empty((surfX.shape[0], 2, surfX.shape[1]))
        surf[:, 0] = surfX
        surf[:, 1] = self.zeroSurface(surfY)

        return surf

    def readStress(self, inFile):
        """Read vertical normal stresses from a file.

        Params:
            :param inFile: str
                string to file containing vertical normal stresses at each timestep

        Returns:
            :return vertNormStress: np.array
                numpy array containing vertical normal stresses at each timestep
        """
        fType = inFile[-3:]
        if fType == 'txt':
            vertNormStress = -np.genfromtxt(inFile, delimiter=',')
        else:
            hf = h5py.File(inFile)
            vertNormStress = -hf['sigmayy'][:, int(self.lithH / self.ystp)]
            hf.close()

        return vertNormStress

    def zeroSurface(self, yCoords):
        """Convert depth of lithospheric markers to relative changes in depths.

        Parameters:
            :param yCoords: np.array
                numpy array containing y coordinates of the surface deflection profile.

        Returns:
            :return zeroSurface: np.array
                numpy array containing relative movements of the surface of the lithosphere.
        """
        zeroSurface = np.array(yCoords.shape)
        # subtract original depth to give relative change
        zeroSurface = yCoords - yCoords[0]

        return zeroSurface

    def rmsError(self, curve_1, curve_2):
        """Compute the rms error between two curves.

        Parameters:
            :param curve_1: np.array
                numpy array containing deflection values
            :param curve_2: np.array
                numpy array containing deflection values

        Returns:
            :return error: float
                float containing rms error between the two curves.
        """
        return np.sqrt(np.sum((curve_1-curve_2)**2))


if __name__ == '__main__':
    main()
