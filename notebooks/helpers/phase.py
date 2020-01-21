import copy

from matplotlib import pyplot
import numpy as np


def phasePlanePlot(tracein, plot=True):
    copytrace = copy.copy(tracein)
    time = copytrace["t"].max() - copytrace["t"].min()
    hertz = copytrace["t"].shape[0] / time
    dT = 1 / hertz

    vgradient = np.gradient(copytrace["v"], copytrace["t"])  # Is this correct?

    if plot is True:
        fig = pyplot.figure(figsize=(20, 12))  # Default figsize is (8,6)
        ax = pyplot.subplot(131)
        voltage = copytrace["v"]
        ax.plot(voltage, vgradient)
        ax.set_xlim(-80, 80)
        ax.set_ylim(-50, 500)
        ax.set_title("Phase Plot")
        ax.set_ylabel('dt/dt [mV]')
        ax.set_xlabel('v [mV]')


