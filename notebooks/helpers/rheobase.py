import numpy as np
from matplotlib import pyplot
from neuron import h

from helpers import generic
from helpers.generic import stimulate


def calculateRheobase(cellbuilder, start= -0.2, precisiondigits = 2, intialincrease = 0.025 ,ranget = 3, notlower=True, plot=True, apthreshold = 30, doprint = False):
    
    guess = start
    increase = intialincrease
    innerrun = 20
    latestunder = None
    for run in range(0, precisiondigits):
        
        if doprint: print("Intial guess",guess)
        thisrange = abs(guess)*ranget/10**(run+1)
        lower = guess - thisrange
        updatedlower = guess if lower < start and notlower else lower
        parammaps = [{"delay":100,"dur":500,"amp": updatedlower + thisrange*i} for i in range(0,innerrun)]
        if doprint: print("Running in range ", thisrange, "from ",parammaps[0]["amp"], "to ", parammaps[innerrun-1]["amp"])
        runparams = None

        for key, param in enumerate(parammaps):

            trace = generic.stimulate(cellbuilder, param)
            runparams = param
            if doprint: print(key, innerrun, param["amp"], trace["aps"])
            if (trace["aps"].shape[0] > 0):
                if doprint: print(trace["aps"])
                print("Detected first AP at", runparams["amp"] ," [mA]")
                break
                
            latestunder = param
            if key == innerrun - 1:
                raise Exception("We did not detect anything")
        
        guess = latestunder["amp"]
        increase = increase/10**run

    undertrace = stimulate(cellbuilder, latestunder, apthreshold)
    overtrace = stimulate(cellbuilder, runparams, apthreshold)

    if plot is True:
        fig = pyplot.figure(figsize=(15, 9))  # Default figsize is (8,6)
        ax = pyplot.subplot(121)
        ax.plot(undertrace["t"], undertrace["v"], label="First Undertrace Amp {0:.3f}".format(undertrace["params"]["amp"]))
        ax.plot(overtrace["t"], overtrace["v"], label="First Ap Trace Amp {0:.3f}".format(overtrace["params"]["amp"]))
        ax.set_xlabel('time (ms)')
        ax.set_ylabel('mV')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax.set_title("Voltage Graphs of Aps " + cellbuilder.__name__ + " Cell")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=4)
        pyplot.show()


        
    return (undertrace, overtrace)

