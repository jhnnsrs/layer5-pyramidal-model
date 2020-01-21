import numpy as np
from matplotlib import pyplot
from neuron import h

from helpers.generic import stimulate


def calculateRheobase(cellbuilder, start= -0.2, precisiondigits = 1, intialincrease = 0.025 ,ranget = 3, notlower=True, plot=True):
    
    guess = start
    increase = intialincrease
    innerrun = 20
    latestunder = None
    for run in range(0, precisiondigits):
        
        print("Intial guess",guess)
        thisrange = abs(guess) *ranget/10**run
        lower = guess - thisrange
        updatedlower = guess if lower < start and notlower else lower
        up = guess + thisrange
        print("Running in range ", thisrange, "from ",updatedlower, "to ", up)
        parammaps = [{"delay":100,"dur":500,"amp": updatedlower + (thisrange/innerrun)*i} for i in range(0,innerrun)]
        runparams = None
        for param in parammaps:
            cell = cellbuilder()
            singlepulse = h.IClamp(cell.soma(0.5))
            singlepulse.delay = param["delay"]
            singlepulse.dur = param["dur"]
            singlepulse.amp = param["amp"]

            apc = h.APCount(cell.soma(0.5))
            apc.thresh = 30

            v_vec = h.Vector()
            a_vec = h.Vector() # Action potential vector
            t_vec = h.Vector()  
            v2_vec = h.Vector()# Time stamp vector
            v_vec.record(cell.soma(0.5)._ref_v)
            t_vec.record(h._ref_t)
            apc.record(a_vec)

            h.tstop = 800
            h.run()
            aps = np.array(a_vec)
            runparams = param
            if (aps.shape[0] > 0):
                print("Detected first AP at", runparams["amp"] ," [mv]")
                
                break
                
            latestunder = param
        
        guess = runparams["amp"]
        increase = increase/10**run

    undertrace = stimulate(cellbuilder, latestunder)
    overtrace = stimulate(cellbuilder, runparams)

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


    