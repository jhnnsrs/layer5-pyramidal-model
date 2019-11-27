from neuron import h, gui
from matplotlib import pyplot
import numpy as np

def calculateRheobase(cellbuilder, start= -0.2, precisiondigits = 3, intialincrease = 0.025 ,ranget = 3, notlower=True):
    
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
        
        
    return (runparams, latestunder)


    