from neuron import h, gui
from matplotlib import pyplot
import numpy as np
import pandas as pd
import copy

def stimulate(cellbuilder, param):
    cell = cellbuilder()
    singlepulse = h.IClamp(cell.soma(0.5))
    singlepulse.delay = param["delay"]
    singlepulse.dur = param["dur"]
    singlepulse.amp = param["amp"]
    
    apc = h.APCount(cell.soma(0.5))
    apc.thresh = 30

    v_vec = h.Vector()
    a_vec = h.Vector()# Membrane potential vector
    t_vec = h.Vector()  
    v2_vec = h.Vector()# Time stamp vector
    v_vec.record(cell.soma(0.5)._ref_v)
    t_vec.record(h._ref_t)
    apc.record(a_vec)
    h.tstop = 800
    h.run()
    return {"v": np.array(v_vec), "t": np.array(t_vec) ,"aps": np.array(a_vec)}

def frameTrace(trace):
    exclude_keys = ['t', 'v/t',"v","aps"]
    new_d = {k: trace[k] for k in set(list(trace.keys())) - set(exclude_keys)}
    return pd.DataFrame.from_dict(new_d, orient="index")
        
        

def timetoframe(params, timepoint):
    time = params["t"].max() - params["t"].min()
    hertz = params["t"].shape[0]/time
    return int(hertz * timepoint)
    
    
def frametotime(timeline, frame):
    framecount = timeline.shape[0]
    starttime = timeline.min()
    endtime = timeline.max()
    timecount = endtime - starttime
    timepoint = starttime + (timecount/framecount)*frame
    return timepoint
        
    
def calculateInputResistance(cellbuilder, params = None):
    if not params:
        params = [{"delay":100,"dur":500,"amp": 0 + i*0.01} for i in range(0,7)]
    
    outputs = []
        
    for param in params:
        cell = cellbuilder()
        singlepulse = h.IClamp(cell.soma(0.5))
        singlepulse.delay = param["delay"]
        singlepulse.dur = param["dur"]
        singlepulse.amp = param["amp"]

        apc = h.APCount(cell.soma(0.5))
        #apc.thresh = 40

        v_vec = h.Vector()
        a_vec = h.Vector()# Membrane potential vector
        t_vec = h.Vector()  
        v2_vec = h.Vector()# Time stamp vector
        v_vec.record(cell.soma(0.5)._ref_v)
        t_vec.record(h._ref_t)
        apc.record(a_vec)
        
        
        h.tstop = 800
        h.run()
        outputs.append({"v": np.array(v_vec), "t": np.array(t_vec) ,"aps": np.array(a_vec), "params": param})
     
    maxValues = []
    for traces in outputs:
        if traces["aps"].shape[0] > 1 : 
            print("This Trace Contains APS do not regard", traces["params"])
            break
            
        start = traces["params"]["delay"]
        duration = traces["params"]["dur"]
        
        maxValues.append([traces["v"][timetoframe(traces,start):timetoframe(traces,start+duration)].max(),traces["params"]["amp"]])
    
    resistances = []
    baseline = maxValues[0][0]
    for pair in maxValues[1:]:
        vmax, amp = pair
        resistances.append((vmax - baseline)/amp)
        
    return np.array(resistances)

    
def isolateTime(params, timepoint, area = 5):
    
    newtrace = copy.copy(params)
    start = timetoframe(params, timepoint - area)
    end = timetoframe(params, timepoint + area)
    
    newtrace["t"] = params["t"][start:end]
    newtrace["v"] = params["v"][start:end]
    return newtrace


def plotTrace(trace):

    pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
    pyplot.plot(trace["t"], trace["v"])
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')
    pyplot.show()
    
def plotTraces(traces):

    pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
    for trace in traces:
        pyplot.plot(trace["t"], trace["v"])
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')
    pyplot.show()
    