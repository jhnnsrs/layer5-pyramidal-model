import copy

import numpy as np
import pandas as pd
from matplotlib import pyplot
from neuron import h


def stimulate(cellbuilder, param, apthreshold = 30, synmodifier = 0.03, interval = 7):
    cell = cellbuilder()
    singlepulse = h.IClamp(cell.soma(0.5))
    singlepulse.delay = param["delay"]
    singlepulse.dur = param["dur"]
    singlepulse.amp = param["amp"]
    
    apc = h.APCount(cell.soma(0.5))
    apc.thresh = apthreshold


    netstim0 = h.NetStim()
    netstim0.interval = interval
    netstim0.number = int(param["dur"] / interval)
    netstim0.start = param["delay"]

    netstim0.noise = 0

    synapses = []
    connections = []
    for synapse in cell.c.synapses:
        if synapse["type"] == "Exp2Syn":
            syn = h.Exp2Syn(cell.ais(synapse["loc"]))
            syn.tau1 = synapse["tau1"]
            syn.tau2 = synapse["tau2"]
            syn.e = synapse["e"]
            synapses.append(syn)
            conn = h.NetCon(netstim0, syn)
            conn.weight[0] = synapse["weight"] * synmodifier

            connections.append(conn)

    v_vec = h.Vector()
    a_vec = h.Vector()# Membrane potential vector
    t_vec = h.Vector()
    v_vec.record(cell.soma(0.5)._ref_v)
    t_vec.record(h._ref_t)
    apc.record(a_vec)
    h.tstop = 1000
    h.run()
    return {"v": np.array(v_vec), "t": np.array(t_vec) ,"aps": np.array(a_vec), "cellbuilder": cellbuilder, "params": param}

def frameTrace(trace):
    exclude_keys = ['t', 'v/t',"v","aps"]
    new_d = {k: trace[k] for k in set(list(trace.keys())) - set(exclude_keys)}
    return pd.DataFrame.from_dict(new_d, orient="index")
        
        

def timetoframe(timeline, timepoint):
    time = timeline.max() - timeline.min()
    hertz = timeline.shape[0]/time
    return int(hertz * timepoint)
    
    
def frametotime(timeline, frame):
    framecount = timeline.shape[0]
    starttime = timeline.min()
    endtime = timeline.max()
    timecount = endtime - starttime
    timepoint = starttime + (timecount/framecount)*frame
    return timepoint
        
    
def calculateInputResistance(cellbuilder, params = None, mult=-0.005, iterations=10 , startoffset=300, endoffset=100 , delay=100, duration=500, plot=True):
    if not params:
        params = [{"delay":delay,"dur":duration,"amp": 0 + mult*i} for i in range(0,iterations)]
    
    outputs = []
        
    for param in params:
        outputs.append(stimulate(cellbuilder, param, apthreshold=30))
     
    values = []
    for traces in outputs:
        if traces["aps"].shape[0] > 1 : 
            print("This Trace Contains APS do not regard", traces["params"])
            break
            
        startx = traces["params"]["delay"]
        durationx = traces["params"]["dur"]
        
        vcurvestable = traces["v"][timetoframe(traces["t"],startx+startoffset):timetoframe(traces["t"],startx+durationx-endoffset)]
        vcurvemaxdef = traces["v"][timetoframe(traces["t"],startx):timetoframe(traces["t"],startx+startoffset)]
        stablevalue = vcurvestable.mean()
        maxvalue = vcurvemaxdef.min() if mult < 0 else vcurvemaxdef.max()
        values.append([stablevalue,maxvalue,traces["params"]["amp"]])
    
    
    resistances = []
    baselinestable = values[0][0]
    others = values[1:]
    for pair in others:
            stablev, maxvalue, amp = pair
            resistances.append([amp,abs((stablev - baselinestable))/abs(amp),abs((maxvalue - baselinestable))/abs(amp)])
            
    
    if plot is True:
        fig = pyplot.figure(figsize=(15,9)) # Default figsize is (8,6)
        ax = pyplot.subplot(121)
        for trace in outputs:
            ax.plot(trace["t"], trace["v"], label="Amp {0:.3f}".format(trace["params"]["amp"]))
        ax.axvline(delay+startoffset, label="Start Mean", linestyle="-")
        ax.axvline(delay+duration-endoffset, label="End Mean" , linestyle="-")
        ax.set_xlabel('time (ms)')
        ax.set_ylabel('mV')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax.set_title("Voltage Graphs of " + cellbuilder.__name__ + " Cell")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=4)
        pyplot.show()

        

        ax2 = pyplot.subplot(122)
        ax2.plot(
            [item[0] for item in resistances],
            [item[1] for item in resistances], label="Stable"
            )
        ax2.plot(
            [item[0] for item in resistances],
            [item[2] for item in resistances], label="MaxDef"
            )

        resistancemeanstable = np.array([item[1] for item in resistances]).mean()
        resistancemeanmax = np.array([item[2] for item in resistances]).mean()
        ax2.set_title("Input Resistances of " + cellbuilder.__name__ + " Cell")
        ax2.axhline(resistancemeanstable, color="green", label="Stable Mean {0:.2f} mOhm".format(resistancemeanstable))
        ax2.axhline(resistancemeanmax, color="red", label="MaxDef Mean {0:.2f} mOhm".format(resistancemeanmax))
        ax2.set_xlabel('mA')
        ax2.set_ylabel('mOhm')

        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=2)

    
    
    return np.array(resistances)

    
def isolateTime(params, timepoint, area = 5):
    
    newtrace = copy.copy(params)
    start = timetoframe(params["t"], timepoint - area)
    end = timetoframe(params["t"], timepoint + area)
    
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



    