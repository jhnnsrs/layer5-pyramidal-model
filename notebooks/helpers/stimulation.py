from neuron import h, gui
from matplotlib import pyplot
import numpy as np

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
