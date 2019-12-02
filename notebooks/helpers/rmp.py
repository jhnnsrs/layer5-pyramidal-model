from neuron import h, gui
from matplotlib import pyplot
import numpy as np
from cells import layer_v_with_ais_test 
from neuron import h, gui
from .record import *
import numpy as np
import copy
from importlib import reload
from .rheobase import *
from .generic import *

def find_first_nearest_index_for_value(array, value):
        idx = (np.abs(array - value)).argmin()
        return idx

    
def calculateRMP(cellbuilder,ampstep = 0.025,ms = 15,steps = 20,delay = 100,dur = 500):
    foundtrace = None
                        
    baselinetrace = { "dur": dur, "delay": delay, "amp": 0}
    baselinetrace = stimulate(cellbuilder, baselinetrace)
    rmp = baselinetrace["v"][timetoframe(baselinetrace["t"],delay): timetoframe(baselinetrace["t"],delay+100)].mean()

    return rmp
            