from neuron import h, gui
from matplotlib import pyplot


def attach_current_clamp(cell, delay=5, dur=1, amp=.1, loc=1):
    """Attach a current Clamp to a cell.

    :param cell: Cell object to attach the current clamp.
    :param delay: Onset of the injected current.
    :param dur: Duration of the stimulus.
    :param amp: Magnitude of the current.
    :param loc: Location on the dendrite where the stimulus is placed.
    """
    stim = h.IClamp(cell.dend(loc))
    stim.delay = delay
    stim.dur = dur
    stim.amp = amp
    return stim

def set_recording_vectors(cell):
    """Set soma, dendrite, and time recording vectors on the cell.

    :param cell: Cell to record from.
    :return: the soma, dendrite, and time vectors as a tuple.
    """
    soma_v_vec = h.Vector()   # Membrane potential vector at soma
    dend_v_vec = h.Vector()   # Membrane potential vector at dendrite
    t_vec = h.Vector()        # Time stamp vector
    soma_v_vec.record(cell.soma(0.5)._ref_v)
    dend_v_vec.record(cell.dend(0.5)._ref_v)
    t_vec.record(h._ref_t)
    return soma_v_vec, dend_v_vec, t_vec

def simulate(tstop=25):
    """Initialize and run a simulation.

    :param tstop: Duration of the simulation.
    """
    h.tstop = tstop
    h.run()

def show_output(soma_v_vec, dend_v_vec, t_vec, new_fig=True):
    """Draw the output.

    :param soma_v_vec: Membrane potential vector at the soma.
    :param dend_v_vec: Membrane potential vector at the dendrite.
    :param t_vec: Timestamp vector.
    :param new_fig: Flag to create a new figure (and not draw on top
            of previous results)
    """
    if new_fig:
        pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
    soma_plot = pyplot.plot(t_vec, soma_v_vec, color='black')
    dend_plot = pyplot.plot(t_vec, dend_v_vec, color='red')
    pyplot.legend(soma_plot + dend_plot, ['soma', 'dend(0.5)'])
    pyplot.xlabel('time (ms)')
    pyplot.ylabel('mV')
