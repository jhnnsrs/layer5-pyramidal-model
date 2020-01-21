from matplotlib import pyplot

from helpers import generic
from helpers.generic import timetoframe


    
def calculateRMP(cellbuilder,startoffset=100, endoffset=100,delay = 100,dur = 500, plot=False):
    foundtrace = None
                        
    baselinetraceparams = { "dur": dur, "delay": delay, "amp": 0}
    baselinetrace = generic.stimulate(cellbuilder, baselinetraceparams)
    rmp = baselinetrace["v"][timetoframe(baselinetrace["t"],delay+startoffset): timetoframe(baselinetrace["t"],delay+dur-endoffset)].mean()

    if plot is True:
        fig = pyplot.figure(figsize=(15, 9))  # Default figsize is (8,6)
        ax = pyplot.subplot(121)
        ax.plot(baselinetrace["t"], baselinetrace["v"], label="Amp {0:.3f}".format(baselinetrace["params"]["amp"]))
        ax.axvline(delay + startoffset, label="Start Mean", linestyle="-", color="#00FFFF")
        ax.axvline(delay + dur - endoffset, label="End Mean", linestyle="-", color="#00FFFF")
        ax.axhline(rmp, label=f"RMP {rmp:.3f}", linestyle="-", color="#FF00FF")
        ax.set_xlabel('time (ms)')

        ax.set_ylabel('mV')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax.set_title("RMP of " + cellbuilder.__name__ + " Cell")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=4)
        pyplot.show()

    return rmp
            