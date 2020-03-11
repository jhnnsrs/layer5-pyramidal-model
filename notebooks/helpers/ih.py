from matplotlib import pyplot
from neuron import h
from sklearn.linear_model import LinearRegression
import numpy as np

from helpers import generic
from helpers.generic import timetoframe


def calculateIhSagTraces(cellbuilder):
    baselineparams = {"dur": 500, "delay": 100, "amp": 0}
    baseline = generic.stimulate(cellbuilder, baselineparams)

    params = [{"dur": 500, "delay": 100, "amp": 0 - i * 0.025} for i in range(1, 10)]
    traces = []

    for param in params:
        traces.append(generic.stimulate(cellbuilder, param))

    return [baseline, traces]


def calculateStableAndMinVoltage(baseline, traces, startoffset=200, endoffset=100, plot=True):
    results = []
    baseline = baseline
    baselineamp = baseline["params"]["amp"]
    baselinedelay = baseline["params"]["delay"]
    baselinedur = baseline["params"]["dur"]

    baselinestablebegin = generic.timetoframe(baseline["t"], baselinedelay + startoffset)
    baselinestableend = generic.timetoframe(baseline["t"], baselinedelay + baselinedur - endoffset)
    delayoffset = generic.timetoframe(baseline["t"], baselinedelay)

    baselineminvoltage = baseline["v"][delayoffset:].min()
    baselinestablevoltage = baseline["v"][baselinestablebegin:baselinestableend].mean()

    for trace in traces:
        amp = trace["params"]["amp"]
        delay = trace["params"]["delay"]
        dur = trace["params"]["dur"]

        stablebegin = generic.timetoframe(trace["t"], delay + startoffset)
        stableend = generic.timetoframe(trace["t"], delay + dur - endoffset)
        minvoltage = trace["v"][delayoffset:].min()  # As we are sending negativ impulses

        stablevoltage = trace["v"][stablebegin:stableend].mean()

        results.append({"amp": amp,
                        "min": minvoltage,
                        "stable": stablevoltage,
                        "mintobase": minvoltage - baselineminvoltage,
                        "stabletobase": stablevoltage - baselinestablevoltage
                        })

    if plot is True:

        ystable = [-el["stabletobase"] for el in results]
        xmin = [-el["mintobase"] for el in results]

        ymin = [-el["mintobase"] for el in results]
        xamp = [-el["amp"] for el in results]

        # This will plot a function with slope 1
        minmin = np.array(xmin).max()
        maxmin = np.array(xmin).min()
        x = np.linspace(minmin, maxmin, 100)

        from matplotlib import pyplot
        fig = pyplot.figure(figsize=(20, 8))  # Default figsize is (8,6)

        ax2 = pyplot.subplot(131)  # Default figsize is (8,6)
        ax2.plot(xamp, ymin, marker="o")
        ax2.set_title("Voltage Deflection")
        ax2.set_xlabel('injected current [mA]')
        ax2.set_ylabel('min Voltage [mV]')

        model = LinearRegression()
        lx = np.array(xmin).reshape((-1, 1))
        ly = np.array(ystable)
        model.fit(lx, ly)
        r_sq = model.score(lx, ly)
        y_pred = model.intercept_ + model.coef_ * xmin

        ax = pyplot.subplot(132)
        ax.plot(xmin, ystable, marker="o", label="Ih Slope")
        ax.plot(xmin, y_pred, marker="o", label=f"Linear Slope {model.coef_}")
        ax.set_title("Ih Slope")
        ax.plot(x, x, label="No Ih Slope x=y")
        ax.set_xlabel('max voltage def [mV]')
        ax.set_ylabel('stable voltage def [mV]')
        ax.legend()

        ax3 = pyplot.subplot(133)
        ax3.plot(baseline["t"], baseline["v"], label="Baseline".format(baseline["params"]["amp"]))
        for trace in traces:
            ax3.plot(trace["t"], trace["v"], label="Amp {0:.3f}".format(trace["params"]["amp"]))
        ax3.axvline(baselinedelay + startoffset, label="Start Mean", linestyle="-")
        ax3.axvline(baselinedelay + baselinedur - endoffset, label="End Mean", linestyle="-")
        ax3.set_xlabel('time (ms)')
        ax3.set_ylabel('mV')
        box3 = ax3.get_position()
        ax3.set_position([box3.x0, box3.y0 + box3.height * 0.1, box3.width, box3.height * 0.9])
        ax3.set_title("Voltage Graphs of")
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=4)
        pyplot.show()

    return results


def calculateInputResistanceAndIhSag(cellbuilder, traces= None, params = None, mult=-0.025, iterations=10 , startoffset=300, endoffset=100 , delay=100, duration=500, plot=True):
    if not params:
        params = [{"delay": delay, "dur": duration, "amp": 0 + mult * i} for i in range(0, iterations)]


    if traces is None:
        traces = []
        for param in params:
                traces.append(generic.stimulate(cellbuilder, param))


    baseline = traces[0]

    baselineamp = baseline["params"]["amp"]
    baselinedelay = baseline["params"]["delay"]
    baselinedur = baseline["params"]["dur"]

    baselinestablebegin = generic.timetoframe(baseline["t"], baselinedelay + startoffset)
    baselinestableend = generic.timetoframe(baseline["t"], baselinedelay + baselinedur - endoffset)
    delayoffset = generic.timetoframe(baseline["t"], baselinedelay)

    baselineminvoltage = baseline["v"][delayoffset:].min()
    baselinestablevoltage = baseline["v"][baselinestablebegin:baselinestableend].mean()

    others = traces[1:]

    results = []
    for trace in others:
        amp = trace["params"]["amp"]
        delay = trace["params"]["delay"]
        dur = trace["params"]["dur"]

        stablebegin = generic.timetoframe(trace["t"], delay + startoffset)
        stableend = generic.timetoframe(trace["t"], delay + dur - endoffset)
        minvoltage = trace["v"][delayoffset:].min()  # As we are sending negativ impulses

        stablevoltage = trace["v"][stablebegin:stableend].mean()

        results.append({"amp": amp,
                        "min": minvoltage,
                        "stable": stablevoltage,
                        "mintobase": minvoltage - baselineminvoltage,
                        "stabletobase": stablevoltage - baselinestablevoltage
                        })


    resistances = []
    for item in results:
        amp = item["amp"]
        mintobase = item["mintobase"]
        stabletobase = item["stabletobase"]
        resistances.append(
            [amp, abs((stabletobase)) / abs(amp), abs((mintobase)) / abs(amp)])

    ystable = [-el["stabletobase"] for el in results]
    xmin = [-el["mintobase"] for el in results]

    ymin = [-el["mintobase"] for el in results]
    xamp = [-el["amp"] for el in results]

    # This will plot a function with slope 1
    minmin = np.array(xmin).max()
    maxmin = np.array(xmin).min()
    x = np.linspace(minmin, maxmin, 100)

    resistancemeanstable = np.array([item[1] for item in resistances]).mean()
    resistancemeanmax = np.array([item[2] for item in resistances]).mean()

    model = LinearRegression()
    lx = np.array(xmin).reshape((-1, 1))
    ly = np.array(ystable)
    model.fit(lx, ly)
    r_sq = model.score(lx, ly)
    y_pred = model.intercept_ + model.coef_ * xmin

    if plot is True:

        from matplotlib import pyplot
        fig = pyplot.figure(figsize=(22, 8))  # Default figsize is (8,6)

        ax2 = pyplot.subplot(141)  # Default figsize is (8,6)
        ax2.plot(xamp, ymin, marker="o")
        ax2.set_title(f"Voltage Deflection of {cellbuilder.__name__} Cell")
        ax2.set_xlabel('injected current [mA]')
        ax2.set_ylabel('min Voltage [mV]')


        ax = pyplot.subplot(142)
        ax.plot(xmin, ystable, marker="o", label="Ih Slope")
        ax.plot(xmin, y_pred, marker="o", label=f"Linear Slope {model.coef_}")
        ax.set_title(f"Ih Slope of {cellbuilder.__name__} Cell")
        ax.plot(x, x, label="No Ih Slope x=y")
        ax.set_xlabel('max voltage def [mV]')
        ax.set_ylabel('stable voltage def [mV]')
        ax.legend()

        ax3 = pyplot.subplot(143)
        ax3.plot(baseline["t"], baseline["v"], label="Baseline".format(baseline["params"]["amp"]))
        for trace in traces:
            ax3.plot(trace["t"], trace["v"], label="Amp {0:.3f}".format(trace["params"]["amp"]))
        ax3.axvline(baselinedelay + startoffset, label="Start Mean", linestyle="-")
        ax3.axvline(baselinedelay + baselinedur - endoffset, label="End Mean", linestyle="-")
        ax3.axvline(baselinedelay, label="Start Max/Def", linestyle="-")
        ax3.set_xlabel('time (ms)')
        ax3.set_ylabel('mV')
        box3 = ax3.get_position()
        ax3.set_position([box3.x0, box3.y0 + box3.height * 0.1, box3.width, box3.height * 0.9])
        ax3.set_title(f"Voltage Graphs of {cellbuilder.__name__} Cell")
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=4)

        ax4 = pyplot.subplot(144)
        ax4.plot(
            [item[0] for item in resistances],
            [item[1] for item in resistances], label="Stable"
        )
        ax4.plot(
            [item[0] for item in resistances],
            [item[2] for item in resistances], label="MaxDef"
        )

        ax4.set_title("Input Resistances of " + cellbuilder.__name__ + " Cell")
        ax4.axhline(resistancemeanstable, color="green", label="Stable Mean {0:.2f} mOhm".format(resistancemeanstable))
        ax4.axhline(resistancemeanmax, color="red", label="MaxDef Mean {0:.2f} mOhm".format(resistancemeanmax))
        ax4.set_xlabel('mA')
        ax4.set_ylabel('mOhm')

        box = ax4.get_position()
        ax4.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=2)


        pyplot.show()



    return {"res_stablemean": resistancemeanstable, "res_maxmean": resistancemeanmax, "ih_slope": model.coef_}