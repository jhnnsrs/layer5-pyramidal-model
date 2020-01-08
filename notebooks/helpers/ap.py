from .generic import *

def find_first_nearest_index_for_value(array, value):
        idx = (np.abs(array - value)).argmin()
        return idx


def calculateAPThresholdAndDerivative(tracein):
    ''' Finds the first AP in a trace and checks for its Properties'''
    copytrace = copy.copy(tracein)
    time = copytrace["t"].max() - copytrace["t"].min()
    hertz = copytrace["t"].shape[0]/time
    dT = 1/hertz

    gradient = np.gradient(copytrace["v"],dT) # Is this correct?
    over = gradient > 20
    apstart = np.argmax(over)
    
    if over[apstart] == False: 
        print("ERROR: No AP detected")
        return False
    
    apstarttime = frametotime(copytrace["t"], apstart)
    apstartvalue = copytrace["v"][apstart]
    
    
    # Cutting the Trace to begin with the AP
    apVtrace = copytrace["v"][apstart:]
    apTtrace = copytrace["t"][apstart:]
    
    
    
    # Finding the First Peak after the Threshold
    peakindex = np.argmax(apVtrace)
    peaktime = frametotime(apTtrace, peakindex)
    peakValue = apVtrace[peakindex]
    
    # Finding the Max Slop in between the APStart and the Peak
    slopeV = apVtrace[:peakindex]
    slopeT = apTtrace[:peakindex]
    slopeVgradient = np.gradient(slopeV,dT)
    maxslopeindex = np.argmax(slopeVgradient)
    maxslopetime = frametotime(slopeT, maxslopeindex)
    maxslopevalue = slopeVgradient[maxslopeindex]
    maxslopepoint = slopeV[maxslopeindex]
    
    # Finding the first Dump after the Peak
    apAfterPeakVtrace = apVtrace[peakindex:]
    apAfterPeakTtrace = apTtrace[peakindex:]
    apAfterPeakGradient = np.gradient(apAfterPeakVtrace,dT)
    zero_crossings = np.where(np.diff(np.sign(apAfterPeakGradient)))[0]
 
    hypermaxindex = zero_crossings[0]
    hypermaxtime = frametotime(apAfterPeakTtrace, hypermaxindex)
    hypermaxvalue = apAfterPeakVtrace[hypermaxindex]
    
    

    # Set parameter of ouput
    trace = copytrace
    trace["v/t"] = gradient
    trace["apStartIndex"] = apstart
    trace["firstOver20Time"] = apstarttime
    trace["firstOver20Value"] = apstartvalue
    trace["peakTime"] = peaktime
    trace["peakValue"] = peakValue
    trace["maxSlopeTime"] = maxslopetime
    trace["maxSlopePoint"] = maxslopepoint
    trace["maxSlopeValue"] = maxslopevalue
    trace["hyperMaxTime"] = hypermaxtime
    trace["hyperMaxValue"] = hypermaxvalue
    trace["apEndIndex"] = peakindex + apstart + hypermaxindex #Maybe not 100% accurate
    trace["Ap AHP"] =  abs(apstartvalue) - abs(hypermaxvalue)
    trace["AP MaxSlope"] = maxslopevalue
    
    
    middleAPValue = apstartvalue + (peakValue - apstartvalue)/2
    tenPercentAPValue = apstartvalue + (peakValue - apstartvalue)/10
    ninetyPercentAPValue = apstartvalue + 9*(peakValue - apstartvalue)/10
    

    firstIndexClosestToMiddleApValue= find_first_nearest_index_for_value(apVtrace[:peakindex],middleAPValue)
    firstIndexClosestToTenPercentApValue= find_first_nearest_index_for_value(apVtrace[:peakindex],tenPercentAPValue)
    firstIndexClosestToNinetyPercentValue= find_first_nearest_index_for_value(apVtrace[:peakindex],ninetyPercentAPValue)
    secondIndexClosestToMiddleApValue = find_first_nearest_index_for_value(apVtrace[peakindex:],middleAPValue)
    secondIndexClosestToTenPercentApValue = find_first_nearest_index_for_value(apVtrace[peakindex:],tenPercentAPValue)
    secondIndexClosestToNinetyPercentApValue = find_first_nearest_index_for_value(apVtrace[peakindex:],ninetyPercentAPValue)
    
    midApUpTime = frametotime(apTtrace[:peakindex], firstIndexClosestToMiddleApValue)
    midApDownTime = frametotime(apTtrace[peakindex:], secondIndexClosestToMiddleApValue)
    midApUpValue = apVtrace[:peakindex][firstIndexClosestToMiddleApValue]
    midApDownValue = apVtrace[peakindex:][secondIndexClosestToMiddleApValue]
    
    tenPercentUpTime = frametotime(apTtrace[:peakindex], firstIndexClosestToTenPercentApValue)
    tenPercentDownTime = frametotime(apTtrace[peakindex:], secondIndexClosestToTenPercentApValue)
    ninetyPercentUpTime = frametotime(apTtrace[:peakindex], firstIndexClosestToNinetyPercentValue)
    ninetyPercentDownTime = frametotime(apTtrace[peakindex:], secondIndexClosestToNinetyPercentApValue)
    
    tenPercentUpValue = apVtrace[:peakindex][firstIndexClosestToTenPercentApValue]
    tenPercentDownValue = apVtrace[peakindex:][secondIndexClosestToTenPercentApValue]
    ninetyPercentUpValue = apVtrace[:peakindex][firstIndexClosestToNinetyPercentValue]
    ninetyPercentDownValue = apVtrace[peakindex:][secondIndexClosestToNinetyPercentApValue]
    
    trace["midApUpTime"] = midApUpTime
    trace["midApUpValue"] = midApUpValue
    trace["midApDownTime"] = midApDownTime
    trace["midApDownValue"] = midApDownValue
    trace["tenPercentUpValue"] = tenPercentUpValue
    trace["tenPercentDownValue"] = tenPercentDownValue
    trace["ninetyPercentUpValue"] = ninetyPercentUpValue
    trace["ninetyPercentDownValue"] = ninetyPercentDownValue
    trace["tenPercentUpTime"] = tenPercentUpTime
    trace["tenPercentDownTime"] = tenPercentDownTime
    trace["ninetyPercentUpTime"] = ninetyPercentUpTime
    trace["ninetyPercentDownTime"] = ninetyPercentDownTime
    trace["AP Rise Time"] = ninetyPercentUpTime - tenPercentUpTime
    trace["AP Fall Time"] = tenPercentDownTime - ninetyPercentDownTime
    trace["AP HW"] = midApDownTime - midApUpTime
    
    
    
    
    return trace


def plotApParameters(test, deriv=False):
    
    fig = pyplot.figure(figsize=(20,12)) # Default figsize is (8,6)
    ax = pyplot.subplot(131)
    ax.plot(test["t"], test["v"])
    ax.plot(test["firstOver20Time"],test["firstOver20Value"], marker='+', color= "red", markersize=15, label = "AP Start {:.2f} [ms] {:.2f} [mv]".format(test["firstOver20Time"],test["firstOver20Value"]))
    #pyplot.axvline(x=test["firstOver20Time"], color = "green", linestyle= "--",label ="AP start {:.2f} [ms]".format(test["firstOver20Time"]))
    ax.axhline(y=test["peakValue"], color = "yellow", linestyle= "--",label ="AP Peak {:.2f} [mv]".format(test["peakValue"]))
    ax.axhline(y=test["hyperMaxValue"], color = "red", linestyle= "--",label ="AP AHP {:.2f} [mv]".format(test["Ap AHP"]))
    ax.set_xlabel('time (ms)')
    ax.set_ylabel('mV | del Mv / del ms')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    ax.set_title("AP Parameters of " + test["cellbuilder"].__name__ + " Cell")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
              fancybox=True, shadow=True, ncol=2)

    ax2 = pyplot.subplot(132)

    ax2.plot(test["t"][test["apStartIndex"]:test["apEndIndex"]],test["v"][test["apStartIndex"]:test["apEndIndex"]])
    ax2.plot([test["midApUpTime"],test["midApDownTime"]], [test["midApUpValue"],test["midApDownValue"]], label ="AP HW  {:.2f} [ms]".format(test["AP HW"]))
    ax2.plot(test["maxSlopeTime"],test["maxSlopePoint"], marker='x', markersize=15, label = "AP max slope  {:.2f} [mV/ms]".format(test["maxSlopeValue"]))
    ax2.plot(test["firstOver20Time"],test["firstOver20Value"], marker='+',markersize=15, label = "AP Start {:.2f} [ms] {:.2f} [mv]".format(test["firstOver20Time"],test["firstOver20Value"]))
    #pyplot.axvline(x=test["firstOver20Time"], color = "green", linestyle= "--",label ="AP start {:.2f} [ms]".format(test["firstOver20Time"]))
    ax2.axhline(y=test["peakValue"], color = "yellow", linestyle= "--",label ="AP Peak {:.2f} [mv]".format(test["peakValue"]))
    ax2.axhline(y=test["hyperMaxValue"], color = "red", linestyle= "--",label ="AP AHP {:.2f} [mv]".format(test["Ap AHP"]))
    ax2.axhline(y=test["tenPercentUpValue"], color = "green", linestyle= "-.",label ="AP 10% {:.2f} [mv]".format(test["tenPercentUpValue"]))
    ax2.axhline(y=test["ninetyPercentUpValue"], color = "blue", linestyle= "--",label ="AP 90% {:.2f} [mv]".format(test["ninetyPercentUpValue"]))
    ax2.set_xlabel('time (ms)')
    ax2.set_ylabel('mV | del Mv / del ms')

    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    ax2.set_title("AP Isolation of " + test["cellbuilder"].__name__ + " Cell")
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
              fancybox=True, shadow=True, ncol=2)
    
    if deriv:
        ax3 = pyplot.subplot(133)

        ax3.plot(test["t"],test["v/t"])
        ax3.plot(test["maxSlopeTime"],test["maxSlopeValue"], marker='x', markersize=15, label = "AP max slope {:.2f} [mV/ms]".format(test["maxSlopeValue"]))
        ax3.plot(test["firstOver20Time"],20, marker='+',markersize=15, label = "AP Start {:.2f} [ms] {:.2f} [mv]".format(test["firstOver20Time"],test["firstOver20Value"]))
        #pyplot.axvline(x=test["firstOver20Time"], color = "green", linestyle= "--",label ="AP start {:.2f} [ms]".format(test["firstOver20Time"]))
        ax3.axhline(y=test["maxSlopeValue"], color = "yellow", linestyle= "--",label ="AP Peak {:.2f} [mv]".format(test["peakValue"]))
        ax3.axhline(y=20, color = "red", linestyle= "--",label ="AP Threshold Defintion {0} [mv/ms]".format(20))
        ax3.set_xlabel('time (ms)')
        ax3.set_ylabel(r'$\frac{\Delta  mV}{\Delta  ms}$')

        box = ax3.get_position()
        ax3.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax3.set_title("First Derivative of " + test["cellbuilder"].__name__ + " Cell")
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),

                   fancybox=True, shadow=True, ncol=1)
    pyplot.show()


    
def findFirstApWithinMS(cellbuilder,ampstep = 0.025,ms = 15,steps = 20,delay = 100,dur = 500):
    foundtrace = None
                        
    baselinetrace = { "dur": dur, "delay": delay, "amp": 0}
    baselinetrace = stimulate(cellbuilder, baselinetrace)
    timetrace = baselinetrace["t"]
    for i in range(0,steps):
        under = { "dur": dur, "delay": delay, "amp": ampstep * i}
        trace = stimulate(cellbuilder, under)
        if trace["aps"].size > 0:
            if trace["aps"][0] < delay + ms:
                foundtrace = trace
                break

    return foundtrace
            