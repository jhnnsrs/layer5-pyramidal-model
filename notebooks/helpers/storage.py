import pickle

def pickleTrace(trace, append="", name=None, withdoc=False):
    if not name:
        name = trace["cellbuilder"].__name__ 
        if withdoc:
            name += " - " + trace["cellbuilder"].__doc__
            
    filename = name
    filename += " - Amp " + str(trace["params"]["amp"])
    filename += append
    filename += ".pickle"
    path = "pickles/{0}".format(filename)
    with open(path, "wb") as handle:
        pickle.dump(trace, handle)
            
    return path

def unpickleTrace(path):
    with open(path, "rb") as handle:
        b = pickle.load(handle)
    return b