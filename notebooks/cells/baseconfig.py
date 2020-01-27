

class BaseConfig(object):

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            self.__setattr__(key,value)

        # DISTRIBUTIONS:

        self.distributions = {}
        self.synapsList = {}

        self.setSynapses()
        self.setDistributions()

    def setSynapses(self):
        raise NotImplementedError

    def setDistributions(self):
        raise NotImplementedError

    def plotDistribution(self, name="nax_ais"):
        import matplotlib.pyplot as plt

        plt.plot(np.linspace(0,100,100),self.distributions[name])
        plt.show()

    def plotDistributions(self, normalized=False):
        import matplotlib.pyplot as plt

        for key, value in self.distributions.items():
            plt.plot(np.linspace(0, 100, 100), value if not normalized else value/value.max(), label=f"{key}")

        plt.legend()
        plt.show()

    @property
    def synapses(self):
        return self.synapsList