# simplified corticospinal cell model (6 compartment) according to Neymotin 2017

from neuron import h
from math import exp, log

from cells.baseconfig import BaseConfig

class Config(BaseConfig):
    Vrest = -88.5366550238
    celsius = 34.0  # for in vitro opt

    # geom properties
    somaL = 48.4123467666
    somaDiam = 28.2149102762
    aisL = 32.00
    axonL = 594.292937602 - aisL
    axonDiam = 1.40966286462
    aisDiam = axonDiam
    apicL = 261.904636003
    apicDiam = 1.5831889597
    bdendL = 299.810775175
    bdendDiam = 2.2799248874

    # passive properties
    axonCap = 1.01280903702
    somaCap = 1.78829677463
    apicCap = 1.03418636866
    bdendCap = 1.89771901209
    rall = 114.510490019
    axonRM = 3945.2107187
    somaRM = 18501.7540916
    apicRM = 10751.193413
    bdendRM = 13123.00174

    # Na, K reversal potentials calculated from BenS internal/external solutions via Nernst eq.
    p_ek = -104.0  # these reversal potentials for in vitro conditions
    p_ena = 42.0
    # h-current
    erev_ih = -37.0  # global
    gbar_h = 0.000140956438043
    h_gbar_tuft = 0.00565  # mho/cm^2 (based on Harnett 2015 J Neurosci)

    # d-current
    gbar_kdmc = 0.000447365630734
    kdmc_gbar_axonm = 20

    # spiking currents
    gbar_nax = 0.0345117294903
    nax_gbar_axonm = 5.0
    gbar_kdr = 0.0131103978049
    kdr_gbar_axonm = 5.0

    # A few kinetic params changed vis-a-vis kdr.mod defaults:
    kdr_vhalfn = 11.6427471384
    gbar_kap = 0.0898600246397
    kap_gbar_axonm = 5.0

    # A few kinetic params changed vis-a-vis kap.mod defaults:
    kap_vhalfn = 32.7885075379
    kap_tq = -52.0967985869
    kap_vhalfl = -59.7867409796  # global!!

    # other ion channel parameters
    cal_gcalbar = 4.41583533572e-06
    can_gcanbar = 4.60717910591e-06
    calginc = 1.0
    h_lambda = 325.0
    kBK_gpeak = 5.09733585163e-05
    kBK_caVhminShift = 43.8900261407
    cadad_depth = 0.119408607923
    cadad_taur = 99.1146852282

    def setDistributions(self):
        pass


###############################################################################
# SPI6 Cell
###############################################################################
class SPI6(object):
    "Simplified Corticospinal Cell Model"

    def __init__(self, config= Config(), x=0, y=0, z=0, ID=0):
        self.x, self.y, self.z = x, y, z
        self.ID = ID
        self.c = config
        self.all_sec = []
        self.set_morphology()
        self.insert_conductances()
        self.set_props()
        self.calc_area()

    def add_comp(self, name):
        self.__dict__[name] = h.Section(name=name)  # ,cell=self)
        self.all_sec.append(self.__dict__[name])

    def calc_area(self):
        self.total_area = 0
        self.n = 0
        for sect in self.all_sec:
            self.total_area += h.area(0.5, sec=sect)
            self.n += 1

    def initial_params(self):
        h.celsius = self.c.celsius
        h.erev_neyih =self.c.erev_ih

    def set_morphology(self):

        self.soma = h.Section(name="soma", cell=self)
        self.axon = h.Section(name="axon", cell=self)
        self.ais = h.Section(name="ais", cell=self)

        self.Bdend = h.Section(name="Bdend", cell=self)
        self.Adend1 = h.Section(name="Adend1", cell=self)
        self.Adend2 = h.Section(name="Adend2", cell=self)
        self.Adend3 = h.Section(name="Adend3", cell=self)


        self.apic = [self.Adend1, self.Adend2, self.Adend3]
        self.basal = [self.Bdend]
        self.alldend = [self.Adend1, self.Adend2, self.Adend3, self.Bdend]
        self.all_sec = [ self.soma, self.ais, self.axon] + self.apic + self.basal + self.alldend


        self.initial_params()
        # Set Geometry here
        self.set_geom()


        self.ais.connect(self.soma, 0.0, 0.0)
        self.axon.connect(self.soma, 1.0, 0.0)
        self.Bdend.connect(self.soma, 0.5, 0.0)  # soma 0.5 to Bdend 0
        self.Adend1.connect(self.soma, 1.0, 0.0)
        self.Adend2.connect(self.Adend1, 1.0, 0.0)
        self.Adend3.connect(self.Adend2, 1.0, 0.0)

    def set_geom(self):
        self.axon.L = self.c.axonL
        self.axon.diam = self.c.axonDiam
        self.soma.L = self.c.somaL
        self.soma.diam = self.c.somaDiam
        self.ais.L = self.c.aisL
        self.ais.diam = self.c.aisDiam

        for sec in self.apic: sec.L, sec.diam = self.c.apicL, self.c.apicDiam

        self.Bdend.L = self.c.bdendL
        self.Bdend.diam = self.c.bdendDiam

    def activeoff(self):
        for sec in self.all_sec: sec.gbar_nax = sec.gbar_kdr = sec.gbar_kap = 0.0

    def set_axong(self):
        axon = self.axon
        axon.gbar_neykdmc = self.c.gbar_kdmc * self.c.kdmc_gbar_axonm
        axon.gbar_neynax = self.c.gbar_nax * self.c.nax_gbar_axonm
        axon.gbar_neykdr = self.c.gbar_kdr * self.c.kdr_gbar_axonm
        axon.gbar_neykap = self.c.gbar_kap * self.c.kap_gbar_axonm

    def set_calprops(self, sec):
        sec.gcalbar_neycal = self.c.cal_gcalbar
        sec.gcanbar_neycan = self.c.can_gcanbar
        sec.gpeak_neykBK = self.c.kBK_gpeak
        sec.caVhmin_neykBK = -46.08 +self.c.kBK_caVhminShift
        sec.depth_neycadad = self.c.cadad_depth
        sec.taur_neycadad = self.c.cadad_taur

    def set_somag(self):
        sec = self.soma
        sec.gbar_neyih = self.c.gbar_h  # Ih
        self.set_calprops(sec)
        sec.gbar_neykdmc = self.c.gbar_kdmc

    def set_bdendg(self):
        sec = self.Bdend
        sec.gbar_neyih = self.c.gbar_h  # Ih
        self.set_calprops(sec)

    def set_apicg(self):
        h.distance(0, 0.5, sec=self.soma)  # middle of soma is origin for distance
        self.nexusdist = nexusdist = 300.0
        self.h_gbar_tuftm = h_gbar_tuftm = self.c.h_gbar_tuft / self.c.gbar_h
        self.h_lambda = h_lambda = nexusdist / log(h_gbar_tuftm)
        for sec in self.apic:
            self.set_calprops(sec)
            for seg in sec:
                d = h.distance(seg.x, sec=sec)
                if d <= nexusdist:
                    seg.gbar_neyih = self.c.gbar_h * exp(d / h_lambda)
                else:
                    seg.gbar_neyih = self.c.h_gbar_tuft
        self.apic[1].gcalbar_neycal = self.c.cal_gcalbar * self.c.calginc  # middle apical dend gets more iL

    # set properties
    def set_props(self):
        self.set_geom()
        # cm - can differ across locations
        self.axon.cm = self.c.axonCap
        self.soma.cm = self.c.somaCap
        self.Bdend.cm = self.c.bdendCap
        for sec in self.apic: sec.cm = self.c.apicCap
        # g_pas == 1.0/rm - can differ across locations
        self.axon.g_pas = 1.0 / self.c.axonRM
        self.soma.g_pas = 1.0 / self.c.somaRM
        self.Bdend.g_pas = 1.0 / self.c.bdendRM
        for sec in self.apic: sec.g_pas = 1.0 / self.c.apicRM
        for sec in self.all_sec:
            sec.ek = self.c.p_ek  # K+ current reversal potential (mV)
            sec.ena = self.c.p_ena  # Na+ current reversal potential (mV)
            sec.Ra = self.c.rall;
            sec.e_pas = self.c.Vrest  # passive
            sec.gbar_neynax = self.c.gbar_nax  # Na
            sec.gbar_neykdr = self.c.gbar_kdr  # KDR
            sec.vhalfn_neykdr = self.c.kdr_vhalfn  # KDR kinetics
            sec.gbar_neykap = self.c.gbar_kap  # K-A
            sec.vhalfn_neykap = self.c.kap_vhalfn  # K-A kinetics
            sec.vhalfl_neykap = self.c.kap_vhalfl
            sec.tq_neykap = self.c.kap_tq
        self.set_somag()
        self.set_bdendg()
        self.set_apicg()
        self.set_axong()

    def insert_conductances(self):
        for sec in self.all_sec:
            sec.insert('k_ion')
            sec.insert('na_ion')
            sec.insert('pas')  # passive
            sec.insert('neynax')  # Na current
            sec.insert('neykdr')  # K delayed rectifier current
            sec.insert('neykap')  # K-A current
        for sec in [self.Adend3, self.Adend2, self.Adend1, self.Bdend, self.soma]:
            sec.insert('neyih')  # h-current
            sec.insert('ca_ion')  # calcium channels
            sec.insert('neycal')  # cal_mig.mod
            sec.insert('neycan')  # can_mig.mod
            sec.insert('neycadad')  # cadad.mod - calcium decay
            sec.insert('neykBK')  # kBK.mod - ca and v dependent k channel
        for sec in [self.soma, self.axon, self.ais]:
            sec.insert('neykdmc')  # K-D current in soma & axon only
            sec.insert("pas")

#
def prmstr(p, s, fctr=2.0, shift=5.0):
    if p == 0.0:
        print(s, '=', str(p - shift), str(p + shift), str(p), 'True')
    elif p < 0.0:
        print(s, '=', str(p * fctr), str(p / fctr), str(p), 'True')
    else:
        print(s, ' = ', str(p / fctr), str(p * fctr), str(p), 'True')



