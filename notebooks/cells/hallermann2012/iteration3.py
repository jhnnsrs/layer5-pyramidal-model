# simplified corticospinal cell model (6 compartment) according to Neymotin 2017

from neuron import h
from math import exp, log
import numpy as np

# Parameters
from cells.baseconfig import BaseConfig


class Config(BaseConfig):
    celsius = 33
    Ri = 80
    Cm = 1.0
    Rm = 15000
    v_init = -85
    Spinescale = 1.5

    # Channel Densities
    # Distribution Arguments
    na_max = 5000
    na_ais=7000
    na_node=7000
    na_collat=500
    na_soma=500
    na_dend=20
    na_myelin=40

    max_distance_apic=300# distance where gbar_na has reached na_dend
    max_distance_basal=100 #corresponding distance for basal dendrites

    vShift_na=10 # affecting activation and inactivation
    vShift_inact_na=10 # affects only inactivation
    q10_na=3
    q10h_na=3

    vShift_nax=10
    vShift_inact_nax=10
    q10_nax=3.0
    q10h_nax=3.0


    length_constant_Kv_and_Kv1=	80 #length constant of the exponential decay in um of Kv and Kv1 according to Keren et al., 2009, JPhysiol, 587:1413-37
    Kv_soma=100
    Kv_dend=3

    Kv1_dend=0.3
    Kv1_ais=2000
    Kv1_soma=100
    Kv1_collat=400

    vShift_Kv1=10
    vShift_inact_Kv1=15

    Kv7_soma=1
    Kv7_dend=1
    Kv7_ais=7

    ca_reduce_fac=0.1

    gca_dend=2.0*ca_reduce_fac
    gca_soma=2.0*ca_reduce_fac
    gca_ais=0.0

    #Ih Related Stuff
    gbar_ih_soma = 0.0002
    gbar_ih_axon = 0.0005
    gh_ih = -0.0002
    eh_ih = -45
    A = 0.0002 #0.000429 before, and also the J Neurosci value I aim for a ~20 Mohm input resistance
    tau = 0.003087

    git2_ais=0.0
    git2_dend=2.0*ca_reduce_fac
    git2_apic=6.0*ca_reduce_fac
    git2_soma=2.0*ca_reduce_fac

    gkca_soma=1.0*ca_reduce_fac
    gkca_dend=1.0*ca_reduce_fac
    gkca_ais=1.0*ca_reduce_fac

    gCa_HVA_apic_hot_fac=1 #i.e. no Ca hot spot//ca
    gCa_LVAst_apic_hot_fac=1 #it2
    gCa_hot_start=685
    gCa_hot_end=885

    spinescale=1.5
    sheaths=10 #number of myelin sheaths

    # geom properties
    somaL = 48.4123467666
    somaDiam = 28.2149102762
    aisL = 20.292937602
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
    p_ek = -98  # these reversal potentials for in vitro conditions
    p_ena = 55
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
        self.distributions["nax_ais"] = np.concatenate(
            (np.linspace(0, self.na_ais, 40), np.linspace(self.na_ais, self.na_ais, 60)))  # Sloping Wave
        self.distributions["na_ais"] = np.concatenate(
            (np.linspace(self.na_soma, self.na_max, 20), np.linspace(self.na_max, 0, 20),
             np.linspace(0, 0,
                         60)))  # Sloping Wave #TODO: Check if it is really wanted to set expression zero

        self.distributions["Kv1_ais"] = np.linspace(self.Kv1_soma, self.Kv1_ais, 100)  # Sloping Wave
        self.distributions["Kv7_ais"] = np.concatenate((np.linspace(self.Kv7_soma, self.Kv7_soma, 40),
                                                        np.linspace(self.Kv7_soma, self.Kv7_ais, 60)))

    def setSynapses(self):
        pass




###############################################################################
# SPI6 Cell
###############################################################################
class SPI6(object):
    "Kole Model in Action"

    def __init__(self,config = Config(), x=0, y=0, z=0, ID=0):
        self.x, self.y, self.z = x, y, z
        self.c = config
        self.ID = ID
        self.all_sec = []
        self.set_morphology()
        self.set_biophysics()
        #self.enhanceCalcium()


    def set_morphology(self):

        self.soma = h.Section(name="soma", cell=self)
        self.axon = h.Section(name="axon", cell=self)
        self.ais = h.Section(name="ais", cell=self)

        # TODO: Implement AIS Here
        self.Bdend = h.Section(name="Bdend", cell=self)
        self.Adend1 = h.Section(name="Adend1", cell=self)
        self.Adend2 = h.Section(name="Adend2", cell=self)
        self.Adend3 = h.Section(name="Adend3", cell=self)


        self.apic = [self.Adend1, self.Adend2, self.Adend3]
        self.basal = [self.Bdend]
        self.alldend = [self.Adend1, self.Adend2, self.Adend3, self.Bdend]
        self.all_sec = [ self.soma, self.axon, self.ais] + self.apic + self.basal + self.alldend 


        # Set Geometry here
        self.set_geom()


        self.ais.connect(self.soma, 0.0, 0.0)
        self.axon.connect(self.ais, 1.0, 0.0)
        self.Bdend.connect(self.soma, 0.5, 0.0)  # soma 0.5 to Bdend 0
        self.Adend1.connect(self.soma, 1.0, 0.0)
        self.Adend2.connect(self.Adend1, 1.0, 0.0)
        self.Adend3.connect(self.Adend2, 1.0, 0.0)


    def set_biophysics(self):
        self.setAllSec()

        # Will Overwrite Default Values
        self.setAIS()
        self.setDendrites()
        self.setSoma()

        return


    def set_geom(self):
        self.axon.L = self.c.axonL
        self.ais.L = self.c.aisL
        self.ais.diam = self.c.aisDiam
        self.axon.diam = self.c.axonDiam
        self.soma.L = self.c.somaL
        self.soma.nseg = 100
        self.soma.diam = self.c.somaDiam

        for sec in self.apic: sec.L, sec.diam, sec.nseg = self.c.apicL, self.c.apicDiam, 100

        self.Bdend.L = self.c.bdendL
        self.Bdend.nseg = 100
        self.Bdend.diam = self.c.bdendDiam


        
    def setDendrites(self):
        # Insert Channels
        for sec in self.alldend:
            sec.insert("na")
            sec.insert("Kv")
            sec.insert("Kv7")
            sec.insert("kca")
            sec.insert("cad")
            sec.insert("ca")
            sec.insert("it2")
            sec.insert("Kv1")
            sec.insert("ih")
            

        g_pas_dend=1/(self.c.Rm/self.c.spinescale)
        cm_dend=self.c.Cm*self.c.spinescale

        # Set Params for Channels
        for sec in self.alldend:
            sec.g_pas = g_pas_dend
            sec.cm = cm_dend
            sec.gbar_Kv=self.c.Kv_dend * self.c.spinescale
            sec.gbar_Kv7=self.c.Kv7_dend
            sec.gbar_kca=self.c.gkca_dend * self.c.spinescale
            sec.gbar_ca = self.c.gca_dend*self.c.spinescale
            sec.gcabar_it2 = self.c.git2_dend*self.c.spinescale
            sec.gbar_Kv1 = self.c.Kv1_dend

        def fromtodistance(origin_segment, to_segment):
            h.distance(0, origin_segment.x, sec=origin_segment.sec)
            return h.distance(to_segment.x, sec=to_segment.sec)

        basaldistribution = {}
        basaldistribution["sec"] = []
        basaldistribution["ih"] = []
        basaldistribution["Kv"] = []
        basaldistribution["Kv1"] = []
        basaldistribution["na"] = []
        basaldistribution["d"] = []
        for sec in self.basal: # middle of soma is origin for distance
            for seg in sec:
                h.distance(0, 0,5, sec=self.soma)
                d = h.distance(seg.x, sec=sec)
                seg.gbar_ih = (self.c.gh_ih+(self.c.A*exp(self.c.tau*d)))*self.c.spinescale
                seg.gbar_Kv = (self.c.Kv_dend+(self.c.Kv_soma-self.c.Kv_dend)*exp(-d/self.c.length_constant_Kv_and_Kv1))*self.c.spinescale
                seg.gbar_Kv1 =(self.c.Kv1_dend+(self.c.Kv1_soma-self.c.Kv1_dend)*exp(-d/self.c.length_constant_Kv_and_Kv1))*self.c.spinescale
                if d >= self.c.max_distance_basal:
                    seg.gbar_na = self.c.na_dend * self.c.spinescale
                else:
                    seg.gbar_na = (self.c.na_dend + (self.c.na_soma-self.c.na_dend)*(1-(d/self.c.max_distance_basal)))*self.c.spinescale

                basaldistribution["sec"].append(seg.x)
                basaldistribution["d"].append(d)
                basaldistribution["ih"].append(seg.gbar_ih)
                basaldistribution["Kv"].append(seg.gbar_Kv)
                basaldistribution["Kv1"].append(seg.gbar_Kv1)
                basaldistribution["na"].append(seg.gbar_na)


        for key, values in basaldistribution.items():
            self.c.distributions[f"{key}_basal"] = np.array(basaldistribution[key])


        for sec in self.apic:
            #Constants
            sec.gbar_Kv7 = self.c.Kv7_soma
            sec.gcabar_it2 = self.c.git2_apic*self.c.spinescale

            h.distance(0, 0.5, sec=self.soma)  # middle of soma is origin for distance
            for seg in sec:
                h.distance(0, 0,5, sec=self.soma)
                d = h.distance(seg.x, sec=sec)
                seg.gbar_ih = (self.c.gh_ih+(self.c.A*exp(self.c.tau*d)))*self.c.spinescale
                seg.gbar_Kv = (self.c.Kv_dend + (self.c.Kv_soma - self.c.Kv_dend) * exp(-d / self.c.length_constant_Kv_and_Kv1)) * self.c.spinescale
                seg.gbar_Kv1 = (self.c.Kv1_dend + (self.c.Kv1_soma - self.c.Kv1_dend) * exp(-d / self.c.length_constant_Kv_and_Kv1)) * self.c.spinescale
                if d >=self.c. max_distance_apic:
                    seg.gbar_na = self.c.na_dend * self.c.spinescale
                else:
                    seg.gbar_na = (self.c.na_dend + (self.c.na_soma-self.c.na_dend)*(1-(d/self.c.max_distance_basal)))*self.c.spinescale
                    
            
    def setSoma(self):
        sec = self.soma
        sec.insert("na")
        sec.insert("Kv")
        sec.insert("Kv7")
        sec.insert("kca")
        sec.insert("cad")
        sec.insert("ca")
        sec.insert("it2")
        sec.insert("Kv1")
        sec.insert("ih")
        
        sec.gbar_na = self.c.na_soma
        sec.gbar_Kv = self.c.Kv_soma
        sec.gbar_ca = self.c.gca_soma
        sec.gbar_kca = self.c.gkca_soma
        sec.gcabar_it2 = self.c.git2_dend
        sec.gbar_Kv7 = self.c.Kv7_soma
        sec.gbar_Kv1 = self.c.Kv1_soma

        #Ih Related
        sec.gbar_ih = self.c.gbar_ih_soma

    def setAxon(self):
        sec = self.axon
        sec.insert("nax")
        sec.insert("Kv1")
        sec.insert("Kv7")
        sec.insert("ih")

        sec.g_pas = 1/(self.c.Rm * (self.c.sheaths + 1))
        sec.cm = self.c.Cm/(self.c.sheaths + 1)

        sec.gbar_nax = self.c.na_myelin
        sec.gbar_Kv7 = self.c.Kv7_soma
        sec.gbar_Kv1 = 20 #TODO: Magic Number


        #Ih Related
        sec.gbar_ih = self.c.gbar_ih_axon
        
    
    def setAIS(self):
        sec = self.ais
        sec.insert("nax")
        sec.insert("na")
        sec.insert("Kv1")
        sec.insert("Kv7")
        sec.insert("cad")
        sec.insert("ca")
        sec.insert("it2")
        sec.insert("kca")

        # Constant Distribution
        sec.gbar_ca = self.c.gca_ais
        sec.gcabar_it2 = self.c.git2_ais
        sec.gbar_kca = self.c.gkca_ais

        # Distribution Arguments
        sec.nseg = 20

        distribution_nax_ais = self.c.distributions["nax_ais"]
        distribution_na_ais = self.c.distributions["na_ais"]

        distribution_Kv1_ais = self.c.distributions["Kv1_ais"]
        distribution_Kv7_ais = self.c.distributions["Kv7_ais"]

        for indexm, seg in enumerate(sec):
            xpercentage = int((seg.x)*100)
            seg.gbar_nax = distribution_nax_ais[xpercentage]
            seg.gbar_na = distribution_na_ais[xpercentage]
            seg.gbar_Kv1 = distribution_Kv1_ais[xpercentage]
            seg.gbar_Kv7 = distribution_Kv7_ais[xpercentage]
        

    def setAllSec(self):
        # Insert Conductances
        for sec in self.all_sec:
            sec.insert('pas')  # passive
            sec.insert("charge_")
            sec.insert("it2") # TODO: No reference for this in Hallermann
        
        # Set Parameters
        for sec in self.all_sec:
            sec.Ra = self.c.Ri
            sec.e_pas = self.c.v_init  # passive
            sec.g_pas = 1.0/self.c.Rm
            sec.cm = self.c.Cm
            sec.ena = self.c.p_ena
            sec.ek = self.c.p_ek

    def enhanceCalcium(self):

        #Calcium Enhancement (Larkum et al 1999)
        for sec in self.all_sec:
            sec.vh1_it2 = 56
            sec.vh2_it2 = 415
            sec.ah_it2 = 30
            sec.v12m_it2 = 45
            sec.v12h_it2 = 65
            sec.am_it2 = 3
            sec.vshift_it2 = 10
            sec.vm1_it2 = 50
            sec.vm2_it2 = 125
            
        for sec in self.all_sec:
            if h.ismembrane('ca_ion', sec=sec):
                sec.eca = 140
                h.ion_style("ca_ion",0,1,0,0,0,sec=sec)
                sec.vshift_ca = 8
                
        for sec in self.all_sec:
            sec.caix_kca  = 0.7 #// Ca-sensitivity of Kca channel
            sec.Ra_kca = 1 #// Activation rate Kca 
            sec.Rb_kca = 2.5 #// inactivation rate Kca, important ! 
            sec.taur_cad = 80 #// Calcium extrusion rate, in ms

        # possibility to implement hot - zone
        for sec in self.apic:
            h.distance(0, 0.5, sec=self.soma)  # middle of soma is origin for distance
            for seg in sec:
                d = h.distance(seg.x, sec=sec)
                if d > self.c.gCa_hot_start and d < self.c.gCa_hot_end:
                    sec.gbar_ca = self.c.gca_dend*self.c.gCa_HVA_apic_hot_fac*self.c.spinescale
                    sec.gcabar_it2 = self.c.git2_apic*self.c.gCa_LVAst_apic_hot_fac*self.c.spinescale
                




