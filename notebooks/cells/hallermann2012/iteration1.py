# simplified corticospinal cell model (6 compartment) according to Neymotin 2017

from neuron import h
from math import exp, log
import numpy as np

# Parameters
celsius = 33
Ri = 80
Cm = 1.0
Rm = 15000     
v_init = -85
Spinescale = 1.5

# Channel Densities
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
axonL = 594.292937602
axonDiam = 1.40966286462
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


###############################################################################
# SPI6 Cell
###############################################################################
class SPI6(object):
    "Simplified Corticospinal Cell Model"

    def __init__(self, x=0, y=0, z=0, ID=0):
        self.x, self.y, self.z = x, y, z
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
        self.all_sec = [ self.soma, self.axon] + self.apic + self.basal + self.alldend


        # Set Geometry here
        self.set_geom()


        self.axon.connect(self.soma, 0.0, 0.0)
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
        self.axon.L = axonL
        self.axon.diam = axonDiam
        self.soma.L = somaL
        self.soma.diam = somaDiam

        for sec in self.apic: sec.L, sec.diam = apicL, apicDiam

        self.Bdend.L = bdendL
        self.Bdend.diam = bdendDiam


        
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
            

        g_pas_dend=1/(Rm/spinescale)
        cm_dend=Cm*spinescale 

        # Set Params for Channels
        for sec in self.alldend:
            sec.g_pas = g_pas_dend
            sec.cm = cm_dend
            sec.gbar_Kv=Kv_dend * spinescale
            sec.gbar_Kv7=Kv7_dend 
            sec.gbar_kca=gkca_dend * spinescale
            sec.gbar_ca = gca_dend*spinescale
            sec.gcabar_it2 = git2_dend*spinescale
            sec.gbar_Kv1 = Kv1_dend

        for sec in self.basal:
            h.distance(0, 0.5, sec=self.soma)  # middle of soma is origin for distance
            for seg in sec:
                d = h.distance(seg.x, sec=sec)
                seg.gbar_ih = (gh_ih+(A*exp(tau*d)))*spinescale
                seg.gbar_Kv = (Kv_dend+(Kv_soma-Kv_dend)*exp(-d/length_constant_Kv_and_Kv1))*spinescale
                seg.gbar_Kv1 =(Kv1_dend+(Kv1_soma-Kv1_dend)*exp(-d/length_constant_Kv_and_Kv1))*spinescale
                if d >= max_distance_basal:
                    sec.gbar_na = na_dend * spinescale
                else:
                    sec.gbar_na = (na_dend + (na_soma-na_dend)*(1-(d/max_distance_basal)))*spinescale

        for sec in self.apic:
            #Constants
            sec.gbar_Kv7 = Kv7_soma
            sec.gcabar_it2 = git2_apic*spinescale

            h.distance(0, 0.5, sec=self.soma)  # middle of soma is origin for distance
            for seg in sec:
                d = h.distance(seg.x, sec=sec)
                seg.gbar_ih = (gh_ih+(A*exp(tau*d)))*spinescale
                seg.gbar_Kv = (Kv_dend + (Kv_soma - Kv_dend) * exp(-d / length_constant_Kv_and_Kv1)) * spinescale
                seg.gbar_Kv1 = (Kv1_dend + (Kv1_soma - Kv1_dend) * exp(-d / length_constant_Kv_and_Kv1)) * spinescale
                if d >= max_distance_apic:
                    sec.gbar_na = na_dend * spinescale
                else:
                    sec.gbar_na = (na_dend + (na_soma-na_dend)*(1-(d/max_distance_basal)))*spinescale
            
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
        
        sec.gbar_na = na_soma
        sec.gbar_Kv = Kv_soma
        sec.gbar_ca = gca_soma
        sec.gbar_kca = gkca_soma
        sec.gcabar_it2 = git2_dend
        sec.gbar_Kv7 = Kv7_soma
        sec.gbar_Kv1 = Kv1_soma

        #Ih Related
        sec.gbar_ih = gbar_ih_soma

    def setAxon(self):
        sec = self.axon
        sec.insert("nax")
        sec.insert("Kv1")
        sec.insert("Kv7")
        sec.insert("ih")

        sec.g_pas = 1/(Rm * (sheaths + 1))
        sec.cm = Cm/(sheaths + 1)

        sec.gbar_nax = na_myelin
        sec.gbar_Kv7 = Kv7_soma
        sec.gbar_Kv1 = 20 #TODO: Magic Number


        #Ih Related
        sec.gbar_ih = gbar_ih_axon
        
    
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
        sec.gbar_ca = gca_ais
        sec.gcabar_it2 = git2_ais
        sec.gbar_kca = gkca_ais

        # Distribution Arguments
        sec.nseg = 20
        na_max = 5000

        distribution_nax_ais = np.concatenate((np.linspace(0, na_ais, 40), np.linspace(na_ais, na_ais, 60)))  # Sloping Wave
        distribution_na_ais = np.concatenate((np.linspace(na_soma, na_max, 20), np.linspace(na_max, 0, 20),
                                              np.linspace(0, 0, 60)))  # Sloping Wave #TODO: Check if it is really wanted to set expression zero

        distribution_Kv1_ais = np.linspace(Kv1_soma, Kv1_ais, 100)  # Sloping Wave
        distribution_Kv7_ais = np.concatenate((np.linspace(Kv7_soma, Kv7_soma, 40), np.linspace(Kv7_soma, Kv7_ais, 60)))  # Sloping Wave

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
            sec.Ra = Ri
            sec.e_pas = v_init  # passive
            sec.g_pas = 1.0/Rm
            sec.cm = Cm
            sec.ena = p_ena
            sec.ek = p_ek

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
                if d > gCa_hot_start and d < gCa_hot_end:
                    sec.gbar_ca = gca_dend*gCa_HVA_apic_hot_fac*spinescale
                    sec.gcabar_it2 = git2_apic*gCa_LVAst_apic_hot_fac*spinescale
                




