# simplified corticospinal cell model (6 compartment) according to Neymotin 2017
from neuron import h
from math import exp, log
import numpy as np
# Goals
# AIS needs to have at least 50 fold higher concentration of Na channels
# AIS needs to have the D-Type Potassium channels Kv1.1 (slowly inactivating potassium current)

# M-Current (Kv7.2/7.3) KCNQ2 and KCNQ3 NICHT MUSKARINERG
# which activates when an excitatory stimulus depolarizes the neuron toward spike threshold, repolarizes the membrane back toward resting potential, and suppresses firing (Rogawski, 2000; Wang et al., 1998). IM is a non-inactivating K+ current that is active at subthreshold membrane potentials and has been implicated in dampening repetitive firing and general neuronal excitability (Brown and Passmore, 2009)




# Initial Set
Vrest = -75.33842520399517 # Calculated from the resting Potential through rmp.caluclateRMP()
h.celsius = 25.0  # for in vitro opt # 25


celcius = 37
Ri = 100
Cm = 0.9
Rm = 150000
v_init= -88
spinescale = 2


na_myelin=80 # in dendrites pS/um2
Km_soma=5 # units in pS/um2
Km_axon=50

Kv_soma=20 #(HVA Kv current, units in pS/um2)
Kv_axon=2000

Kv1_soma=0.01 #0.01 mOhm/cm2 = 100 pS/um2 (LVA-Kv current)
Kv1_axon=0.20 

spinescale=2

vshift_na=10 #provides AP threshold of ~ -60 mV at soma 
vshift_nax=10

# Geometric Properties 

somaL = 48.4123467666
somaDiam = 28.2149102762
axonL = 594.292937602
axonDiam = 1.40966286462

aisL = 31.2 #TODO: Set from experimental data
aisDiam = axonDiam

# Adjusted Axon Length
adjustedAxonL = axonL - aisL
adjustedAxonDiam = axonDiam

# Goals
# AIS needs to have at least 50 fold higher concentration of Na channels
# AIS needs to have the D-Type Potassium channels  (slowly inactivating potassium current)



# Conductance of Na1.6 at the AIS ( N 200 | um² , conductance single 15pS) = 3000pS in Microscopy  (Colbert and Johnston, 1996), but electrophysoliogy measurement N2-4 / um² 


# Ih (hyperpolarization activated cation current) (HCN channel) SUFFIX ih
# Age depended increa in dendritic desnitry (Atkinson and Williams, 2009)

# Na T (transient sodium current) SUFFIX nax

# Na P (persistent sodium current)

# CaLVA (low voltage-activated calcium current) (T-type calcium channel)

# CaHVA (high voltage- activated calcium current) (L, N, P/Q and R types)  SUFFIX cal (LType) SUFFIX (can)

# K A (transient A-type potassium current) SUFFIX kap
#s. A high-density of both transient and sustained KV channels was observed in all apical dendritic compartments. 


# K-D (delayed rectifier potassium current (sustained efflux)  SUFFIX kdmc SUFFIX kdr

# K-M (muscarinic-activated potassium current) 

# BK (big conductance calcium-activated potassium current) (Kole: Kca1.1) SUFFIX KBK
# is the Kca1.1 in koles model acitvated by calcium

# SK (small conductance calcium-activated potassium current.)

# DISTRIBUTIONS

# Dendrite (0 is soma, 1 is tuft)
i_h_distribution_dend = np.linspace(0,1,100)**2 # Exponential till the apical
na_t_distribuition_dend = np.linspace(1,1,100) #Constant
na_t_distribuition_dend = np.linspace(1,1,100) #Constant
ca_lva_distribuition_dend = np.linspace(1,1,100) #Constant
ca_hva_distribuition_dend = np.linspace(1,1,100) #Constant


k_a_distribution_dend = np.linspace(1,0,100) # Slope Down # transient A-type potassium current 
k_d_distribuition_dend = np.linspace(1,0,100) # Slope Down  #KD—delayed rectifier potassium current
bk_distribuition_dend = np.linspace(1,1,100) #Constant 
sk_distribuition_dend = np.linspace(1,1,100) #Constant # is a subfamily of cacium acitvated channels

# Ais (1 is soma, 0 is end of AIS)
i_h_distribution_ais = np.linspace(0,0,100) # Not there
na_t_distribution_ais = np.concatenate((np.linspace(0,1,50),np.linspace(1,0,50))) # Triangle Wave
na_p_distribution_ais = np.concatenate((np.linspace(0,1,50),np.linspace(1,0,50))) # Triangle Wave
ca_lva_distribuition_ais = np.linspace(1,1,100) #Constant (not known)
ca_hva_distribuition_ais = np.linspace(1,1,100) #Constant (not known)

k_a_distribution_ais = np.linspace(1,0,100) # Slope Down 
k_d_distribution_ais = np.concatenate((np.linspace(0,1,50),np.linspace(1,0,50))) # Triangle Wave  
k_m_distribution_ais = np.concatenate((np.linspace(0,1,50),np.linspace(1,0,50))) # Triangle Wave 
bk_distribution_ais = np.linspace(1,1,100) #Constant (not known)
sk_distribution_ais = np.linspace(1,1,100) #Constant (not known)


# KA Section (Harnett et all)





# Adjusted Axon Length
adjustedAxonL = axonL - aisL
adjustedAxonDiam = axonDiam

apicL = 261.904636003
apicDiam = 1.5831889597
bdendL = 299.810775175
bdendDiam = 2.2799248874

# passive properties 
axonCap = 1.01280903702
aisCap = 1.01280903702
somaCap = 1.78829677463
apicCap = 1.03418636866
bdendCap = 1.89771901209

# .Ra Axial resistance (normal rat subthalamic nuceol 123 ohm-cm) #TODO: Literature
rall = 114.510490019


#axonRM = 3945.2107187
#aisRM = 3945.2107187
#somaRM = 18501.7540916
#apicRM = 10751.193413
#bdendRM = 13123.00174

axonRM = aisRM = somaRM = apicRM = bdendRM = 13000

# eSodium reversal potential calculated on T(295.15K) z(1) x{out}(124mM NaCL + 1.25 mM NaH2Po4 + 26 mM NaH2Co3) x{in}(10 mM Na-phosphocreatine + 0.3 mM Na-GTP)
p_ena = 68.16371 #mV
# Potassium reversal potential calculated on T(295.15K) z(1) x{out}(2.5mM KCl) x{in}(4mM KCL + 35K-Glucose) ATTENTION: Buffered with K-Oh to ph7.2
p_ek = -97.55644 #mV



# Na, K reversal potentials calculated from BenS internal/external solutions via Nernst eq.
#p_ek = -104.0  # these reversal potentials for in vitro conditions
#p_ena = 42.0


# Ih Channels Conductance 680 fS, 550 N / m2, expontentially increasing # modeled by a tenfold increased single channel conductance (6.8ps) but normal densitry (kole 2006)

# h-current
h.erev_ih = -47.0  # global
gbar_h = 0.000140956438043 * 20
h_gbar_tuft = 0.00565  # mho/cm^2 (based on Harnett 2015 J Neurosci)



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

    def set_morphology(self):

        self.soma = h.Section(name="soma", cell=self)
        self.axon = h.Section(name="axon", cell=self)

        # TODO: Implement AIS Here
        self.Bdend = h.Section(name="Bdend", cell=self)
        self.Adend1 = h.Section(name="Adend1", cell=self)
        self.Adend2 = h.Section(name="Adend2", cell=self)
        self.Adend3 = h.Section(name="Adend3", cell=self)

        self.ais = h.Section(name="ais", cell=self)


        self.apic = [self.Adend1, self.Adend2, self.Adend3]
        self.basal = [self.Bdend]
        self.alldend = [self.Adend1, self.Adend2, self.Adend3, self.Bdend]
        self.all_sec = [ self.soma, self.ais, self.axon ] + self.apic + self.basal + self.alldend


        # Set Geometry here
        self.set_geom()

        self.ais.connect(self.soma, 0.0, 0.0)
        self.axon.connect(self.ais, 1.0, 0.0)
        self.Bdend.connect(self.soma, 0.5, 0.0)  # soma 0.5 to Bdend 0
        self.Adend1.connect(self.soma, 1.0, 0.0)
        self.Adend2.connect(self.Adend1, 1.0, 0.0)
        self.Adend3.connect(self.Adend2, 1.0, 0.0)

    def set_geom(self):
        self.axon.L = adjustedAxonL
        self.axon.diam = adjustedAxonDiam

        self.ais.L = aisL
        self.ais.diam = aisDiam
        self.soma.L = somaL
        self.soma.diam = somaDiam

        for sec in self.apic: sec.L, sec.diam = apicL, apicDiam

        self.Bdend.L = bdendL
        self.Bdend.diam = bdendDiam

    def activeoff(self):
        for sec in self.all_sec: sec.gbar_nax = sec.gbar_kdr = sec.gbar_kap = 0.0
            
    def set_aisg(self):
        ais = self.ais
        ais.gbar_kdmc = gbar_kdmc * kdmc_gbar_aism
        ais.gbar_nax = gbar_nax * nax_gbar_aism
        ais.gbar_kdr = gbar_kdr * kdr_gbar_aism
        ais.gbar_kap = gbar_kap * kap_gbar_aism

    def set_axong(self):
        axon = self.axon
        axon.gbar_kdmc = gbar_kdmc * kdmc_gbar_axonm
        axon.gbar_nax = gbar_nax * nax_gbar_axonm
        axon.gbar_kdr = gbar_kdr * kdr_gbar_axonm
        axon.gbar_kap = gbar_kap * kap_gbar_axonm

    def set_calprops(self, sec):
        sec.gcalbar_cal = cal_gcalbar
        sec.gcanbar_can = can_gcanbar
        sec.gpeak_kBK = kBK_gpeak
        sec.caVhmin_kBK = -46.08 + kBK_caVhminShift
        sec.depth_cadad = cadad_depth
        sec.taur_cadad = cadad_taur

    def set_somag(self):
        sec = self.soma
        sec.gbar_ih = gbar_h  # Ih
        self.set_calprops(sec)
        sec.gbar_kdmc = gbar_kdmc

    def set_bdendg(self):
        sec = self.Bdend
        sec.gbar_ih = gbar_h  # Ih
        self.set_calprops(sec)

    def set_apicg(self):
        h.distance(0, 0.5, sec=self.soma)  # middle of soma is origin for distance
        self.nexusdist = nexusdist = 300.0
        self.h_gbar_tuftm = h_gbar_tuftm = h_gbar_tuft / gbar_h
        self.h_lambda = h_lambda = nexusdist / log(h_gbar_tuftm)
        for sec in self.apic:
            self.set_calprops(sec)
            for seg in sec:
                d = h.distance(seg.x, sec=sec)
                if d <= nexusdist:
                    seg.gbar_ih = gbar_h * exp(d / h_lambda)
                else:
                    seg.gbar_ih = h_gbar_tuft
        self.apic[1].gcalbar_cal = cal_gcalbar * calginc  # middle apical dend gets more iL

    # set properties
    def set_props(self):
        self.set_geom()
        # cm - can differ across locations
        self.axon.cm = axonCap
        self.ais.cm = aisCap
        self.soma.cm = somaCap
        self.Bdend.cm = bdendCap
        for sec in self.apic: sec.cm = apicCap
        # g_pas == 1.0/rm - can differ across locations
        self.axon.g_pas = 1.0 / axonRM
        self.ais.g_pas = 1.0 / aisRM
        self.soma.g_pas = 1.0 / somaRM
        self.Bdend.g_pas = 1.0 / bdendRM
        for sec in self.apic: sec.g_pas = 1.0 / apicRM
        for sec in self.all_sec:
            sec.ek = p_ek  # K+ current reversal potential (mV)
            sec.ena = p_ena  # Na+ current reversal potential (mV)
            sec.Ra = rall; # Axial resistance
            sec.gbar_nax = gbar_nax  # Na
            sec.gbar_kdr = gbar_kdr  # KDR
            sec.vhalfn_kdr = kdr_vhalfn  # KDR kinetics
            sec.gbar_kap = gbar_kap  # K-A
            sec.vhalfn_kap = kap_vhalfn  # K-A kinetics
            sec.vhalfl_kap = kap_vhalfl
            sec.tq_kap = kap_tq
        self.set_somag()
        self.set_bdendg()
        self.set_apicg()
        self.set_axong()
        self.set_aisg()

    def insert_conductances(self):
        for sec in self.all_sec:
            sec.insert('Kv1')
            sec.insert('Km')
            sec.insert('Kv')  
            sec.insert('pas')  
        for sec in [self.Adend3, self.Adend2, self.Adend1, self.Bdend, self.soma]:
            sec.insert('na')  # h-current
            sec.insert('ca_ion')  # calcium channels
            sec.insert('cal')  # cal_mig.mod
            sec.insert('can')  # can_mig.mod
            sec.insert('cadad')  # cadad.mod - calcium decay
            sec.insert('kBK')  # kBK.mod - ca and v dependent k channel
        for sec in [self.soma, self.axon, self.ais]:
            sec.insert('kdmc')  # K-D current in soma & axon only
            sec.insert("pas")

#
def prmstr(p, s, fctr=2.0, shift=5.0):
    if p == 0.0:
        print(s, '=', str(p - shift), str(p + shift), str(p), 'True')
    elif p < 0.0:
        print(s, '=', str(p * fctr), str(p / fctr), str(p), 'True')
    else:
        print(s, ' = ', str(p / fctr), str(p * fctr), str(p), 'True')



