# -*- coding: utf-8 -*-
'''
@author: stransky
'''

from mrj.io.swc import SwcParser
from mp.model.synapse import Synapse

from mrj.model.experiment import IClamp
from mrj.model.experiment import Experiment
from mrj.model.experiment import Neuron_passive_parameter
from mrj.model.experiment import plot

swc_parser  = SwcParser()
m_swc   = swc_parser.parse("../../../Data/test.swc") #slow swc_parser.parse("../../../Data/MSO_Neurone/P09/ohne_axon_SkeletonTree_0015_2-1G.swc")
m_swc.neuron_create() # initialisiert die Synapse in NEURON
compartment = m_swc.compartments[5]

s = Synapse(compartment, 0.5, syntimes=[2.71])
s.neuron_create()

clamp = IClamp(compartment=compartment, position=0.5, delay=1,amplitude=-10,duration=3)
neuron_passive_parameter = \
    Neuron_passive_parameter(Ra=1, g=0.004, e=-60, nseg=10)

f = Experiment(m_swc, clamp, neuron_passive_parameter, 
               description='passive channels')

f.neuron_create()
r_in, tau_eff = f.run_simulation(duration=8,dt=0.001)
tau_eff_fit     = f.tau_fit()
#TODO: ?why can linfit fail?
import math
if math.isnan(tau_eff_fit):
    tau_eff_fit = 1
print (tau_eff, tau_eff_fit)
#print (r_in, tau_eff, f.tau_fit());
#f.plot_fit(r_in, tau_eff, tau_eff_fit, save=False)
plot( [f] );#[e] works, (e) is flattend
