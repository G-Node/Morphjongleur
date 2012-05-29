# -*- coding: utf-8 -*-
'''
@author: stransky
'''

class Neuron_passive_parameter(object):

    def __init__(self, Ra=200, g=0.004, e=-60, neuron_passive_parameter_key=None):
        '''
        Ra Axial Resistivity (Ra) [Ohm * cm]
        g  Passive conductivity (pas.g) [S/cm^2]
        e  Passive reversal potential (pas.e) [mV]
        '''
#TODO: unit change#TODO: rho =80, g = 3 
        self.neuron_passive_parameter_key   = neuron_passive_parameter_key
        self.Ra         = Ra
        self.g          = g
        self.e          = e

    def neuron_create(self, sec):
        '''
        Ra Axial Resistivity (Ra) [Ohm * cm]
        g  Passive conductivity (pas.g) [S/cm^2]
        e  Passive reversal potential (pas.e) [mV]
        '''
        sec.Ra = self.Ra;       # * 1e2
        sec.insert("pas")
        for seg in sec:
            seg.pas.g = self.g; # * 1e2 ??
            seg.pas.e = self.e; # * 1e3

    def __repr__(self):
            return 'Neuron_passive_parameter(Ra=%f, g=%f, e=%f, neuron_passive_parameter_key=%s)' % (
                self.Ra,
                self.g,
                self.e,
                str(self.neuron_passive_parameter_key)
            )

    def __str__(self):
        return """<Neuron_passive_parameter(Ra=%f, g=%f [S/cm^2], e=%f, neuron_passive_parameter_key = %s)>""" % (
                self.Ra, 
                self.g, 
                self.e, 
                str(self.neuron_passive_parameter_key)
        )
