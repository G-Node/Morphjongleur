# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import morphjongleur.util.auto_string

class Clamp(object):
    '''
    abstract represention of a clamp in an experiment
    '''

    def __init__(self, compartment, position=0.5):
        '''
        Constructor for clamp.
        '''
        self.compartment = compartment
        self.position = position
#        self.compartment_id = compartment.compartment_id
#        self.morphology = compartment._morphology

    def __str__(self):
        return '%s(compartment=<%s>,position=%f)' % (
            self.__class__.__name__, 
            str(self.compartment), self.position
        );

    def __repr__(self):
        return '%s(compartment=<%s>,position=%f)' % (
            self.__class__.__name__, 
            self.compartment.__repr__(), self.position
        );

class VClamp(Clamp):
    '''
    Constructor for VClamp
    represents VClamp in an Experiment
    Implementing a voltage clamp electrode
    '''

    def __init__(self, compartment, position=0.5):#
        '''
        Constructor for VClamp.
        '''
        super(VClamp, self).__init__(compartment, position)
#        self.delay = delay
#        self.amplitude = amplitude
#        self.duration = duration

    def neuron_create(self):
        import neuron
        # Locate the electrode at the given position in given of the compartment
        target = self.compartment.neuron_h_Section
        self.neuron_h_VClamp = neuron.h.VClamp( (target)( self.position ))

        self.neuron_clamp   = self.neuron_h_VClamp


class IClamp(Clamp):
    '''
    represents IClamp in an Experiment
    Implementing a current clamp electrode
    '''
    def __init__(self, compartment, position=0.5, 
                 amplitude=-1e-9, delay=0e-3, duration=3e-3):
        '''
        Constructor for IClamp.
        
        Default values:
        amplitude=-1e-9 [A]
        delay    = 0e-3 [s]
        duration = 3e-3 [s]
        '''
        import sys
        if amplitude > 1e-7 or amplitude < -1e-7:
            print >> sys.stderr, "current unit is [A], no [nA]"
        if delay > 1e-1 :
            print >> sys.stderr, "delay unit is [s], not [ms]: "+str(delay)
        if duration > 1e-2:
            print >> sys.stderr, "duration unit is [s], not [ms]"+str(duration)
        super(IClamp, self).__init__(compartment, position)
        self.amplitude = amplitude
        self.delay = delay
        self.duration = duration

    def neuron_create(self):
        '''
        neuron units:
        amp   [nA]
        delay [ms]
        dur   [ms]
        '''
        import neuron

        # Locate the electrode at the given position in given of the compartment
        target = self.compartment.neuron_h_Section
        self.neuron_h_IClamp = neuron.h.IClamp( (target)( self.position ))
        #stim = neuron.h.IClamp(soma(0.5))
        
        # Setting recording paradigm
        self.neuron_h_IClamp.amp   = self.amplitude * 1e9
        self.neuron_h_IClamp.delay = self.delay     * 1e3
        self.neuron_h_IClamp.dur   = self.duration  * 1e3

        self.neuron_clamp   = self.neuron_h_IClamp

    def __str__(self):
        return '%s(compartment=%i,position=%f, amplitude=%g A,delay=%g s,duration=%g s)' % (
            self.__class__.__name__, 
            self.compartment.compartment_id, self.position,
            self.amplitude, self.delay, self.duration
        );

    def __repr__(self):
        return '%s(compartment=<%s>,position=%f, amplitude=%g A,delay=%g s,duration=%g s)' % (
            self.__class__.__name__, 
            self.compartment.__repr__(), self.position,
            self.amplitude, self.delay, self.duration
        );

class PatternClamp(IClamp):
        
    import numpy
    '''
    represents IClamp in an Experiment
    Implementing a current clamp electrode
    '''
    def __init__(self, compartment, position=0.5, 
                 amplitude=-1e-9, function=numpy.sin, delta_t=1e-3,
                 delays=[], durations=[], default_duration=3e-3):
        '''
        Constructor for SinusClamp.
        
        Default values:
        delay = 0.001 [s]
        amp   = -1e-9 [A]
        function    = function t -> y (should be in [-1,1])
        dur   = 0.003 [s]
        delta_t  = 0.001 [s]
        delays   = [[s]]
        durations  = [[s]]
        default_durations  = [s]
        '''
        self.iclamps    = []
        self.compartment = compartment
        self.position = position
        #self.compartment_id = compartment.compartment_id
        self.amplitude = amplitude
        self.function = function
        self.delta_t = delta_t
        self.delays = delays
        self.delays   = map(float, delays) #numpy type not orm mappable: not list(delays) & = [] leads to shadowing 
        self.durations = map(float, durations) 
#TODO: why defined twice? 
        import numpy    
        for i in range(len(delays)):
            if len(durations) <= i:
                self.durations.append(default_duration)
            for t in numpy.arange(self.delays[i], self.durations[i], self.delta_t):
                self.iclamps.append( 
                    IClamp(
                        compartment=self.compartment, position=self.position, 
                        amplitude=amplitude*self.function(t), 
                        delay=t, duration=delta_t
                    )
                )

    def neuron_create(self):
        self.neuron_clamp   = []
        for iclamp in self.iclamps:
            iclamp.neuron_create()
            self.neuron_clamp.append(iclamp.neuron_clamp)

    def __str__(self):
        for c in self.iclamps:
            print c
        return '%s(%s)' % (
            self.__class__.__name__, 
            str(self.iclamps)
        );

    def __repr__(self):
        return '%s(%s)' % (
            self.__class__.__name__, 
            map(self.iclamps,repr())
        );

@morphjongleur.util.auto_string.auto_string
class Clamp_groups(object):
    pass
