'''
Created on 30.03.2012

@author: stransky
'''
import numpy
import morphjongleur.model.clamp
import morphjongleur.model.synapse


class SinusClamp(morphjongleur.model.clamp.PatternClamp):
    '''
    represents IClamp in an Experiment
    Implementing a current clamp electrode
    '''
    def __init__(self, compartment, position=0.5, 
                 amplitude=-1e-9, frequency=265, delta_t=1e-4,
                 delay=0e-3, duration=3e-3):
        '''
        Constructor for SinusClamp.
        
        Default values:
        delay = 0.001 [s]
        amp   = -1 [A]
        dur   = 0.003 [s]
        frequency   = 265 [Hz=1/s]
        delta_t  = 0.001 [s]
        '''
        self.frequency = frequency
        self.amplitude = amplitude
        super(SinusClamp, self).__init__(
            compartment=compartment, position=position,
            amplitude=amplitude, function=self.signal, delta_t=delta_t,
            delays=[delay], default_duration=duration
        )

    def signal(self,t):
        import numpy
        return numpy.sin(self.frequency * t * 2*numpy.pi )

    def __str__(self):
        return '%s(compartment=%s,position=%f, amplitude=%g A,frequency=%f Hz,delta_t=%g s)' % (
            self.__class__.__name__, 
            self.compartment, self.position,
            self.amplitude, self.frequency, self.delta_t
        );

    def __repr__(self):
        return '%s(compartment=<%s>,position=%f, amplitude=%g,frequency=%f,delta_t=%g)' % (
            self.__class__.__name__, 
            self.compartment.__repr__(), self.position,
            self.amplitude, self.frequency, self.delta_t
        );


class Pattern_generators(object):
    '''
    generates "airborne vibration" measured at the antennae
    '''


    def __init__(self, frequency=265, delay=0e-3, duration=20e-3, amplitude=1):
        '''
        frequency[Hz=1/s]
        delay[s]
        duration[s]
        amplitude[A oder dB]
        '''
        self.frequency  = frequency
        self.delay      = delay
        self.duration   = duration
        self.amplitude  = amplitude

    def iclamp(self, compartment, position=0.5):
        return morphjongleur.model.clamp.IClamp(
                compartment=compartment, position=position, 
                delay=self.delay, amplitude=self.amplitude, duration=self.duration
            )

    def synapses(self, compartment, position=0.5):
        return morphjongleur.model.synapse.Synapse_MSO(
                compartment, position=0.5, 
                e = -70, gmax = 0.1, tau = 0.1, syntimes=_times(f=self.frequency, duration=self.duration, delay=self.delay)
                
            )
    
def _times(f=265, sigma=0, fireing_rate=1, duration=20e-3, delay=0e-3):
    full_circle =  2*numpy.pi # 1.0 #
    n   = int(numpy.ceil(duration * f))
    #mu  = full_circle/ f  # winkelgeschwindigkeit
    if sigma == 0 :
        return numpy.arange(delay, duration, 1.0/f)    # numpy.array( range(delay, n)  ) / float(f) + delay

    #print "%f, %f, %f"%(mu, sigma, n)
    if numpy.isinf(sigma):
        thetas = numpy.random.uniform(0, full_circle, n)
    else: 
        thetas = numpy.random.normal(0, sigma, n) #http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.normal.html#numpy.random.normal
    thetas = numpy.mod(thetas, full_circle)
    for i in range( len(thetas) ):
        if numpy.random.uniform() < fireing_rate:
            if thetas[i] > full_circle/2.0:    #  |\/ -> /|\ 
                thetas[i]  -= full_circle #TODO? and i > 0 to but prevent < 0
            thetas[i]   = (thetas[i] / full_circle + i) / f + delay 
    return thetas   #syntimes
