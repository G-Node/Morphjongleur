# -*- coding: utf-8 -*-
'''
@author: stransky
@see http://www.neuron.yale.edu/neuron/static/docs/help/neuron/neuron/mech.html#AlphaSynapse
@see http://www.neuron.yale.edu/neuron/faq#playevents 
@see http://www.neuron.yale.edu/phpbb/viewtopic.php?f=28&t=2117
'''
import numpy

class Synapse(object):

    def __init__(self, compartment, position=0.5, e = -70, gmax = 0.1, tau = 0.1, syntimes=[]):#, mso_experiment_key=None, neuron_key=None
        '''
        Constructor for Alpha-Synapse.
        e
        g
        tau [ms]
        synimes [s]
        '''
#        self.morphology         = compartment._morphology
        self.compartment= compartment
        self.position = position
        self.neuron_h_AlphaSynapse = None
        self.e          = e
        self.gmax       = gmax
        self.tau        = tau
        self.syntimes   = map(float, syntimes) #numpy type not orm mappable: not list() & = [] leads to shadowing 

        if not vars(compartment).has_key('synapse') or compartment.synapse == None:
            compartment.synapse  = []
        compartment.synapse.append( self )#if not stored -> not in NEURON!

    def neuron_create(self):
        '''
        Alpha-Synapse
        maybe Synapse & IClamp should inherit from point process
        '''
        import neuron
        #=======================================================================
        #TODO: http://web.mit.edu/neuron_v7.1/doc/help/neuron/neuron/classes/netcon.html
        # source = neuron.h.Section()
        # target = neuron.h.Section()#self.compartment.neuron_h_Section
        # syn = neuron.h.ExpSyn(target(0.5))
        # nc = neuron.h.NetCon(source(0.5)._ref_v, syn)
        # syntimes = [1.0,2.0,3.0]
        # def loadqueue():
        #    for syntime in syntimes:
        #        nc.event(syntime)
        # 
        # fih = neuron.h.FInitializeHandler("loadqueue()")
        #=======================================================================
        if self.neuron_h_AlphaSynapse == None:
            self.neuron_h_AlphaSynapse  = []
        for syntime in self.syntimes:
                neuron_h_AlphaSynapse = \
                    neuron.h.AlphaSynapse(self.compartment.neuron_h_Section(self.position))
                neuron_h_AlphaSynapse.e       = self.e
                neuron_h_AlphaSynapse.onset   = syntime * 1000 # ms
                neuron_h_AlphaSynapse.gmax    = self.gmax
                neuron_h_AlphaSynapse.tau     = self.tau
                self.neuron_h_AlphaSynapse.append( neuron_h_AlphaSynapse )

    def remove(self):
        '''
        NEURON removes synapses, when their python reference is deleted !!!
        '''
        self.compartment.synapse = None

    def __str__(self):
        return '<%s(compartment=%i,position=%f,e=%f,gmax=%f,tau=%f,syntimes=%s)>' % (
            self.__class__.__name__, 
            self.compartment.compartment_id, self.position, 
            self.e, self.gmax, self.tau, str(self.syntimes)
        )

    def __repr__(self):
        return '%s(compartment=%s,position=%f,e=%f,gmax=%f,tau=%f,syntimes=%s))' % (
            self.__class__.__name__, 
            self.compartment.__repr__(), self.position,
            self.e, self.gmax, self.tau, str(self.syntimes)
        )

class Synapse_MSO(Synapse):
    def __init__(self, compartment, position=0.5, e=-70, gmax=0.1, tau=0.1, 
                 f=500, r=1, fireing_rate=1, duration=50e-3, delay=0e-3):#, mso_neuron_key=None
        '''
        f    carrier frequency Hz = 1/s
        duration [s]
        fireing_rate [0,1]
        delay [s]
        '''
#        self.mso_neuron_key = mso_neuron_key
        self.f  = f #soll
        self.r  = r
        self.fireing_rate   = fireing_rate
        self.duration       = duration
        self.delta_time     = delay

        super(Synapse_MSO, self).__init__(compartment, position=0.5, 
            e = e, gmax = gmax, tau = tau,
            syntimes=getSynapticActivity(f, r=r, fireing_rate=fireing_rate, duration=duration, delay=delay)
        )

        self.vector_strength= vector_strength(f, self.syntimes) #ist

    def similar(self, compartment=None):
        """
        similar properties, different syntimes
        """
        if compartment == None:
            compartment = self.compartment
        return Synapse_MSO(compartment,self.position,self.e,self.gmax,self.tau, self.f, self.r, self.fireing_rate, self.duration, self.delta_time) 

    def __str__(self):
        return '<%s(compartment=%i,position=%f,e=%f,gmax=%f,tau=%f,f=%f,r=%f,fireing_rate=%f,duration=%f,delta_time=%f,vector_strength=%f,%s)>' % (
            self.__class__.__name__, 
            self.compartment.compartment_id, self.position, 
            self.e, self.gmax, self.tau, 
            self.f, self.r, self.fireing_rate, self.duration, self.delta_time, self.vector_strength,
            str(self.syntimes)
        )

    def __repr__(self):
        return '%s(compartment=%s,position=%f,e=%f,gmax=%f,tau=%f,f=%f,r=%f,fireing_rate=%f,duration=%f,delta_time=%f)' % (
            self.__class__.__name__, 
            self.compartment.__repr__(), self.position,
            self.e, self.gmax, self.tau, 
            self.f, self.r, self.fireing_rate, self.duration, self.delta_time
        )

def getSynapticActivity(f=500, r=1, fireing_rate=1, duration=50e-3, delay=0e-3):
    return _getSynapticActivity(f, getSigma(r), fireing_rate, duration, delay)

def _getSynapticActivity(f=500, sigma=0, fireing_rate=1, duration=50e-3, delay=0e-3):
    '''
    f    carrier frequency Hz = 1/s
    duration [s]
    fireing_rate [0,1]
    delay [s]
    return fireing times [s]
    '''
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

def getSigma(r=0.5):
    '''
    r    vector_strength [0,1]
    '''
    assert 0 <= r and r <= 1
    if r == 0:
        return numpy.Infinity
    if r == 1:
        return 0
    #import scipy.special
    #return - 0.5 * scipy.special.erf( r - 1.235176 ) + 0.5
    #sigma = 0.37343167819938 * numpy.log (1.0/r -1 ) + 1.212942707929927    # logistic
    #sigma = numpy.exp( numpy.log( 1.0/r - 1 ) / 6.918950 ) / 0.377030 - 1.475555 # custom
    sigma = ( ( 0.193228 / numpy.log( 2.016684 / r - 1) - 0.002929) / 0.226976 )**(1/-1.690475) # customost
    if sigma > 0:   return sigma
    else:           return 0

def vector_strength(f, spikes):
    '''
    abweichung von sum euklid
    '''
    n = len(spikes)
    sum_x = 0
    sum_y = 0
    for i in range( n ):
        theta = numpy.mod(2 * numpy.pi * f * spikes[i], 2*numpy.pi) # f = 1/T -> omega = delta_t * f
        x   = numpy.cos(theta)
        y   = numpy.sin(theta)
        sum_x += x
        sum_y += y
    r = numpy.sqrt( sum_x**2 +sum_y**2 )/n
    #sum_x , um_y -> theta
    return r

def plot_spike_tains(f, duration, r, fireing_rate, n, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
    import matplotlib.pyplot
    matplotlib.pyplot.xlabel('time [ms]')
    matplotlib.pyplot.ylabel('repetition#')
    for i in range(n):
        spike_train = getSynapticActivity(f=f, r=r, fireing_rate=fireing_rate, duration=duration)
        spike_train *= 1000 #ms
        x = numpy.zeros( len(spike_train) ) + i
        matplotlib.pyplot.scatter(spike_train, x, marker='+')
    matplotlib.pyplot.axis([-0.01*duration*1000, 1.01*duration*1000, -0.01*n, 1.01*n]);
    matplotlib.pyplot.legend(["carrier frequency = %i Hz\nduration               = %i ms\nvector strength    = %.2f\nfireing rate           = %.2f" % (f,duration*1000,r,fireing_rate) ], loc='best');
    matplotlib.pyplot.grid()
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format)
    matplotlib.pyplot.show()
    matplotlib.pyplot.close()

def plot_r(f=500, d=100e-3, dr=0.01, picture_file=None, picture_formats=['png', 'pdf', 'svg']):#x_axis='r',
    """
    ds > 1/f !
    """ 
    import matplotlib.pyplot
    i = 0
    rs = []
    sigmas = []
    ys = []
    print "r_soll ->\tsigma ->\tr"
    datas   = []
    for r in numpy.arange(0, 1+dr, dr) :
        for t in [0] :
            print "%f\t" %(r),
            sigma = getSigma(r)
            print "%f\t" % (sigma),
            rs.append(r)
            sigmas.append(sigma)
            v = getSynapticActivity(f=f, r=r, fireing_rate=1, duration=d, delay=t)
            #print v
            #matplotlib.pyplot.scatter(v, numpy.zeros( len(v) ) + i  )
            r = vector_strength(f, v)
            print "%f" % (r)
            ys.append(r)
            i = i+1
            datas.append([sigma,r])
    numpy.savetxt("../../../Data/%.1f_%f@%i.dat" % (getSigma(dr),dr,int(f*d)), datas) 

    matplotlib.pyplot.figure()
    matplotlib.pyplot.xlabel('sigma')
    matplotlib.pyplot.ylabel('measured vector strength')
    matplotlib.pyplot.xlim(0, getSigma(dr))
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.grid()
    matplotlib.pyplot.scatter(sigmas,ys, marker='x', color='black')#, basex=10, basey=10, ls="-"
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'sigma_'+str(getSigma(dr))+'_'+str(int(f*d))+'.'+picture_format,format=picture_format)
    else:
        matplotlib.pyplot.show()

    matplotlib.pyplot.figure()
    matplotlib.pyplot.xlabel('aimed vector strength')
    matplotlib.pyplot.ylabel('measured vector strength')
    #matplotlib.pyplot.legend(["based on %i examples / dot" % (f*d) ], loc='best');
    matplotlib.pyplot.xlim(0, 1)
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.grid()

    matplotlib.pyplot.scatter(rs,ys, marker='x', color='black')
    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'_'+str(dr)+'_'+str(int(f*d))+'.'+picture_format,format=picture_format)
    else:
        matplotlib.pyplot.show()

    matplotlib.pyplot.close('all')
    datas   = numpy.ndarray((len(datas),2), buffer=numpy.array(datas),dtype=float)
    return datas

if __name__ == "__main__":
    #plot_spike_tains(f=500, duration=35e-3, r=0.64, fireing_rate=1, n=50, picture_file='../../../doc/spike_train_500_64_1_50')
    #plot_spike_tains(f=500, duration=35e-3, r=0.93, fireing_rate=1, n=50, picture_file='../../../doc/spike_train_500_93_1_50')
    plot_r(dr=0.001, picture_file='../../../doc/r_lp_')
    #plot_r(f=1000, d=100.0, dr=0.0001, picture_file='../../../doc/r_lp_')#_10000_1e-3')Data/0-5_100000.dat'
