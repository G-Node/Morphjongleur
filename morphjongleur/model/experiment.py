# -*- coding: utf-8 -*-
'''
@author: stransky
'''
#import numpy;#scipy big -> numpy    http://www.scipy.org/Getting_Started
#import scipy
#import math;
#from morphjongleur.model.pydesignlib import Neuron_MSO, Simulation
import morphjongleur.model.neuron_passive

class Experiment(object):
    '''
    Experiment using Morphjokey
    '''

    def __init__(self, morphology, recordingpoints, clamps, 
        neuron_passive_parameter=morphjongleur.model.neuron_passive.Neuron_passive_parameter(),
        duration=5e-3, dt=1e-6, nseg=10,
        description=''):#experiment_key=None,
        '''
        morphology = Morphology
        clamp = IClamp
        experiment duration [s]
        dt    temporal resolution [s]
        '''
        self.morphology = morphology
        self.recordingpoints     = []
        try:
            self.recordingpoints.extend(recordingpoints)
        except TypeError:
            self.recordingpoints.append(recordingpoints)
        for recordingpoint in self.recordingpoints:
            recordingpoint.experiment  = self
        self.clamps     = []
        self.iclamps    = []
        self.vclamps    = []
        try:
            self.clamps.extend(clamps)
        except TypeError:
            self.clamps.append(clamps)
        for clamp in self.clamps:
            if isinstance(clamp, morphjongleur.model.clamp.IClamp):
                self.iclamps.append(clamp)
            elif isinstance(clamp, morphjongleur.model.clamp.VClamp):
                self.iclamps.append(clamp)
            else: #TODO: Patternclamp
                raise Exception("unknown clamp")
            clamp.experiment  = self
        self.nseg       = nseg  #TODO: NEURON !
        self.neuron_passive_parameter = neuron_passive_parameter
        self.duration   = duration
        self.dt         = dt
        self.description= description

    def __repr__(self):
        # Optional: for str(autoobj) or print autoobj    #has %s in %self.clamp, self.morphology
        return 'Experiment(morphology=%s,recordingpoints=%s,clamps=%s, neuron_passive_parameter=%s, experiment_key=%s, description=%s)' % (
                self.morphology.__repr__(), self.recordingpoints.__repr__(), self.clamps.__repr__(), 
                self.neuron_passive_parameter.__repr__(),
                str(self.experiment_key), str(self.description) 
        );
    def __str__(self):
        return """<Experiment(
 morphology     = %s
 recordingpoints         = %s
 clamps         = %s
 neuron_passive_parameter = %s 
 experiment_key = %s 
 description    = %s
%s)>""" % (
                str(self.morphology), str(self.recordingpoints), str(self.clamps), 
                str(self.neuron_passive_parameter),
                str(self.experiment_key if vars().has_key('experiment_key') else ''), str(self.description),
                str( str(self._results[0]) if False else '')#self.__dict__.has_key('_results') 
        );
        
    def neuron_create(self):
        import neuron
        """
        set passive properties
        Ra Axial Resistivity (Ra): 200 [Ohm * cm]
        g  Passive conductivity (pas.g) = 0.003 [S/cm^2]
        e  Passive reversal potential (pas.e) = -60 [mV]

        duration experiment duration [ms]
        dt  temporal resolution [ms]
        """
        # Implementing a current clamp electrode
        # ======================================
#TODO: unit change
        self.morphology.neuron_create()
#TODO: ??? 
#        self.cell = Neuron_MSO();
#        self.sim = Simulation(self.cell);
#        self.sim.set_IClamp();
#        self.cell.set_passive_parameters(gp=0.004, E=-60, rho=1);
        
        # An important part of NEURON-design is the following:
        # In order to increase the resolution of `neuron.h.Section()`s (=compartments),
        # each neuron_h_Section is divided into `nseg` segments of equal length.
        # Some parameters are provided for a entire neuron_h_Section, others for a segment:
        #TODO: nseg!
        
        #self.E  = E;#TODO wird E immer gesetzt?
        #if(use_passive_channes):
        #    self.set_passive_parameters(gp, E, Ra);

        for sec in neuron.h.allsec():
            self.neuron_passive_parameter.neuron_create(sec)
        
        for clamp in self.clamps:
            clamp.neuron_create()
            self.stim   = clamp.neuron_clamp

        for recordingpoint in self.recordingpoints:
            recordingpoint.neuron_create()

        neuron.h.dt = self.dt   * 1e3
        neuron.h.finitialize(self.neuron_passive_parameter.e);
        self.neuron_duration    = self.duration * 1e3

        neuron.init();

    def run_simulation(self):
        import neuron
        self.neuron_create()
        neuron.run(self.neuron_duration);

class RecordingPoint(object):#Clamp
    '''
    represents a RecordingPoint in an Experiment
    '''

    def __init__(self, compartment, position=0.5, experiment=None):
        '''
        Constructor for RecordingPoint.
        '''
        self.compartment= compartment
        self.position   = position
        self.experiment = experiment

    def neuron_create(self):
        import neuron
        # Record Time
        self.rec_t = neuron.h.Vector();
        self.rec_t.record(neuron.h._ref_t);
        # Record Voltage
        self.rec_v = neuron.h.Vector();
        self.rec_v.record( ( self.compartment.neuron_h_Section )( self.position )._ref_v);

    def get_voltage_trace(self, delay=0):
        if(len(self.rec_t) == 0):
            raise Exception("run experiment first");
        return VoltageTrace(t=self.rec_t,v=self.rec_v, delay=delay, recordingpoint=self)

    def get_tau_fit(self, iclamp):
        if(len(self.rec_t) == 0):
            raise Exception("run experiment first");
        self.tau_fit   = TauFit(t=self.rec_t,v=self.rec_v, iclamp=iclamp, recordingpoint=self)
        return self.tau_fit



class Experiment_info(object):
    pass

class Experiment_groups(object):
    pass

#class Result(object):
#    pass

class VoltageTrace(object):
    def __init__(self, t,v, delay=0, recordingpoint=None):
        assert len(t) == len(v)
        self.recordingpoint = recordingpoint
        self.delay  = delay
        self.t      = list(t)
        self.v      = list(v)

        if len(self.t) != len(self.v):
            import sys
            print >> sys.stderr, "WARNING: len(self.rec_t) %i != %i len(self.rec_v)" %(len(self.rec_t),len(self.rec_v))

        self.v_min = float('+inf')
        self.v_max = float('-inf')
        for i in xrange(len(self.t)):
            self.t[i]   /= 1e3 # ms -> s
            if(self.t[i] > delay):
                if(self.v[i] < self.v_min):
                    self.t_min  = self.t[i]
                    self.v_min  = self.v[i]
                if(self.v[i] > self.v_max):
                    self.t_max  = self.t[i]
                    self.v_max  = self.v[i]
        
        self.amplitude = (self.v_max - self.v_min)/2.
        self.duration  = abs(self.t_max - self.t_min)

    def plot(self, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import matplotlib.pyplot
        matplotlib.pyplot.xlabel("t [ms]")
        matplotlib.pyplot.ylabel("U [mV]")
        #matplotlib.pyplot.xlim(0, self.recordingpoint.experiment.duration*1000)
        #matplotlib.pyplot.ylim(min(self.recordingpoint.experiment.neuron_mso.e, self.recordingpoint.experiment.e), max(self.recordingpoint.experiment.neuron_mso.e, self.recordingpoint.experiment.e))
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.plot(self.t, self.v)

        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=True)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

    def plot_ensemble(self, voltage_traces=[], title=None, xlim=None, ylim=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        """
        needs a list of voltage traces
        """
        import matplotlib.pyplot
        #figure = matplotlib.pyplot.figure()
        titles  = []

        for vt in voltage_traces:
            matplotlib.pyplot.plot(vt.t, vt.v, label="%s @ %i" % (vt.recordingpoint.experiment.description, self.recordingpoint.compartment.compartment_id) )
            titles.append(vt.recordingpoint.experiment.description)
#TODO:  matplotlib.pyplot.axis('image');
#        xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis();
#        matplotlib.pyplot.axis([int(xmin), int(xmax + 0.5), int(ymin - 1.01), int(ymax + 0.01)]);#int(-0.5) = 0
    
        matplotlib.pyplot.ylabel('Voltage [mV]')
        if ylim != None:
            matplotlib.pyplot.ylim(ylim)
        matplotlib.pyplot.xlabel('Time [ms]')
        if xlim != None:
            matplotlib.pyplot.xlim(xlim)
        matplotlib.pyplot.legend(titles, loc='best');
        if title    == None:
            title   = self.recordingpoint.experiment.description
        matplotlib.pyplot.title(title);

        matplotlib.pyplot.grid();

        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=True)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()


    def __str__(self):
        return '<%s(%s: experiment_key=%s, t_min=%f [ms],v_min=%f [mV], t_max=%f [ms],v_max=%f [mV], len(t)=%i,len(v)=%i)>' % (
            self.__class__.__name__, str(self.result_key if vars().has_key('result_key') else ''), 
            str(self.experiment_key if vars().has_key('experiment_key') else ''),
            self.t_min,self.v_min,
            self.t_max,self.v_max,
            len(self.t),len(self.v)
        )

    def __repr__(self):
        return '%s(t_min=%f [ms],v_min=%f, t_max=%f,v_max=%f, t=%s,v=%s, result_key=%s,morphology_key=%s)' % (
            self.__class__.__name__, 
            self.t_min,self.v_min,
            self.t_max,self.v_max,
            str(self.t),str(self.v), 
            str(self.result_key if vars().has_key('result_key') else ''), str(self.experiment_key if vars().has_key('experiment_key') else ''),
        )

class TauFit(object):
    
    def __init__(self, t,v, iclamp, recordingpoint=None):#r_in, tau_eff, tau_eff_fit=0,
        assert len(t) == len(v)
        self.recordingpoint = recordingpoint
        self.iclamp     = iclamp   
        self.t              = t
        self.v              = v
        self.r_in       = float(self.R_in())
        self.tau_eff    = float(self.tau_eff())
        self.tau_eff_fit= float(self.tau_fit())
        self.g_max          = float('nan')#TODO:g_max

    def get_recording(self, start=0, end=float('+Infinity')):
        """
        return time, voltage arrays
        """
        import numpy
        time = numpy.array(self.t)
        voltage = numpy.array(self.v)
        start_index = int(start/self.recordingpoint.experiment.dt);#floor
        if(end==float('+Infinity')):
            end_index = None;
        else:
            end_index   = int(numpy.ceil(end/self.recordingpoint.experiment.dt)) + 1;
        time    = time[start_index:end_index];
        voltage = voltage[start_index:end_index];
        assert len(time) == len(voltage)
        return time, voltage;

    def get_tau_eff(self):
        """
        original version
        """
        import numpy
        time, voltage = self.get_recording();
        vsa = numpy.abs(voltage-voltage[0]); #vsa: voltage shifted and absolut
        v_max = numpy.max(vsa);
        exp_val = (1-1/numpy.e) * v_max; # 0.6321 * v_max
        ix_tau = numpy.where(vsa > ( exp_val ));
        if len(ix_tau[0]) == 0:
            import sys
            print >> sys.stderr, "IClamp seems to be used as sensor. use RecordingPoint instead -> get_tau_eff := 0"
            return 0
        tau = time[ix_tau[0][0]] - self.iclamp.delay;
        return tau;

    def tau_eff(self):
        """
        tau_eff (effective time constant): time after stimulus onset when the voltage trace reaches 1 - 1/e of its maximal change (~63%).
        """
        import numpy
        time, voltage = self.get_recording(self.iclamp.delay,self.iclamp.delay+self.iclamp.duration)
        end   = self.v[int(numpy.ceil((self.iclamp.delay+self.iclamp.duration)/self.recordingpoint.experiment.dt))];
        start = self.v[int(self.iclamp.delay/self.recordingpoint.experiment.dt)];#floor
        u_eff = (end - start)*(1 - numpy.exp(-1)) + start; # 1/e * v_max
        if not ((start < u_eff and u_eff < end) or (start > u_eff and u_eff > end)):
            import sys
            print >> sys.stderr, "IClamp seems to be used as sensor. use RecordingPoint instead -> get_tau_eff := 0"
            return 0
        
        #binary search
        start=0;
        end=len(voltage)-1;
        current = (end+start)/2;
        while(start+1 < end):#necessary, if value not contained
            if(voltage[current] < u_eff):
                if(voltage[start] < voltage[current]):
                    start = current;
                else:
                    end = current;
            elif(voltage[current] > u_eff):
                if(voltage[current] < voltage[end]):
                    end = current;
                else:
                    start = current;
            else:
                break;
            current = (end+start)/2;
        
        #slightly prefer shallow right side  
        if(numpy.abs(voltage[start] - u_eff) < numpy.abs(u_eff - voltage[end])):
            current = start;
        else:
            current = end;
        
        assert len(voltage) == len(time);
        tau = time[current] - self.iclamp.delay;
        if tau != self.get_tau_eff():
            import sys
            print >> sys.stderr, "τ_eff: %f is more accurate than %f" % (tau, self.get_tau_eff());
        return tau;
    
    def get_R_in(self):
        """
        original version
        """
        import numpy
        _, voltage = self.get_recording();
        volt_diff = max(voltage) - min(voltage);    # wenn es kurz ist?
        r_in = numpy.abs(float(volt_diff / self.iclamp.neuron_h_IClamp.amp));
        return r_in;
    
    def R_in(self):
        """
        R_in (input resistance): dV(max)/dI
        Is it extendable for more clamps?
        """
        import numpy
        u_diff = self.v[int(numpy.ceil((self.iclamp.delay+self.iclamp.duration)/self.recordingpoint.experiment.dt))] - self.v[int(self.iclamp.delay/self.recordingpoint.experiment.dt)];    #TODO: experiment too short
        i_diff   = self.iclamp.neuron_h_IClamp.amp - 0;
        r_in = float(u_diff) / i_diff; #physics -> sign is postivie
        assert r_in > 0;
        #TODO: assert r_in == self.get_R_in()	#print "r_in %f =?= %f get_R_in" %(r_in, self.get_R_in())#only holds for single experiments 
        return r_in;

    def tau_lin_fit(self, dt=1e-6):
        import scipy
        #import matplotlib.pyplot;
        #geht!! trotz fehlendem tau_eff a = (numpy.log(numpy.e - 1) - 1)/tau_eff;
        xdata, ydata  = self.get_recording(self.iclamp.delay, self.iclamp.delay+self.iclamp.duration);
        ydata = ydata - ydata[-1] - (ydata[-1]-ydata[-2]);
        ydata = scipy.log(ydata);
        #matplotlib.pyplot.plot(xdata, ydata);
        polycoeffs = scipy.polyfit(xdata, ydata, 1)
        yfit = scipy.polyval(polycoeffs, xdata)

        #matplotlib.pyplot.plot(xdata, yfit);
        #matplotlib.pyplot.show();

        #e^(at) = e^(-t/tau)
        return -1.0/(polycoeffs[0]);

    def tau_fit(self, fit_method='leastsq'):
        import neuron
        dhl = None#TODO:
        """
        Returns dictionary:
        - `p_lsq`: [tau, offset_0, asymptote]
        """
        return self.tau_lin_fit();
        #import datahandlinglib as dhl
        import scipy.optimize;
        time, voltage  = self.get_recording();
        t1_idx = (time >= self.clamp.delay);
        t2_idx = (time < (self.clamp.duration + self.clamp.delay));
        lidx_frame = t1_idx * t2_idx;#True * True = 1 , sonst 0
        time_frame = time[lidx_frame]
        voltage_frame = voltage[lidx_frame]
        tau0 = 0.5
        offset_00 = -10
        asymptote0 = -60
        ps0 = [tau0, offset_00, asymptote0]
        if fit_method == 'leastsq':
            plsq = scipy.optimize.leastsq(dhl.exp_error, ps0, args=(voltage_frame,time_frame-self.clamp.delay))[0]
        elif fit_method == 'fmin': 
            plsq = scipy.optimize.fmin(dhl.exp_error_2, ps0, args=(voltage_frame,time_frame-self.clamp.delay))
        else:
            import sys
            print >> sys.stderr, "ERROR: given fit-method not recognized"
            return {}
        fit_frame = dhl.exp_func_m(time_frame-self.clamp.delay, *plsq)
        return {
            'p_lsq': plsq,
            'time_frame': time_frame,
            'fit_frame': fit_frame}

    def plot(self, dt=0.001, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        """
        only one clamp allowed
        show and save impossible
        """
        import numpy
        import matplotlib.pyplot;
        t = numpy.array(self.t);
        u = numpy.array(self.v);
        matplotlib.pyplot.plot(t, u);#semilogy
        matplotlib.pyplot.xlabel("Time [ms]");
        matplotlib.pyplot.ylabel("Voltage [mV]");
        #matplotlib.pyplot.title('IClamp + fit');
        matplotlib.pyplot.axis('image');
        xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis();

        start = u[int(numpy.ceil((self.iclamp.delay+self.iclamp.duration)/self.experiment.dt))];
        end   = u[int(self.iclamp.delay/self.experiment.dt)];
        print str(end) +" - "+ str(start);
        t = numpy.arange(self.iclamp.delay, self.iclamp.delay+self.iclamp.duration+dt/2.0, dt);
        u = (end-start) * numpy.exp(- (t-self.iclamp.delay)/self.tau_eff ) + start;
        matplotlib.pyplot.plot(t, u);
        matplotlib.pyplot.axis([int(xmin), int(xmax + 0.5), int(ymin - 1.01), int(ymax + 1.01)]);
        xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis();

        print "τ_eff: "+ str(self.tau_eff+self.iclamp.delay) +", "+ str( (end-start)*(1-1/numpy.e)+start);
        pheight  = numpy.abs(start-end)/(ymax-ymin);
        pmin     = (min(start,end)-ymin)/(ymax-ymin);
        matplotlib.pyplot.axvline(x=self.iclamp.delay+self.iclamp.duration, ymin=pmin, ymax=pmin+pheight, label="$\Delta u$", color='r');
        matplotlib.pyplot.annotate("$\Delta u$", xy=(self.iclamp.delay+self.iclamp.duration-0.5, (start+end)/2.0) );
        matplotlib.pyplot.axvline(x=self.tau_eff+self.iclamp.delay, color='m');
        matplotlib.pyplot.axhline(y=(end-start)/numpy.e+start, color='m');
        matplotlib.pyplot.annotate("$(1 - \\frac{1}{e}) \cdot \Delta u + u_0$", xy=(self.tau_eff+self.iclamp.delay+0.2, (end-start)/numpy.e+start))
        matplotlib.pyplot.annotate("$\\tau_{eff} + t_{delay}$", xy=(self.tau_eff+self.iclamp.delay+0.2, ymin+0.05))
        matplotlib.pyplot.legend([
                "experimental $R_{in} = %.2f \\quad \\tau_{eff} = %.2f$" % (self.r_in, self.tau_eff), 
                "exponential fit $%.2f \\cdot e^{-\\frac{t - %.1f}{%.2f}} %.2f$" %  (self.r_in, self.iclamp.delay, self.tau_eff_fit, start)
            ], loc='best');

        #matplotlib.pyplot.axis('image');
        matplotlib.pyplot.grid(True, which='both', color='lightgrey')

        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format, transparent=True)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

    def __repr__(self):
        # Optional: for str(autoobj) or print autoobj    #has %s in %self.clamp, self.morphology
        return '<Result(r_in=%f, tau_eff=%f, tau_eff_fit=%f)>' % (
                self.r_in,
                self.tau_eff,
                self.tau_eff_fit,
                #str( self.tau_eff_fit  if self.tau_eff_fit != None else '') 
        )
            
    def __str__(self):
        return '%s\t%f\t%f\t%f' % (self.experiment.description, self.r_in, self.tau_eff, self.tau_eff_fit)

def test():
    import pydesignlib;
    cell = pydesignlib.MSO_Neuron();
    sim = pydesignlib.Simulation(cell);
    sim.set_IClamp();
    cell.set_passive_parameters(gp=0.004, E=-60, rho=200);
    cell.change_diameter(3);
    cell.change_length(150);
    sim.go();
    Rin = sim.get_Rin();
    tau_eff = sim.get_tau_eff();
    time, voltage = sim.get_recording();
    print (time, voltage);
    return Rin, tau_eff;

if __name__ == "__main__":
    import morphjongleur.model.neuron_passive
    import morphjongleur.model.morphology
    import morphjongleur.model.clamp

    #from neuron import gui; # For an overview (gui menu bar): Tools -> Model View -> 1 real cell -> root...
    #print test();

    ms = []
    ms.append( morphjongleur.model.morphology.Star(5,8) )
    ms.append( morphjongleur.model.morphology.Star(1,1,True, dendrite_diameter=3) )
    ms.append( morphjongleur.model.morphology.Star(1,1,False) )
    vts  = []
    tfs  = []
    for m in ms:
        r = RecordingPoint(m.getCompartment(0), position=0.5)
        c = morphjongleur.model.clamp.IClamp(m.getCompartment(0), position=0.5, amplitude=-1e-9,delay=1e-3,duration=3e-3)
                #simulation: , gp=0.004, E=-60, Ra=200, soma_Ra=1, soma_L=40, soma_diam=20, soma_nseg=1, dendrite_Ra=1, dendrite_length=150, dendrite_diameter=3, dendrite_nseg=1);
    
        e = Experiment(
            morphology=m, 
            recordingpoints=r, 
            clamps=c, 
            neuron_passive_parameter=morphjongleur.model.neuron_passive.Neuron_passive_parameter(Ra=1, g=0.004, e=-60), 
            duration=5e-3, dt=1e-4,
            description='passive channels'
        )
        assert len(e.recordingpoints) == 1
        assert len(e.clamps) == 1
        
        e.run_simulation()

        tau_fit = r.get_tau_fit(iclamp=c, recordingpoint=r)
        print(tau_fit)
        tau_fit.plot(e.dt, pic_dir='/tmp/', picture_formats=['png'])
        vts.append( r.get_voltage_trace() )
        tfs.append( tau_fit )

    vts[0].plot_ensemble( voltage_traces=vts, picture_file='/tmp/stars', picture_formats=['png'])
