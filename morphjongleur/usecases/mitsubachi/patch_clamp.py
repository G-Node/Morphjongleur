# -*- coding: utf-8 -*-
'''
@author: stransky
'''
#from mrj.model.pydesignlib import MSO_Neuron, Simulation
from mrj.io.database import Database
from mrj.model.neuron_passive import *
import mrj.model.clamp
import mrj.util.pattern_generator
import mrj.model.experiment
import mrj.orm.experiment 
#from neuron import gui; # For an overview (gui menu bar): Tools -> Model View -> 1 real cell -> root...

def generate(morphologies=[16,17,18,19]):    
    for m in morphologies:
    #  try:
        db = Database(
            db_name='postgresql://hal08.g-node.pri/morphjongleur',
            exec_role='morphjokey_admin',
            exec_path='mitsubachi')
        
        # must be mapped before Object is created
        #mapping  = mrj.orm.experiment.Mapper( db.engine )
        mapping     = mrj.orm.morphology.Mapper( db.engine )
        mapping.orm_map()
        
        morphology = db.load_morphology( m )
        print morphology 

        iclamp  = mrj.model.clamp.IClamp(morphology.biggest)
        r  = mrj.model.experiment.RecordingPoint(morphology.compartments[len(morphology.compartments)-1])

        experiment  = mrj.model.experiment.Experiment(
                        morphology, 
                        clamps=[iclamp], recording_points=r, 
                        neuron_passive_parameter=Neuron_passive_parameter(Ra=35.4,g=0.001)
                    )
        db.store( experiment )
        print experiment

        experiment.run_simulation()
        voltage_trace   = experiment.get_voltage_trace()
        db.store( voltage_trace )
        print voltage_trace
        
    #  except Exception, e:
    #    import traceback
    #    print traceback.format_exc()

def plot_patch_clamp(experiment_key, picture_file='/tmp/voltage_trace_f-265_iclamp', picture_formats=['png']):
    db = Database(
    db_name='postgresql://hal08.g-node.pri/morphjongleur',
    exec_role='morphjokey_admin',
    exec_path='mitsubachi')
    
    # must be mapped before Object is created
    mapping  = mrj.orm.experiment.Mapper( db.engine )
    mapping.orm_map()
    e = db.load_experiment( experiment_key )
    print e

    print e.voltage_trace
    e.voltage_trace.plot(picture_file=picture_file+str(e.morphology), picture_formats=picture_formats)
    #experiment.plot_fit(r.r_in, r.tau_eff, r.tau_eff_fit)

def sinusclamp(morphology, rounds=3, frequency=265, picture_file='/tmp/voltage_trace_f-265_sinusclamp_', picture_formats=['png']):
    iclamp   = mrj.util.pattern_generator.SinusClamp(compartment=morphology.biggest,
                    frequency=frequency, duration=rounds/float(frequency)
                )
    r  = mrj.model.experiment.RecordingPoint(morphology.compartments[len(morphology.compartments)-1])

    experiment  = mrj.model.experiment.Experiment(
                    morphology=morphology, 
                    clamps=[iclamp], recording_points=r, 
                    neuron_passive_parameter    = Neuron_passive_parameter(Ra=35.4,g=0.001), 
                    duration=rounds/float(frequency)
                )
    print experiment
    experiment.run_simulation()
    voltage_trace   = r.get_voltage_trace()
    print voltage_trace
    voltage_trace.plot(picture_file=picture_file+str(morphology.name), picture_formats=picture_formats)

def sinus_maxima(morphologies, rounds=3, frequency=265, xlim=None, ylim=None, picture_file='/tmp/', picture_formats=['png']):
    import numpy
    maxima  = {}
    minima  = {}
    titles  = {} 
    for morphology in morphologies:
       
        voltage_traces = {}
        cs  = [morphology.biggest]
        cs.extend( morphology.leafs )
        print len(cs)
        
        print cs
        leafs    = {}
        for c in cs:
            leafs[c.radius]    = c
        cs  = leafs.values()#TODO: remove pseudo representative choice
        print len(cs)
        print cs

        maxima[morphology]   = []
        minima[morphology]   = []
        titles[morphology]   = morphology.name
        for c in cs:
            recording_point  = mrj.model.experiment.RecordingPoint(compartment=morphology.biggest)
            print str(cs.index(c)+1)+"/"+str(len(cs))
            iclamp  = mrj.util.pattern_generator.SinusClamp(compartment=c, 
                        frequency=frequency, duration=rounds/float(frequency)
                    )
            neuron_passive_parameter    = Neuron_passive_parameter(Ra=35.4,g=0.001)
            experiment  = mrj.model.experiment.Experiment(
                            morphology, 
                            clamps=[iclamp], recording_points=[recording_point], 
                            neuron_passive_parameter=neuron_passive_parameter, 
                            duration=rounds/float(frequency),
                            dt=numpy.power(10, - (numpy.ceil(numpy.log10(265))+1) ),#one magnitue greater: e.g. 265 -> 1e-4, 
                            description = u"%f µm @ %i" % (c.radius, c.compartment_id)
                        )
            #print experiment
            experiment.run_simulation()
            voltage_trace   = recording_point.get_voltage_trace()
            maxima[morphology].append(voltage_trace.v_max)
            minima[morphology].append(voltage_trace.v_min)
            #print voltage_trace
            voltage_traces[c]    = voltage_trace
        voltage_traces[morphology.biggest].plot_ensemble(voltage_traces.values(), title="Voltagetrace in %s for sinusoidal inputsignal at %i Hz" % (morphology.name, frequency), ylim=ylim, picture_file=picture_file+'voltage_trace_f-'+str(frequency)+'_'+str(morphology.name), picture_formats=picture_formats)
    plot_daempfung(maxima.values(),titles=titles.values(), picture_file=picture_file+'/daempfung_maxima_f-265_'+str(morphology.name), picture_formats=picture_formats)
    plot_daempfung(minima.values(),titles=titles.values(), picture_file=picture_file+'/daempfung_minima_f-265_'+str(morphology.name), picture_formats=picture_formats)
    
    
def plot_daempfung(data = [[1,2,3,4,5],[1,2,3,4,5]], titles=["test", 'affe'], picture_file=None, picture_formats=['png']):
    import matplotlib.pyplot

    matplotlib.pyplot.title('Signaldaempfung')
    
    matplotlib.pyplot.boxplot(data, positions=range(len(data)))#, notch=0, sym='+', vert=1, whis=1.5
    #, color='black'
    matplotlib.pyplot.xticks(range(len(data)), titles)# ,rotation=18

    if(picture_file != None):
        for picture_format in picture_formats:
            matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format)
    else:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()


if __name__ == '__main__':
    #generate()
    #e = plot_patch_clamp(27)
    import mrj.io.swc
    swc   = '../../../data/test.swc'
    #swc   = '../../../data/H060607VB_10_2(zentai).swc'
    #sinusclamp(mrj.model.morphology.Morphology.swc_parse(swc))
    #swcs=['../../../data/test.swc', '../../../data/test.swc']
    swcs=['../../../data/H060602DB_10_2_zentai_.swc','../../../data/H060602VB_10_2_zentai_.swc','../../../data/H060607DB_10_2(zentai).swc','../../../data/H060607VB_10_2(zentai).swc']
    sinus_maxima(mrj.model.morphology.Morphology.swcs_parse(swcs), ylim=(-20, -100), picture_formats=['svg'])
    #plot_daempfung(picture_file='/tmp/d')
