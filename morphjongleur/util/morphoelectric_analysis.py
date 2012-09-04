#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import morphjongleur.model.morphology
import morphjongleur.model.neuron_passive
import morphjongleur.model.clamp
import morphjongleur.model.experiment

def experiment(morphology):#, amplitude
        recording_point  = morphjongleur.model.experiment.RecordingPoint(compartment=morphology.root_biggest_child)
       
        iclamp  = morphjongleur.model.clamp.IClamp(compartment=morphology.root_biggest_child,
                    amplitude=-1e-9, delay=0e-3, duration=3e-3
                )
        neuron_passive_parameter    = morphjongleur.model.neuron_passive.Neuron_passive_parameter(Ra=35.4,g=0.001)
        experiment  = morphjongleur.model.experiment.Experiment(
                        morphology=morphology, 
                        recording_points=[recording_point], clamps=[iclamp], 
                        neuron_passive_parameter=neuron_passive_parameter, 
                        duration=5e-3, dt=1e-4,
                        description = "TauExperiment @ %i\t%s " % (morphology.root_biggest_child.compartment_id, morphology.name),
                    )
        
        experiment.run_simulation()
        #voltage_trace   = experiment.get_voltage_trace(delay=2./frequency)# nach Einschwinphase
        tau_fit = recording_point.get_tau_fit(iclamp, recording_point)
        print "%f %f %f" % (tau_fit.get_R_in(), tau_fit.get_tau_eff(), tau_fit.tau_lin_fit()) 
        return tau_fit

if __name__ == '__main__':
    '''
    Parameter: files, not directories

    cd data/amira/
    pathto/metric_analysis.py *.swc
    '''
    import sys
    import morphjongleur.util.parser.swc

    print "name\tR_in\ttau_eff\ttau_eff_fit"
    for swc in sys.argv[1:]:#['../../data/test.swc']:#
        m   = morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False)
        #print m
        result = experiment(morphology=m)
        print result
        result.plot(picture_file='/tmp/taufit_%s' % (m.name), picture_formats=['png','svg'])

    #sys.exit(1)

    from morphjongleur.util.metric_analysis import MetricAnalysis
    bars= ['dorsal branch','ventral branch']
    xs  = ['$R_{in}$','$\\tau_{eff fit}$']
    results    = [experiment(morphology=morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False)) for swc in sys.argv[1:]]
    v   = {'dorsal branch': {'$R_{in}$':results[2].get_R_in()/results[1].get_R_in()-1,'$\\tau_{eff fit}$':results[2].tau_lin_fit()/results[1].tau_lin_fit()-1},
           'ventral branch':{'$R_{in}$':results[4].get_R_in()/results[3].get_R_in()-1,'$\\tau_{eff fit}$':results[4].tau_lin_fit()/results[3].tau_lin_fit()-1}
           }
    MetricAnalysis.bars_plot(v=v, bars=bars, xs=xs, colors=['#00ff00','#0000ff'], horizontal=True, tex=True, ratio=(16,9), picture_file='/tmp/change_tau', picture_formats=['png','svg'])#, y_label='change: forager / nurse - 1'
