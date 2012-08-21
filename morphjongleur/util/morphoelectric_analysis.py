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
        recording_point  = morphjongleur.model.experiment.RecordingPoint(compartment=morphology.biggest)
       
        iclamp  = morphjongleur.model.clamp.IClamp(compartment=morphology.biggest,
                    amplitude=-1e-9, delay=0e-3, duration=3e-3
                )
        neuron_passive_parameter    = morphjongleur.model.neuron_passive.Neuron_passive_parameter(Ra=35.4,g=0.001)
        experiment  = morphjongleur.model.experiment.Experiment(
                        morphology=morphology, 
                        recording_points=[recording_point], clamps=[iclamp], 
                        neuron_passive_parameter=neuron_passive_parameter, 
                        duration=5e-3, dt=1e-4,
                        description = "TauExperiment @ %i\t%s " % (morphology.biggest.compartment_id, morphology.name),
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
        #result.plot(pic_dir='/tmp/', picture_formats=['png'])
