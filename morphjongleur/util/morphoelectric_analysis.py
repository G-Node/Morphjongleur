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
        recordingpoint  = morphjongleur.model.experiment.RecordingPoint(compartment=morphology.root)
       
        iclamp  = morphjongleur.model.clamp.IClamp(compartment=morphology.root,
                    amplitude=-1e-9, delay=0e-3, duration=3e-3
                )
        neuron_passive_parameter    = morphjongleur.model.neuron_passive.Neuron_passive_parameter(Ra=35.4,g=0.001)
        experiment  = morphjongleur.model.experiment.Experiment(
                        morphology=morphology, 
                        recordingpoints=[recordingpoint], clamps=[iclamp], 
                        neuron_passive_parameter=neuron_passive_parameter, 
                        duration=5e-3, dt=1e-4,
                        description = "TauExperiment @ %i\t%s " % (morphology.root.compartment_id, morphology.name),
                    )
        
        experiment.run_simulation()
        tau_fit = recordingpoint.get_tau_fit(iclamp)
        print "%f %f %f" % (tau_fit.get_R_in(), tau_fit.get_tau_eff(), tau_fit.tau_lin_fit()) 
        return tau_fit

def amplitudes(morphology):
        recordingpoints  = [morphjongleur.model.experiment.RecordingPoint(compartment=recordingpoint) for recordingpoint in morphology.compartments]
       
        iclamp  = morphjongleur.model.clamp.IClamp(compartment=morphology.root,
                    amplitude=-1e-9, delay=0e-3, duration=3e-3
                )
        neuron_passive_parameter    = morphjongleur.model.neuron_passive.Neuron_passive_parameter(Ra=35.4,g=0.001)
        experiment  = morphjongleur.model.experiment.Experiment(
                        morphology=morphology, 
                        recordingpoints=recordingpoints, clamps=[iclamp], 
                        neuron_passive_parameter=neuron_passive_parameter, 
                        duration=5e-3, dt=1e-4,
                        description = "TauExperiment @ %i\t%s " % (morphology.root.compartment_id, morphology.name),
                    )
        
        experiment.run_simulation()
        #voltage_trace   = experiment.get_voltage_trace(delay=2./frequency)# nach Einschwinphase
        return recordingpoints

def taus(morphology):
        recordingpoints  = [morphjongleur.model.experiment.RecordingPoint(compartment=recordingpoint) for recordingpoint in morphology.compartments]
       
        iclamp  = morphjongleur.model.clamp.IClamp(compartment=morphology.root,
                    amplitude=-1e-9, delay=0e-3, duration=3e-3
                )
        neuron_passive_parameter    = morphjongleur.model.neuron_passive.Neuron_passive_parameter(Ra=35.4,g=0.001)
        experiment  = morphjongleur.model.experiment.Experiment(
                        morphology=morphology, 
                        recordingpoints=recordingpoints, clamps=[iclamp], 
                        neuron_passive_parameter=neuron_passive_parameter, 
                        duration=5e-3, dt=1e-4,
                        description = "TauExperiment @ %i\t%s " % (morphology.root.compartment_id, morphology.name),
                    )
        
        experiment.run_simulation()
        taus = [recordingpoint.get_tau_fit(iclamp) for recordingpoint in recordingpoints]
        return taus

if __name__ == '__main__':
    '''
    Parameter: files, not directories

    cd data/amira/
    pathto/metric_analysis.py *.swc
    '''
    import sys
    import numpy
    import morphjongleur.util.parser.swc
    import morphjongleur.util.transformations
    import morphjongleur.model.morphology
    picture_formats = ['png','svg', 'pdf']

    print "direkt vergleich"

    ms = [morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False) for swc in sys.argv[1:]]
    from morphjongleur.util.metric_analysis import MetricAnalysis
    bars= ['dorsal branch','ventral branch']
    xs  = ['$R_{in}$','$\\tau_{\\mathrm{eff}}$']
    results    = [experiment(morphology=m) for m in ms]
    v   = {'dorsal branch': {'$\\tau_{\\mathrm{eff}}$':results[2].tau_lin_fit()/results[1].tau_lin_fit()-1,'$R_{in}$':results[2].get_R_in()/results[1].get_R_in()-1},
           'ventral branch':{'$\\tau_{\\mathrm{eff}}$':results[4].tau_lin_fit()/results[3].tau_lin_fit()-1, '$R_{in}$':results[4].get_R_in()/results[3].get_R_in()-1}
           }
    print v
    MetricAnalysis.bars_plot(v=v, bars=bars, xs=xs, colors=['#00ff00','#0000ff'], horizontal=True, tex=True, ratio=(8,6), xlimit=(-0.4,0.2), y_label='change: forager / nurse - 1', picture_file='/tmp/change_tau', picture_formats=picture_formats)#, y_label='change: forager / nurse - 1'


    # scaled nurse
    #need new morphologies as clamps are already set
    ms = [morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False) for swc in sys.argv[1:]]
    mas = [MetricAnalysis(m) for m in ms]

    scale_nurse_dorsal   = mas[2].total_length / mas[1].total_length
    nurse_dorsal_scaled  = ms[1].scale(scale_nurse_dorsal)
    print "DB total length nurse %f vs. forager %f -scale: %f-> nurse %f" % (mas[1].total_length , mas[2].total_length, scale_nurse_dorsal , MetricAnalysis(nurse_dorsal_scaled).total_length)

    scale_nurse_ventral   = mas[4].total_length / mas[3].total_length
    nurse_ventral_scaled  = ms[3].scale(scale_nurse_ventral)
    print "VB total length nurse %f vs. forager %f -scale: %f-> nurse %f" % (mas[3].total_length , mas[4].total_length, scale_nurse_ventral, MetricAnalysis(nurse_ventral_scaled).total_length)


    results    = [experiment(morphology=m) for m in [ms[0], nurse_dorsal_scaled, ms[2], nurse_ventral_scaled, ms[4]]]
    v   = {'dorsal branch': {'$\\tau_{\\mathrm{eff}}$':results[2].tau_lin_fit()/results[1].tau_lin_fit()-1,'$R_{in}$':results[2].get_R_in()/results[1].get_R_in()-1},
           'ventral branch':{'$\\tau_{\\mathrm{eff}}$':results[4].tau_lin_fit()/results[3].tau_lin_fit()-1, '$R_{in}$':results[4].get_R_in()/results[3].get_R_in()-1}
           }
    print v
    MetricAnalysis.bars_plot(v=v, bars=bars, xs=xs, colors=['#00ff00','#0000ff'], horizontal=True, tex=True, ratio=(8,7), xlimit=(-0.4,0.2), y_label='change: forager / nurse - 1', picture_file='/tmp/change_tau_scaled_nurse', picture_formats=picture_formats)#, y_label='change: forager / nurse - 1'


    # scaled forager
    #need new morphologies as clamps are already set
    ms = [morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False) for swc in sys.argv[1:]]
    mas = [MetricAnalysis(m) for m in ms]

    scale_forager_dorsal   = mas[1].total_length / mas[2].total_length
    forager_dorsal_scaled  = ms[2].scale(scale_forager_dorsal)
    print "DB total length nurse %f vs. forager %f -scale: %f-> forager %f" % (mas[1].total_length , mas[2].total_length, scale_forager_dorsal, MetricAnalysis(forager_dorsal_scaled).total_length)

    scale_forager_ventral   = mas[3].total_length / mas[4].total_length
    forager_ventral_scaled  = ms[4].scale(scale_forager_ventral)
    print "VB total length nurse %f vs. forager %f -scale: %f-> forager %f" % (mas[3].total_length , mas[4].total_length, scale_forager_ventral, MetricAnalysis(forager_ventral_scaled).total_length)


    results    = [experiment(morphology=m) for m in [ms[0], ms[1], forager_dorsal_scaled, ms[3], forager_ventral_scaled]]
    v   = {'dorsal branch': {'$\\tau_{\\mathrm{eff}}$':results[2].tau_lin_fit()/results[1].tau_lin_fit()-1,'$R_{in}$':results[2].get_R_in()/results[1].get_R_in()-1},
           'ventral branch':{'$\\tau_{\\mathrm{eff}}$':results[4].tau_lin_fit()/results[3].tau_lin_fit()-1, '$R_{in}$':results[4].get_R_in()/results[3].get_R_in()-1}
           }
    print v
    MetricAnalysis.bars_plot(v=v, bars=bars, xs=xs, colors=['#00ff00','#0000ff'], horizontal=True, tex=True, ratio=(8,6), xlimit=(-0.4,0.2), y_label='change: forager / nurse - 1', picture_file='/tmp/change_tau_scaled_forager', picture_formats=picture_formats)#, y_label='change: forager / nurse - 1'

    # scaled both
    #need new morphologies as clamps are already set
    ms = [morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False) for swc in sys.argv[1:]]
    mas = [MetricAnalysis(m) for m in ms]

    nurse_dorsal_scaled  = ms[1].scale(scale_nurse_dorsal)
    nurse_ventral_scaled  = ms[3].scale(scale_nurse_ventral)
    print "DB total length scaled nurse %f vs. forager %f" % (MetricAnalysis(nurse_dorsal_scaled).total_length, MetricAnalysis(forager_dorsal_scaled).total_length)

    forager_dorsal_scaled  = ms[2].scale(scale_forager_dorsal)
    forager_ventral_scaled  = ms[4].scale(scale_forager_ventral)
    print "VB total length scaled nurse %f vs. forager %f" % (MetricAnalysis(nurse_ventral_scaled).total_length, MetricAnalysis(forager_ventral_scaled).total_length)


    results    = [experiment(morphology=m) for m in [ms[0], nurse_dorsal_scaled, forager_dorsal_scaled, nurse_ventral_scaled, forager_ventral_scaled]]
    v   = {'dorsal branch': {'$\\tau_{\\mathrm{eff}}$':results[2].tau_lin_fit()/results[1].tau_lin_fit()-1,'$R_{in}$':results[2].get_R_in()/results[1].get_R_in()-1},
           'ventral branch':{'$\\tau_{\\mathrm{eff}}$':results[4].tau_lin_fit()/results[3].tau_lin_fit()-1, '$R_{in}$':results[4].get_R_in()/results[3].get_R_in()-1}
           }
    print v
    MetricAnalysis.bars_plot(v=v, bars=bars, xs=xs, colors=['#00ff00','#0000ff'], horizontal=True, tex=True, ratio=(8,6), xlimit=(-0.4,0.2), y_label='change: forager / nurse - 1', picture_file='/tmp/change_tau_scaled_both', picture_formats=picture_formats)


    sys.exit(0)

    print "name\tR_in\ttau_eff\ttau_eff_fit"
    for swc in sys.argv[1:]:#['../../data/test.swc']:#
        m   = morphjongleur.model.morphology.Morphology.swc_parse(swc, verbose=False)
        print m.name, 
        try:
            result = experiment(morphology=m)
            print result
            result.plot(picture_file='/tmp/taufit_%s' % (m.name), picture_formats=picture_formats)
            
            rps = amplitudes(morphology=m)
            vts = [rp.get_voltage_trace() for rp in rps]
            ts  = taus(morphology=m)

            #for c in m.compartments:
            #    c.color = '0'
            amplitude_max   = float(numpy.max([vt.amplitude for vt in vts]))
            amplitude_min   = 0.
            for vt in vts:
                vt.recordingpoint.compartment.color = "%f" % (vt.amplitude / amplitude_max)
            morphjongleur.model.morphology.Compartment.plot_color([m.compartments], picture_file='/tmp/amplitude_grayscale_%s' % (m.name), picture_formats=picture_formats)

            r_in_max = float(numpy.max([t.r_in for t in ts]))
            r_in_min = 0.
            for t in ts:
                t.recordingpoint.compartment.color = "%f" % (t.r_in / r_in_max)
            morphjongleur.model.morphology.Compartment.plot_color([m.compartments], picture_file='/tmp/r_in_grayscale_%s' % (m.name), picture_formats=picture_formats)
            tau_eff_max     = float(numpy.max([t.tau_eff for t in ts]))
            tau_eff_min     = 0.
            for t in ts:
                t.recordingpoint.compartment.color = "%f" % (t.tau_eff / tau_eff_max)
            morphjongleur.model.morphology.Compartment.plot_color([m.compartments], picture_file='/tmp/tau_eff_grayscale_%s' % (m.name), picture_formats=picture_formats)
            tau_eff_fit_max = float(numpy.max([t.tau_eff_fit for t in ts]))
            tau_eff_fit_min = 0.
            for t in ts:
                t.recordingpoint.compartment.color = "%f" % (t.tau_eff_fit / tau_eff_fit_max)
            morphjongleur.model.morphology.Compartment.plot_color([m.compartments], picture_file='/tmp/tau_eff_fit_grayscale_%s' % (m.name), picture_formats=picture_formats)

        except Exception, e:
            print swc
            import traceback 
            print traceback.format_exc()