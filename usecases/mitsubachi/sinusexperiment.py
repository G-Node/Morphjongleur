'''
Created on 08.05.2012

@author: stransky
'''
import mrj.model.neuron_passive
import mrj.model.clamp
import mrj.util.pattern_generator
import mrj.model.experiment
import mrj.model.sinusresult
#import mrj.dj.sinusresult

def experiment(morphology, compartment, frequency, amplitude):
        vclamp  = mrj.model.clamp.VClamp(compartment=morphology.biggest)
       
        iclamp  = mrj.util.pattern_generator.SinusClamp(compartment=compartment,
                    amplitude=amplitude,
                    frequency=frequency, duration=3./frequency
                )
        neuron_passive_parameter    = mrj.model.neuron_passive.Neuron_passive_parameter(Ra=35.4,g=0.001)
        experiment  = mrj.model.experiment.Experiment(
                        morphology, 
                        vclamps=[vclamp], iclamps=[iclamp], 
                        neuron_passive_parameter=neuron_passive_parameter, 
                        duration=3./frequency,
                        dt=1e-4,
                        description = "SinusExperiment %s %i %f" % (morphology.name, compartment.compartment_id, frequency),
                    )
        experiment.run_simulation()
        voltage_trace   = experiment.get_voltage_trace(delay=2./frequency)# nach Einschwinphase
        return mrj.model.sinusresult.SinusResult(voltage_trace, compartment, frequency, amplitude)

if __name__ == '__main__':
    import mrj.io.swc
    from mrj.io.database import Database
    db = Database(
                db_name='postgresql://hal08.g-node.pri/morphjongleur',
                exec_role='morphjokey_admin',
                exec_path='mitsubachi'
            )
    # must be mapped before Object is created
    import mrj.orm.sinusresult
    import mrj.model.morphology
    mapping = mrj.orm.sinusresult.Mapper( db )
    mapping.orm_map()

    for swc in ['../../../data/H060602DB_10_2_zentai_.swc','../../../data/H060602VB_10_2_zentai_.swc','../../../data/H060607DB_10_2(zentai).swc','../../../data/H060607VB_10_2(zentai).swc']:#['../../../data/test.swc']:
        print swc
        morphology   = mrj.model.morphology.swc_parse(swc)
        cs  = []#morphology.biggest
        cs.extend( morphology.terminal_tips )
        for f in range(10,1001,15):
            print "%i Hz" % (f)
            i = 0
            for c in cs:
                i = i + 1 
                print "%i/%i\t\r" % (cs.index(c)+1, len(cs)), 
                sinusresult = experiment(morphology=morphology, compartment=c, frequency=f, amplitude=-1e-9)
                #sinusresult.save()
                db.store(sinusresult)
