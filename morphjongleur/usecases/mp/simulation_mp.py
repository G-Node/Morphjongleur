# -*- coding: utf-8 -*-
'''
@author: stransky
'''

from mrj.io.database import Database
from mp.model.experiment_mso import Experiment_MSO
from mp.model.neuron_mso import Neuron_MSO
import time
#from sqlalchemy.dialects.postgresql import ARRAY

db = Database(
    db_name='postgresql://hal08.g-node.pri/lehmann',
    exec_role='morphjokey_admin',
    exec_path='morphjokey'
)

# must be mapped before Object is created
import mp.orm.experiment_mso
mapping = mp.orm.experiment_mso.Mapper( db.engine )
mapping.orm_map()

ms          = [500, 501, 505, 506, 507, 512, 513, 514, 515, 516, 517, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 534, 535, 536, 537, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551] # range(500, 551+1) # range(392, 443+1) # [244] #
synapses    = [10]          # range(10, 51, 5)
rs          = [1,0.93,0.64,0.32,0] # [1,0]
fs          = [500]         # 500 Hz [1/s]
duration    = 10e-3         # 50-100 ms
iterations  = 10            # 10
gns         = ['both', 'medial', 'lateral']
description = 'exitatory'
resumepoint = None  # (0,424,1,10,500,'both') #
start_time  = time.time()

#from neuron import gui
#print "morphology type r -> sigma synapsen minU[mV] at[ms] calc_time[s]"
print "resumepoint(iteration, morphology_key, r, #synapses, frequency, group_name)\tinit time + exp time + store time = run time [s]\tResult"
for m in ms:
    for i in range(iterations):
        if resumepoint == None:
            morphology   = db.load_morphology( m )
        for r in rs:
            for n in synapses:
                for f in fs:
                    for gn in gns:
                        run_start   = start_time
                        print "(%i,%i,%f,%i,%f,%s)\t" % (i,m,r,n,f,gn) ,
                        if resumepoint != None:
                            if resumepoint == (i,m,r,n,f,gn):
                                morphology  = db.load_morphology( m )
                                resumepoint = None
                            else: 
                                print "skipped"
                                continue

                        mso_neuron   = Neuron_MSO( morphology )
                        experiment = Experiment_MSO(mso_neuron, f=f, r=r, n=n, group_name=gn, duration=duration, description=description)
                        #print experiment
                        db.store( experiment )
                        experiment.neuron_create()

                        end_time = time.time()
                        print "%f" %(end_time-start_time) ,
                        start_time = end_time
                        
                        experiment.run_simulation()

                        end_time = time.time()
                        experiment.calc_time = end_time-start_time
                        print "+ %f" %(experiment.calc_time) ,
                        start_time = end_time

                        result  = experiment.get_result()
                        #print result
                        db.store( result )
                        ts      = str( result.t )
                        vs      = str( result.v )
                        syn_comp_ids = str( experiment.synapses_on_compartment_id )
                        syn_rs       = str( experiment.rs )
                        
#from mp.model.synapse import getSigma
#import math
#                        sigma   = getSigma(r)
#                        if math.isinf( sigma ):
#                            db.session.execute("INSERT INTO morphjokey.mp_results(morphology_key,dendrite, r,sigma,frequency,n,synapses,rs, v_min,t_min,duration,t,v, description,i,calc_time) VALUES(%i,'%s', %f,'inf',%f,%i,'{%s}','{%s}', %f,%f,%f,'{%s}','{%s}', '%s',%i,%f);" %(m,str(gn),r,       f,n,syn_comp_ids[1:-1],syn_rs[1:-1], result.v_min,result.t_min, duration,ts[1:-1],vs[1:-1], experiment.description,i,end_time-run_start))
#                        else:
#                            db.session.execute("INSERT INTO morphjokey.mp_results(morphology_key,dendrite, r,sigma,frequency,n,synapses,rs, v_min,t_min,duration,t,v, description,i,calc_time) VALUES(%i,'%s', %f,%f   ,%f,%i,'{%s}','{%s}', %f,%f,%f,'{%s}','{%s}', '%s',%i,%f);" %(m,str(gn),r,sigma, f,n,syn_comp_ids[1:-1],syn_rs[1:-1], result.v_min,result.t_min, duration,ts[1:-1],vs[1:-1], experiment.description,i,end_time-run_start))
#                        db.session.commit()

                        mso_neuron.remove_synapses()

                        end_time = time.time()
                        print "+ %f" %(end_time-start_time) ,
                        start_time = end_time

                        print "= %f" %(end_time-run_start) ,
                        print result
