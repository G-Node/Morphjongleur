# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import neuron
#import scipy
#from mrj.model.pydesignlib import MSO_Neuron, Simulation
from mrj.io.database import Database
from mrj.model.neuron_passive import *
from mrj.model.morphology import *
from mrj.model.clamp import *
from mrj.model.experiment import *
from mrj.model.experiment import plot
import mrj.orm.experiment 
from neuron import gui; # For an overview (gui menu bar): Tools -> Model View -> 1 real cell -> root...

db = Database(
    db_name='postgresql://hal08.g-node.pri/lehmann',
    exec_role='morphjokey_admin',
    exec_path='morphjokey')

# must be mapped before Object is created
mapping  = mrj.orm.experiment.Mapper( db.engine )
mapping.orm_map()

morphology = db.load_morphology( 256 )
print morphology 

#morphology.create_tree()
clamp   = IClamp(morphology.root.children[0])

neuron_passive_parameter    = Neuron_passive_parameter(Ra=80)
experiment  = Experiment(morphology, clamp, neuron_passive_parameter)
print experiment
db.store( experiment )

#TODO: change compartments needs group group mapping
#c   = m.getCompartment( 5 )
#c   = m.getCompartment( 7 )

experiment.neuron_create()
experiment.run_simulation()
result   = experiment.get_result()
print result
db.store( result )

e = db.load_experiment( experiment.experiment_key )
print e

e.run_simulation()
r   = e.get_result()
print r

experiment.plot_fit(r.r_in, r.tau_eff, r.tau_eff_fit)
plot( [experiment] )    # [e] works, (e) is flattend
