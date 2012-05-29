# -*- coding: utf-8 -*-
from mrj.io.database import Database
from mrj.model.morphology import *
import os

db = Database(
    db_name='postgresql://hal08.g-node.pri/lehmann',
    exec_role='morphjokey_admin',
    exec_path='morphjokey')
m = db.load_morphology( 291 )
print( m )
print( m.compartments[2] )

mg  = Morphology_groups()
mg.name = 'groupname'
mg.age  = '23'
mg.description = 'cool'

m._groups.append( mg )
m.description  += ' owner:'+str( os.environ['USER'] ) # use your name!

cg  = Compartment_groups()
cg.name = 'primes'
cg.description = 'cool'
cg.type = 'axon'

print( m )    # key not available
print( m.compartments[2] )

#append OR store sufficient to get key
m.compartments[2]._groups.append(cg)
m.compartments[3]._groups.append(cg)
db.store( m )

print( m )    # key available
print( m.compartments[2] )

#identic
#m   = db.load_morphology( m.morphology_key )
#print( m )
#print( m.compartments[2] )
