# -*- coding: utf-8 -*-
from mrj.io.swc import SwcParser
from mrj.io.database import Database

db = Database(
    db_name='postgresql://hal08.g-node.pri/lehmann',
    exec_role='morphjokey_admin',
    exec_path='morphjokey')
swc_parser  = SwcParser(db)

m   = swc_parser.parse('../../../Data/example_morph.swc')
#print( m )
#print( m._compartments[0])

db.store( m )
#db.create_subtree( m )
# morphology_key mapped via sqlalchemy:
print m.morphology_key
m = db.load_morphology( m.morphology_key )

print( m )  # contains ALL auto calculated data from database
#print( m._compartments[0])

#print "\n\n\ncompartments as dicts:"
#for c in m._compartments:
#    print c.info.__dict__

#a = db.session.execute("SELECT compartment_id,compartment_parent_id,children FROM morphjokey.mrj_v_compartments_children WHERE morphology_key = 150;")
#b = a.fetchall()
