# -*- coding: utf-8 -*-
from mrj.io.database import Database
from mrj.model.neuromorph import Morph

db = Database(
    db_name='postgresql://hal08.g-node.pri/lehmann',
    exec_role='morphjokey_admin',
    exec_path='morphjokey')
# must be mapped before Object is created
import mrj.orm.morphology
mapping = mrj.orm.morphology.Mapper( db.engine )
mapping.orm_map()
m_db = db.load_morphology( 460 )

m = Morph(m_db)
print(m.get_surface_volume_ratio())
print(m.db_model.analysis['surface_volume_ratio'])
print(m.db_model.analysis)
print m

m.show_tree()
