# -*- coding: utf-8 -*-
from mrj.io.swc import SwcParser
from mrj.io.database import Database

'''
Created on 29.02.2012

@author: stransky
'''

db = Database(
    db_name='postgresql://hal08.g-node.pri/morphjongleur',
    exec_role='morphjokey_admin',
    exec_path='mitsubachi')
swc_parser  = SwcParser(db)

def import_swcs(swcs):
    for swc in swcs:
        import_swc(swc)

def import_swc(swc):
    m_swc   = swc_parser.parse(swc)
    #m_swc.create_tree()

    print("uploading data according to '%s'." %(swc)) ,
    import time
    start_time = time.time()
    db.store( m_swc )
    db.session.commit()
    end_time = time.time()
    print("in %s seconds." %(end_time-start_time))
    
    print m_swc.morphology_key
    m_swc = db.load_morphology( m_swc.morphology_key )
    print( m_swc )

    return m_swc 


if __name__ == '__main__':
    import_swcs([
         '../../../data/H060602DB_10_2_zentai_.swc', 
         '../../../data/H060602VB_10_2_zentai_.swc', 

         '../../../data/H060607DB_10_2(zentai).swc', 
         '../../../data/H060607VB_10_2(zentai).swc'
    ])
