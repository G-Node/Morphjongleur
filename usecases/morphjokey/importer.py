#!env python
# -*- coding: utf-8 -*-
'''
Created on 22.02.2011

@author: stransky
'''

import time
from mrj.io.swc import SwcParser
from mrj.orm.morphology import *
from mrj.io.database import Database

def import_swc(swcs):
    db = Database(
        db_name   = 'postgresql://hal08.g-node.pri/lehmann',
        exec_role = 'morphjokey_admin',
        exec_path = 'morphjokey'
    )
    swc_parser  = SwcParser()

    for swc in swcs:
        print (swc)
        m_swc   = swc_parser.parse(swc)
        db.store( m_swc )
        print("start to commit data according to '%s'." %(swc))
        tic = time.time()
        db.session.commit()
        toc = time.time()
        print("Done. It took %s seconds." %(toc-tic))
    

if __name__ == "__main__":
    '''
    Parameter: files, not directories
    commited at end -> resetable, but keep ram in mind
    

    cd Data/Morphologies/
    ../../py_morphjokey/usecases/importer.py P??/*
    '''
    import sys
    import_swc(sys.argv[1:])
