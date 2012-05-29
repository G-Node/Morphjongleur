#!env python
# -*- coding: utf-8 -*-
'''
Created on 09.06.2011

@author: stransky
'''

from mrj.io.swc import SwcParser
from mrj.model.morphology import *
from mrj.io.database import Database
from sqlalchemy import Table
from sqlalchemy.sql import select
#from sqlalchemy.dialects.postgresql import ARRAY, INTEGER

db = Database(
    db_name   = 'postgresql://hal08.g-node.pri/lehmann',
    exec_role = 'morphjokey_admin',
    exec_path = 'morphjokey'
)
db_dendrit_info = Database(
    db_name   = 'postgresql://hal09.g-node.pri/bayer',
    exec_role = 'stransky',
    exec_path = 'philipp_data_mso'  # does sometimes not work, so define schema.
)
dendrit_info_table = Table('mso52lat_infos', db_dendrit_info.metadata, autoload=True, 
    #autoload_with=db_dendrit_info.engine,    #recommended for console call, but does not work 
    schema='philipp_data_mso')

swc_parser  = SwcParser(db)

def import_swcs(swcs):
    for swc in swcs:
        import_swc(swc)

def import_swc(swc):
    m_swc   = swc_parser.parse(swc)
    #m_swc.create_tree()

    s = select([dendrit_info_table], dendrit_info_table.c.label == m_swc.name)
    result = s.execute()
    #if len(result) != 1:           print "no information available"

    #generate subtree & add groups
    for row in result:
        for lateral in row[7]:    #'first_dend_lateral_segment_id'
            cg  = Compartment_groups()
            cg.name = 'lateral'
            cg.description = 'below #'+str(lateral)
            cg.type = 'dendrit'
            for c in m_swc.getSubtree( int(lateral) ):
                c._groups.append(cg)
                #print( c )
        for medial in row[8]:    #'first_dend_medial_segment_id integer' 
            cg  = Compartment_groups()
            cg.name = 'medial'
            cg.description = 'below #'+str(medial)
            cg.type = 'dendrit'
            for c in m_swc.getSubtree( int(medial) ):
                c._groups.append(cg)
                #print( c )

        print("uploading data according to '%s'." %(swc)) ,
        import time
        start_time = time.time()
        db.store( m_swc )
        db.session.commit()
        end_time = time.time()
        print("in %s seconds." %(end_time-start_time))
        
        return m_swc 

def import_msos(swcs, description):
    m_swcs  = []

    mgWAxon  = Morphology_groups()
    mgWAxon.name = 'Axon'
    mgWAxon.description = 'with'

    mgWOAxon  = Morphology_groups()
    mgWOAxon.name = 'Axon'
    mgWOAxon.description = 'without'
    
    mas = {}
    
    for swc in swcs:
        age = swc[swc.find('P')+1:swc.find('P')+3]
        if mas.has_key(age):
            mgAge = mas[age]
        else:
            mgAge  = Morphology_groups()
            mgAge.name = 'age'
            mgAge.description = str(age)
            mgAge.age = int(age)
            mas[age]    = mgAge

        m_swc   = import_swc(swc)
        m_swc.description   = description
        m_swc._groups.append( mgAge )

        if swc.find('ohne') == -1:
            m_swc._groups.append( mgWAxon )
        else:
            m_swc._groups.append( mgWOAxon )
        db.session.add( m_swc )

        m_swcs.append( m_swc )
        #print m_swc
    db.session.commit()#TODO: to ensure in same groups!
    return m_swcs
    

if __name__ == "__main__":
    '''
    Parameter: files, not directories
    commited at end -> resetable, but keep ram in mind
    
    #cd /
    python description src/mp/usecases/import_mso.py Data/MSO_Neurone/P??/*
    '''
    import sys
    import_msos(sys.argv[2:], sys.argv[1])

#    ms   = import_msos(['../../../Data/MSO_Neurone/P10/SkeletonTree_0001_1-1G.swc', '../../../Data/MSO_Neurone/P10/SkeletonTree_0001_1-1G.swc', '../../../Data/MSO_Neurone/P10/SkeletonTree_0002_2-3R.swc', '../../../Data/MSO_Neurone/P10/ohne_axon_SkeletonTree_0002_2-1G.swc'], "test")
#    print ms[0]
