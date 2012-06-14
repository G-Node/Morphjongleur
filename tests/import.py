# -*- coding: utf-8 -*-
'''
Created on 29.02.2012

@author: stransky
'''
import datajongleur
import time
import morphjongleur.util.parser.swc
from morphjongleur.orm.datajongleur.morphology import *
from morphjongleur.orm.datajongleur.neuron_single_cell_experiment import *
from morphjongleur.orm.datajongleur.clamp import *
from morphjongleur.orm.datajongleur.neuron import *

session = datajongleur.get_session()

def import_swcs(swcs):
    for swc in swcs:
        import_swc(swc)

def import_swc(swc):
    m_swc   = Morphology.swc_parse(swc)
    #m_swc   = Morphology.specify(m_swc)

    print("uploading data according to '%s'." %(swc)) ,
    start_time = time.time()
    m_swc.save()
    end_time = time.time()
    print("in %s seconds." %(end_time-start_time))
    
    print m_swc.uuid
    uuid = m_swc.uuid
    m_swc2 = Identity.load(uuid)
    print( type(m_swc2), m_swc2 )

    return m_swc


if __name__ == '__main__':
    import sqlalchemy
    print sqlalchemy.__version__
    assert sqlalchemy.__version__ == "0.7.4"    #sudo easy_install "sqlalchemy==0.7.4"
    import_swcs(['/home/stransky/git/mitsubachi/data/amira/Skeletontree_20120315_6-8.swc'])
