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
    root_dir = '/home/philipp/Repositories/G-Node/mitsubachi/'
    #root_dir = '/home/stransky/git/mitsubachi/'
    #root_dir = '../../../'
    import_swcs([
         #root_dir + 'data/H060602DB_10_2_zentai_.swc', 
         #root_dir + 'data/H060602VB_10_2_zentai_.swc', 
         #root_dir + 'data/H060607DB_10_2(zentai).swc', 
         #root_dir + 'data/H060607VB_10_2(zentai).swc'
    ])
