# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import datajongleur
from morphjongleur.orm.datajongleur.morphology import *
from morphjongleur.orm.datajongleur.neuron_single_cell_experiment import *
from morphjongleur.orm.datajongleur.clamp import *
from morphjongleur.orm.datajongleur.neuron import *

session = datajongleur.get_session()

def insert(dir = '/home/stransky/git/mitsubachi/data/'):

    import morphjongleur.util.parser.swc
    morphologies    = []
    keys            = []
    for swc in [
'Skeletontree_20120315_6-8.swc',
'SkeletonTree_20120330_RIM_1-2-3.swc',
'SkeletonTree_20120423_CTRL_1-3.swc',
'SkeletonTree_AEA_20120321_6_7.swc',
'SkeletonTree_AEA_20120327_1-2.swc',
'SkeletonTree_AEA_20120327_5-6_erste Zelle.swc',
'SkeletonTree_AEA_20120327_5-6_zweite Zelle.swc',
'SkeletonTree_AEA_20120327_7-8.swc',
'SkeletonTree_AEA_20120411_7-8-9.swc',
'SkeletonTree_AEA_20120425_5-6.swc',
'SkeletonTree_CTRL_20120315_2-3.swc',
'SkeletonTree_CTRL_20120315_6-8.swc',
'SkeletonTree_CTRL_20120315_CTRL_10-11.swc',
'SkeletonTree_CTRL_20120321_10-11.swc',
'SkeletonTree_CTRL_20120423_1-3.swc',
'SkeletonTree_CTRL_20120423_4-6.swc',
'SkeletonTree_RIM_20120330_1-2-3.swc',
'SkeletonTree_RIM_20120403_1-4.swc',
'SkeletonTree_RIM_20120405_2-3.swc',
'SkeletonTree_RIM_20120405_6-7.swc',
'SkeletonTree_RIM_20120425_3-4.swc',
'SkeletonTree_RIM_20120503_1-2.swc',
'SkeletonTree_RIM_20120503_3-4.swc',
]:
        try: 
            print swc
            morphology   = Morphology.swc_parse( dir+'/'+swc )
            #TODO: print morphology
            morphology.save()
            print morphology.uuid
            keys.append( str(morphology.uuid) )
            morphology   = Morphology.load( morphology.uuid )
            morphologies.append(morphology)
            #morphology.plot_endpoints_histogramm(picture_file='/tmp/amira_endpoints_'+str(morphology.name), picture_formats=['png'])#
        except Exception, e:
            import traceback
            print traceback.format_exc()
            print "%s has no swc format" % (swc)
    print keys
    #morphology.plot_all_properties(morphologies=morphologies, picture_file='/tmp/amira_')

def results(db_ids=[]):
    morphologies    = []
    for uuid in db_ids:
        morphology   = Morphology.load( uuid )
        morphology.info.variable_table(["path_length"])
        morphologies.append(morphology)
        #morphology.plot_endpoints_histogramm(picture_file='/tmp/amira_endpoints_'+str(morphology.name), picture_formats=['png'])#
    morphology.plot_all_properties(morphologies=morphologies, picture_file='/tmp/amira_', picture_formats=['png'])

if __name__ == '__main__':
    #insert('/home/stransky/git/mitsubachi/data/amira/')
    results(['a8177d33-343e-4c97-8eb6-634576e027b7', 'bd13ee3c-e2a4-4203-876d-4b1037e4dd63', 'a13be5c4-6e06-47ad-a3db-99683630a2b9', 'a0f78ee7-e078-48bb-a59e-4849502e2dc8', '4eb7f239-837e-4608-8800-7d61903e91b1', 'c8045bd2-2f19-4d94-ab33-c12a0e0d863e', '5149fe07-7a90-423a-a34e-446333b233ba', 'e03e3be9-b5a7-4001-a19e-8b978502d3a6', '2457c82f-eab2-485b-b251-b87a9ac636b4', '3809ec26-0877-4fcb-ae4d-e401c7a456b3', '27484726-ff0e-4338-8b45-f9eef567c480', 'f6550710-56f7-4005-a9e9-0da5df7e5cc3', 'c53a9508-4039-4087-b140-d084a737d95b', 'd367b616-9931-447e-b8ee-823f696a0a85', '227fb227-d95f-42e9-b279-fa049ff9be2b', '4a4be882-87da-4d62-a78e-7ff51a472f0b', '9f566500-7d4c-458c-a4f2-9a2aabd69b92', '0ff79e12-798f-4217-b7b0-023cedf5b5ef', '526a58fb-8d1c-4465-912b-e58b819d37c9', '89b18884-a8b9-4e0c-90f3-64561207b773', '2617e1f4-4c8e-4049-87fd-c1b3ad98326b'])
