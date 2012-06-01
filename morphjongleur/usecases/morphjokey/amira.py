# -*- coding: utf-8 -*-
'''
@author: stransky
'''

def insert():
    import morphjongleur.util.parser.swc
    import morphjongleur.model.morphology
    import morphjongleur.orm.morphology
    db = Database(
                db_name='postgresql://hal08.g-node.pri/morphjongleur',
                exec_role='morphjokey_admin',#TODO: select does not find table!
                exec_path='mitsubachi'
            )
    # must be mapped before Object is created
    mapping = morphjongleur.orm.morphology.Mapper( db )
    mapping.orm_map()
    
    morphologies    = []
    keys            = []
    for swc in [
'amira/Skeletontree_20120315_6-8.swc',
'amira/SkeletonTree_20120330_RIM_1-2-3.swc',
'amira/SkeletonTree_20120423_CTRL_1-3.swc',
'amira/SkeletonTree_AEA_20120321_6_7.swc',
'amira/SkeletonTree_AEA_20120327_1-2.swc',
'amira/SkeletonTree_AEA_20120327_5-6_erste Zelle.swc',
'amira/SkeletonTree_AEA_20120327_5-6_zweite Zelle.swc',
'amira/SkeletonTree_AEA_20120327_7-8.swc',
'amira/SkeletonTree_AEA_20120411_7-8-9.swc',
'amira/SkeletonTree_AEA_20120425_5-6.swc',
'amira/SkeletonTree_CTRL_20120315_2-3.swc',
'amira/SkeletonTree_CTRL_20120315_6-8.swc',
'amira/SkeletonTree_CTRL_20120315_CTRL_10-11.swc',
'amira/SkeletonTree_CTRL_20120321_10-11.swc',
'amira/SkeletonTree_CTRL_20120423_1-3.swc',
'amira/SkeletonTree_CTRL_20120423_4-6.swc',
'amira/SkeletonTree_RIM_20120330_1-2-3.swc',
'amira/SkeletonTree_RIM_20120403_1-4.swc',
'amira/SkeletonTree_RIM_20120405_2-3.swc',
'amira/SkeletonTree_RIM_20120425_3-4.swc',
'amira/SkeletonTree_RIM_20120503_1-2.swc',
'amira/SkeletonTree_RIM_20120503_3-4.swc'
]:
        try: 
            print swc
            morphology   = morphjongleur.model.morphology.Morphology.swc_parse(swc)
            db.store(morphology)
            print morphology.morphology_key
            keys.append(morphology.morphology_key)
            morphology   = mapping.load_morphology( morphology.morphology_key )
            morphologies.append(morphology)
            morphology.plot_endpoints_histogramm(picture_file='/tmp/amira_endpoints_'+str(morphology.name), picture_formats=['png'])#
        except Exception, e:
            import traceback
            print traceback.format_exc()
            print "%s has no swc format" % (swc)
    print keys
    morphology.plot_all_properties(morphologies=morphologies, picture_file='/tmp/amira_')

def results():
    morphologies    = []
    from morphjongleur.io.database import Database
    import morphjongleur.orm.morphology
    db = Database(
                db_name='postgresql://hal08.g-node.pri/morphjongleur',
                exec_role='morphjokey_admin',#TODO: select does not find table!
                exec_path='mitsubachi'
            )
    # must be mapped before Object is created
    mapping = morphjongleur.orm.morphology.Mapper( db )
    mapping.orm_map()
    for i in [47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]:
        morphology   = mapping.load_morphology(i)
        morphology.info.variable_table(["path_length"])
        morphologies.append(morphology)
        #morphology.plot_endpoints_histogramm(picture_file='/tmp/amira_endpoints_'+str(morphology.name), picture_formats=['png'])#
    morphology.plot_all_properties(morphologies=morphologies, picture_file='/tmp/amira_', picture_formats=['png'])

if __name__ == '__main__':
    #insert()
    results()
