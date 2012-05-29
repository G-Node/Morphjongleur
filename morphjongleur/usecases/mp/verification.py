# -*- coding: utf-8 -*-
'''
@author: stransky
'''

if __name__ == "__main__":
    ###DOES NOT WORK IN SAME CLASS
    #test: 248 , smallest: 392, smallest w/ : 394 , biggest: 433, oldest: 434 
    from mrj.io.database import Database
    
    db = Database(
        db_name='postgresql://hal08.g-node.pri/lehmann',
        exec_role='morphjokey_admin',
        exec_path='morphjokey'
    )

    # must be mapped before Object is created
    import mp.orm.experiment_mso
    mapping = mp.orm.experiment_mso.Mapper( db.engine )
    mapping.orm_map()

    result = db.load_result_mso( 13 )
    #result.experiment_mso.e=-70
    #result.experiment_mso.e=0
    result.plot()#picture_file='../../../doc/voltage_trace_r-0'
    print result.experiment_mso    # takes most of the time
