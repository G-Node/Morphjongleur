# -*- coding: utf-8 -*-
'''
@author: stransky
'''
from mrj.model.neuron_mso import Neuron_MSO
from mrj.model.experiment_mso import Experiment_MSO

if __name__ == "__main__":
    ###DOES NOT WORK IN SAME CLASS
    #test: 248 , smallest: 392, smallest w/ : 394 , biggest: 433, oldest: 434, fastest: 424, slowest: 393
    from mrj.io.database import Database
    import time
    
    db = Database(
        db_name='postgresql://hal08.g-node.pri/lehmann',
        exec_role='morphjokey_admin',
        exec_path='morphjokey'
    )

    # must be mapped before Object is created
    import mrj.orm.experiment_mso
    mapping = mrj.orm.experiment_mso.Mapper( db.engine )
    mapping.orm_map()

#    mapping.create_tables()

    m = db.load_morphology( 424 )
    #from mrj.io.swc import SwcParser
    #swc_parser  = SwcParser()
    #m   = swc_parser.parse('../../../Data/test.swc')

    n   = Neuron_MSO(m)#, exitatoric_medial_dendrits=[4], exitatoric_lateral_dendrits=[8]
    
    s = Experiment_MSO(n, f=200, r=1.00, n=None, duration=10e-3, delay=1e-3, group_name=None, description='test')#, group_name='medial'
    s.neuron_create()
    #db.store( s )
    #print s

    start_time = time.time()
    s.run_simulation()
    end_time = time.time()
    result  = s.get_result()#Result(t=[0.001,0.002,0.003,0.004],v=[-60.0,-60.0,-60.0,-66.0],t_min=0.03,v_min=-66,experiment=experiment)#

    s.calc_time = end_time-start_time
    print("simulation needed %s seconds." %(s.calc_time))

    #db.store( result )
    #print( result )
    #result.plot()
    result.plot(picture_file='../../../doc/voltage_trace_f-200_r-100')##'#1,0.93,0.64,0
