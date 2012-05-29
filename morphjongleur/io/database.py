# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import sqlalchemy.orm

class Database(object):
    """
    This class handles the connection to a database. 
    It contains the following class variables:
    
    * ``engine``
    * ``session``
    * ``metadata``
    """
    
    engine  = None
    session = None
    metadata = None
    
    '''
    classdocs
    '''
    def __init__(self, #TODO: db_name='sqlite:///morphjokey.db'
            db_name     = 'postgresql://hal08.g-node.pri/morphjongleur',
            exec_role   = 'morphjokey_select',
            exec_path   = 'mitsubachi',
            echo        = False,
            autoflush   = False
        ):
        '''
        Constructor
        exec_role    Benutzerrolle
        exec_path    der Standardpfad
        too verbose: echo=True
        policy? default local mysql or public available g-node database
        TODO: mysql compatible with SET ROLE ?
        '''
        #http://www.sqlalchemy.org/docs/05/dbengine.html
        self.engine = sqlalchemy.create_engine(db_name, echo=echo)
        Session = sqlalchemy.orm.sessionmaker(bind=self.engine, autoflush=autoflush)
        self.session = Session()
        if self.engine  == None:
            self.metadata = sqlalchemy.MetaData()
        else:
            self.metadata = sqlalchemy.MetaData(bind=self.engine)
        self.session.execute("SET ROLE %s;" %(exec_role))
        self.session.execute("SET search_path TO %s;" %(exec_path))

    def store(self, o):
        '''
        committed
        '''
        self.session.add( o )
        self.session.commit()
    
    def create_subtree(self, m):
        '''
        TODO: replace by trigger <- fill compartments before morphology 
        '''
        self.session.execute("SELECT morphjokey.delay_create_node_information(%i)" % (m.morphology_key))
        self.session.commit()
