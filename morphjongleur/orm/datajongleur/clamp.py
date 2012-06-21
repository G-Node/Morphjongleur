 # -*- coding: utf-8 -*-
'''
@author: stransky
'''
import sqlalchemy.orm
from datajongleur import Base
from datajongleur.beanbags.models import Identity
from datajongleur.beanbags.models import PREFIX as BB_PREFIX
from datajongleur.beanbags.models import *
import morphjongleur.model.clamp
import morphjongleur.util.auto_string
from morphjongleur.orm.datajongleur.neuron_single_cell_experiment import\
        PointProcess

PREFIX = 'mrj_'


class VClamp(morphjongleur.model.clamp.VClamp, PointProcess):
    __tablename__   = PREFIX + 'v_clamps'
    uuid = sa.Column(
        sqlalchemy.ForeignKey(PointProcess.uuid),
        primary_key=True)


class IClamp(morphjongleur.model.clamp.IClamp, PointProcess):
    __tablename__   = PREFIX + 'i_clamps'
    uuid      = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PointProcess.uuid),
        primary_key=True)
    amplitude       = sqlalchemy.Column('amplitude', sqlalchemy.Float)
    delay           = sqlalchemy.Column('delay', sqlalchemy.Float)
    duration        = sqlalchemy.Column('duration', sqlalchemy.Float)
