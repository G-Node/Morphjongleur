# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import sqlalchemy.orm #from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey
from datajongleur import Base
from datajongleur.beanbags.models import PREFIX as BB_PREFIX
from datajongleur.beanbags.models import *
import morphjongleur.model.neuron_passive
from datajongleur.beanbags.models import Identity

PREFIX = 'mrj_'

class NeuronPassiveParameter(
    morphjongleur.model.neuron_passive.Neuron_passive_parameter,
    Identity): 
    __tablename__   = PREFIX + 'neuron_passive_parameters'
    __mapper_args__ = {'polymorphic_identity': 'NeuronPassiveParameter'}
    uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    Ra              = sqlalchemy.Column('Ra', sqlalchemy.Float)
    g               = sqlalchemy.Column('g', sqlalchemy.Float)
    e               = sqlalchemy.Column('e', sqlalchemy.Float)
