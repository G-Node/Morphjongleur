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

PREFIX = 'mrj_'

clamps_clamp_groups_map = sqlalchemy.Table(
    'clamps_clamp_groups_map',
    Base.metadata,
    sqlalchemy.Column(
        'iclamp_uuid',
        sqlalchemy.ForeignKey(PREFIX + 'iclamps.uuid'),
        primary_key=True),
    sqlalchemy.Column(
        'iclamp_group_key',
        sqlalchemy.ForeignKey(PREFIX + 'iclamps_groups.iclamp_group_key'),
        primary_key=True)

class VClamp(morphjongleur.model.clamp.VClamp, Identity):
    __tablename__   = PREFIX + 'vclamps'
    uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    experiment_uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'experiments.uuid'))
    morphology_uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'morphologies.uuid'))
    compartment_uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'compartments.uuid'))
    position        = sqlalchemy.Column(sqlalchemy.Float)

    @property
    def morphology_uuid(self):
        return compartment

class IClamp(morphjongleur.model.clamp.IClamp, Identity):
    __tablename__   = PREFIX + 'iclamps'
    uuid      = sqlalchemy.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    experiment_key  = sqlalchemy.Column('experiment_key',
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(PREFIX + 'experiments.experiment_key'))
    morphology_key  = sqlalchemy.Column('morphology_key',
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(PREFIX + 'morphologies.morphology_key'))
    compartment_key = sqlalchemy.Column('compartment_key',
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(PREFIX + 'compartments.compartment_key'))
    compartment_id  = sqlalchemy.Column('compartment_id',
        sqlalchemy.Integer)
    position        = sqlalchemy.Column('position', sqlalchemy.Float)
    amplitude       = sqlalchemy.Column('amplitude', sqlalchemy.Float)
    delay           = sqlalchemy.Column('delay', sqlalchemy.Float)
    duration        = sqlalchemy.Column('duration', sqlalchemy.Float)


class PatternClamp(morphjongleur.model.clamp.PatternClamp, Identity):
    __tablename__   = PREFIX + 'vpatternclamps'
    patterngc=   sqlalchemy.Column('patternclamp_key',
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(BB_PREFIX + 'identities.uuid'),
        primary_key=True)
    experiment_key  = sqlalchemy.Column('experiment_key',
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(PREFIX + 'experiments.experiment_key'))
    morphology_key  = sqlalchemy.Column('morphology_key',
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(PREFIX + 'morphologies.morphology_key'))
    compartment_key = sqlalchemy.Column('compartment_key',
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(PREFIX + 'compartments.compartment_key'))
    compartment_id  = sqlalchemy.Column('compartment_id', sqlalchemy.Integer)
    position        = sqlalchemy.Column('position', sqlalchemy.Float)
    delta_t         = sqlalchemy.Column('delta_t', sqlalchemy.Float)
    delay           = sqlalchemy.Column('delay', sqlalchemy.Float)
    duration        = sqlalchemy.Column('duration', sqlalchemy.Float)


class ClampGroups(Identity):
    __tablename__   = PREFIX + 'iclamps_groups'
    uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(BB_PREFIX + 'identities.uuid'),
        primary_key=True)
    name            = sqlalchemy.Column(sqlalchemy.String)
    description     = sqlalchemy.Column(sqlalchemy.String)
    voltagetrace_key= sqlalchemy.Column(
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(PREFIX + 'voltage_traces.voltagetrace_key'))
    morphology_key  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'morphologies.morphology_key'))
    compartment_key = sqlalchemy.Column(
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(PREFIX + 'compartments.compartment_key'))
    compartment_id  = sqlalchemy.Column(sqlalchemy.Integer)
    position        = sqlalchemy.Column(sqlalchemy.Float)
    amplitude       = sqlalchemy.Column(sqlalchemy.Float)
    function        = sqlalchemy.Column(sqlalchemy.String)
    delta_t         = sqlalchemy.Column(sqlalchemy.Float)
    delay           = sqlalchemy.Column(sqlalchemy.Float)
    duration        = sqlalchemy.Column(sqlalchemy.Float)
