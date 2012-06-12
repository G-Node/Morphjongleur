# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import sqlalchemy.orm
from sqlalchemy.orm import relationship
from datajongleur import Base
from datajongleur.beanbags.models import *
import morphjongleur.model.experiment
from datajongleur.beanbags.models import Identity, RegularlySampledSignal
from datajongleur.utils.sa import NumpyType
from morphjongleur.orm.datajongleur.morphology import Morphology, Compartment
from morphjongleur.orm.datajongleur.neuron import NeuronPassiveParameter

PREFIX = 'mrj_'

point_process_point_process_group_maps = sqlalchemy.Table(
    PREFIX + 'point_process_point_process_group_maps',
    Base.metadata,
    sqlalchemy.Column('point_process_uuid',
        sqlalchemy.ForeignKey(PREFIX + 'point_processes.uuid'),
        primary_key=True),
    sqlalchemy.Column('point_process_group_uuid',
        sqlalchemy.ForeignKey(
          PREFIX + 'point_process_groups.uuid'),
        primary_key=True))


class Experiment(morphjongleur.model.experiment.Experiment, Identity):
    __tablename__   = PREFIX + 'experiments'
    uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    morphology_uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'morphologies.uuid'))
    neuron_passive_parameter_uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(
          PREFIX + 'neuron_passive_parameters.uuid'))
    nseg            = sqlalchemy.Column('nseg', sqlalchemy.Integer)
    description     = sqlalchemy.Column('description', sqlalchemy.String)


class PointProcess(Identity):
    __tablename__   = PREFIX + 'point_processes'
    __mapper_args__ = {'polymorphic_identity': 'PointProcess'}
    uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    experiment_uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'experiments.uuid'))
    compartment_uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'compartments.uuid'))
    position        = sqlalchemy.Column(sqlalchemy.Float)

    compartment = relationship(Compartment,
        backref='point_processes',
        primaryjoin='PointProcess.compartment_uuid==Compartment.uuid')
    experiment = relationship(Experiment,
        backref='point_processes',
        primaryjoin='PointProcess.experiment_uuid==Experiment.uuid')

    @property
    def morphology(self):
        return self.compartment.morphology


class PointProcessGroup(Identity):
    __tablename__   = PREFIX + 'point_process_groups'
    uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    name            = sqlalchemy.Column(sqlalchemy.String)
    description     = sqlalchemy.Column(sqlalchemy.String)
    
    point_processes = relationship(PointProcess, 
            backref='point_process_groups',
            secondary=point_process_point_process_group_maps)


class RecordingPoint(morphjongleur.model.experiment.RecordingPoint,
        PointProcess):
    __tablename__   = PREFIX + 'recording_points'
    __mapper_args__ = {'polymorphic_identity': 'RecordingPoint'}
    uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PointProcess.uuid),
        primary_key=True)

class VoltageTrace(RegularlySampledSignal):
    __tablename__ = PREFIX + 'voltage_traces'
    __mapper_args__ = {'polymorphic_identity': __tablename__}
    uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(RegularlySampledSignal.uuid),
        primary_key=True)
    recording_point_uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'recording_points.uuid'))

    recording_point = relationship(RecordingPoint,
        backref='voltage_traces',
        primaryjoin='VoltageTrace.recording_point_uuid==RecordingPoint.uuid')
    

class VoltageTraceInfo(Identity):
    __tablename__   = PREFIX + 'voltage_trace_infos'
    uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True
    )
    voltage_trace_uuid = sa.Column(
        sqlalchemy.ForeignKey(RegularlySampledSignal.uuid),
        sqlalchemy.ForeignKey(PREFIX + 'voltage_traces.uuid'))
    __mapper_args__ = {
            'polymorphic_identity': 'VoltageTraceInfo',
            'inherit_condition': uuid==Identity.uuid}
    voltage_trace = sqlalchemy.orm.relationship(VoltageTrace,
        backref='voltage_trace_infos',
        primaryjoin='VoltageTrace.uuid==VoltageTraceInfo.voltage_trace_uuid')


class VoltageTraceInfoTauFit(VoltageTraceInfo):
    """
    ---> mitsubachi
    Direkter Zugriff auf diese Attribute von 
    """
    __tablename__ = PREFIX + 'voltage_trace_infos_sinus'
    __mapper_args__ = {'polymorphic_identity': __tablename__}
    uuid = sa.Column(
        sqlalchemy.ForeignKey(VoltageTraceInfo.uuid),
        primary_key=True)
    r_in            = sqlalchemy.Column(sqlalchemy.Float)
    tau_eff         = sqlalchemy.Column(sqlalchemy.Float)
    tau_eff_fit     = sqlalchemy.Column(sqlalchemy.Float)


class VoltageTraceInfoMinMax(VoltageTraceInfo):
    """
    ---> mitsubachi
    """
    __tablename__ = PREFIX + 'voltage_trace_infos_min_max'
    __mapper_args__ = {'polymorphic_identity': __tablename__}
    uuid = sa.Column(
        sqlalchemy.ForeignKey(VoltageTraceInfo.uuid),
        primary_key=True)
    t_min           = sqlalchemy.Column(sqlalchemy.Float)
    v_min           = sqlalchemy.Column(sqlalchemy.Float)
    t_max           = sqlalchemy.Column(sqlalchemy.Float)
    v_max           = sqlalchemy.Column(sqlalchemy.Float)


"""
class VoltageTraceAlt(morphjongleur.model.experiment.VoltageTrace, Identity):
    __tablename__   = PREFIX + 'voltage_traces_alt'
    __mapper_args__ = {'polymorphic_identity': 'VoltageTrace'}
    uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    recording_point_uuid  = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'recording_points.uuid'))
    t               = sqlalchemy.Column('t',NumpyType)
    v               = sqlalchemy.Column('v',NumpyType)
    t_min           = sqlalchemy.Column('t_min',
        sqlalchemy.Float)
    v_min           = sqlalchemy.Column('v_min',
        sqlalchemy.Float)
    t_max           = sqlalchemy.Column('t_max', sqlalchemy.Float)
    v_max           = sqlalchemy.Column('v_max', sqlalchemy.Float)
    
    recording_point = relationship(RecordingPoint,
        backref='voltage_traces_alt',
        primaryjoin='VoltageTraceAlt.recording_point_uuid==RecordingPoint.uuid')


class TauFit(morphjongleur.model.experiment.TauFit, Identity):
    __tablename__   = PREFIX + 'tau_fits'
    __mapper_args__ = {'polymorphic_identity': 'TauFit'}
    uuid     = sqlalchemy.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    voltagetrace_uuid    = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'voltage_traces.uuid'))
    r_in            = sqlalchemy.Column('r_in', sqlalchemy.Float)
    tau_eff         = sqlalchemy.Column('tau_eff', sqlalchemy.Float)
    tau_eff_fit     = sqlalchemy.Column('tau_eff_fit', sqlalchemy.Float)
"""
