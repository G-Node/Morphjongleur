# -*- coding: utf-8 -*-
'''
@author: stransky
@see neo.models
'''
#http://www.sqlalchemy.org/docs/orm/mapper_config.html
import sqlalchemy as sa
import sqlalchemy.orm
from sqlalchemy.orm import backref
from datajongleur import Base
from datajongleur.beanbags.models import Identity
import morphjongleur.model.morphology

PREFIX = 'mrj_'

morphologies_morphology_groups_maps = sqlalchemy.Table(
    PREFIX + 'morphologies_morphology_groups_maps',
    Base.metadata,
    sqlalchemy.Column('morphology_uuid',
        sqlalchemy.ForeignKey(PREFIX + 'morphologies.uuid'),
        primary_key=True),
    sqlalchemy.Column('morphology_group_uuid',
        sqlalchemy.ForeignKey(
          PREFIX + 'morphology_groups.uuid'),
        primary_key=True))

compartments_compartment_groups_maps = sqlalchemy.Table(
    PREFIX + 'compartments_compartment_groups_maps',
    Base.metadata,
    sqlalchemy.Column('compartment_uuid',
        sqlalchemy.ForeignKey(PREFIX + 'compartments.uuid'),
        primary_key=True),
    sqlalchemy.Column('compartment_group_uuid',
        sqlalchemy.ForeignKey(
          PREFIX + 'compartment_groups.uuid'),
        primary_key=True))


class Compartment(morphjongleur.model.morphology.Compartment, Identity):
    __tablename__ = PREFIX + 'compartments'
    __mapper_args__ = {'polymorphic_identity': 'Compartment'}
    uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    morphology_uuid = sqlalchemy.Column(
        sqlalchemy.ForeignKey(PREFIX + 'morphologies.uuid'))
    compartment_id  = sqlalchemy.Column('compartment_id',
        sqlalchemy.Integer)
    compartment_parent_id   = sqlalchemy.Column('compartment_parent_id',
        sqlalchemy.Integer)
    radius          = sqlalchemy.Column('radius', sqlalchemy.Float)
    x   = sqlalchemy.Column('x', sqlalchemy.Float)
    y   = sqlalchemy.Column('y', sqlalchemy.Float)
    z   = sqlalchemy.Column('z', sqlalchemy.Float)


class CompartmentInfo(Identity):
    __tablename__   = PREFIX + 'compartment_infos'
    uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True)
    compartment_uuid = sa.Column(
        sqlalchemy.ForeignKey(PREFIX + 'compartments.uuid'))
    __mapper_args__ = {
            'polymorphic_identity': 'CompartmentInfo'}
   

class CompartmentInfoView(CompartmentInfo):
    """
 parent_radius          = %f,
 length                 = %f,
 cylindric_volume       = %f,
 cylindric_lateral_area = %f,
   frustum_length       = %f,
   frustum_volume       = %f,
   frustum_lateral_area = %f,
 #children              = %i
    """
    __tablename__ = PREFIX + 'v_compartments'
    __mapper_args__ = {'polymorphic_identity': 'CompartmentInfoView'}
    uuid = sa.Column(
        sqlalchemy.ForeignKey(PREFIX + 'compartment_infos.uuid'),
        primary_key=True
    )
    compartment_id  = sqlalchemy.Column('compartment_id',
        sqlalchemy.Integer)
    compartment_parent_id   = sqlalchemy.Column('compartment_parent_id',
        sqlalchemy.Integer) 
    radius          = sqlalchemy.Column('radius', sqlalchemy.Float)
    x               = sqlalchemy.Column('x', sqlalchemy.Float)
    y               = sqlalchemy.Column('y', sqlalchemy.Float)
    z               = sqlalchemy.Column('z', sqlalchemy.Float)
    parent_radius   = sqlalchemy.Column('parent_radius', sqlalchemy.Float)
    parent_x        = sqlalchemy.Column('parent_x', sqlalchemy.Float)
    parent_y        = sqlalchemy.Column('parent_y', sqlalchemy.Float)
    parent_z        = sqlalchemy.Column('parent_z', sqlalchemy.Float)
    length          = sqlalchemy.Column('length', sqlalchemy.Float)
    frustum_length  = sqlalchemy.Column('frustum_length', sqlalchemy.Float)
    frustum_volume  = sqlalchemy.Column('frustum_volume', sqlalchemy.Float)
    frustum_lateral_area    = sqlalchemy.Column('frustum_lateral_area',
        sqlalchemy.Float) 
    cylindric_volume        = sqlalchemy.Column('cylindric_volume',
        sqlalchemy.Float)
    cylindric_lateral_area  = sqlalchemy.Column('cylindric_lateral_area',
        sqlalchemy.Float)
    children        = sqlalchemy.Column('children', sqlalchemy.Integer)


class CompartmentGroups(Identity):
    __tablename__ = PREFIX + 'compartment_groups'
    __mapper_args__ = {'polymorphic_identity': 'CompartmentGroups'}
    uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True
    )
    name    = sqlalchemy.Column('name', sqlalchemy.String)
    description = sqlalchemy.Column('description', sqlalchemy.String)
    type    = sqlalchemy.Column('type', sqlalchemy.String)
    
    _compartment  =   sqlalchemy.orm.relation(Compartment,
        secondary=compartments_compartment_groups_maps,
        backref='_groups')


class Morphology(morphjongleur.model.morphology.Morphology, Identity):
    _type_compartment   = Compartment
    __tablename__ = PREFIX + 'morphologies'
    __mapper_args__ = {
            'polymorphic_identity': 'Morphology',
            #'inherit_condition': uuid==Identity.uuid,
            #'primaryjoin': uuid==Identity.uuid,
            }
    uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True
    )
    name        = sqlalchemy.Column('name', sqlalchemy.String)
    file_origin = sqlalchemy.Column('file_origin', sqlalchemy.String)
    description = sqlalchemy.Column('description', sqlalchemy.String)
    datetime_recording  = sqlalchemy.Column('datetime_recording',
        sqlalchemy.String)
    compartments = sqlalchemy.orm.relation(Compartment,
        backref='morphology',
        primaryjoin='Morphology.uuid==Compartment.morphology_uuid'
        )


class MorphologyInfo(Identity):
    __tablename__   = PREFIX + 'morphology_infos'
    uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True
    )
    morphology_uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        sqlalchemy.ForeignKey(PREFIX + 'morphologies.uuid'))
    __mapper_args__ = {
            'polymorphic_identity': 'MorphologyInfo',
            'inherit_condition': uuid==Identity.uuid}
    morphology             = sqlalchemy.orm.relation(Morphology,
        backref='morphology_infos',
        primaryjoin='Morphology.uuid==MorphologyInfo.morphology_uuid')


class MorphologyInfoView(MorphologyInfo):
    """
    creates new table.
    if view is present, it is used.
    """
    __tablename__ = PREFIX + 'v_morphologies'
    __mapper_args__ = {'polymorphic_identity': 'MorphologyInfoView'}
    uuid = sa.Column(
        sqlalchemy.ForeignKey(PREFIX + 'morphology_infos.uuid'),
        primary_key=True)
    name            = sqlalchemy.Column('name', sqlalchemy.String)
    file_origin     = sqlalchemy.Column('file_origin', sqlalchemy.String)
    description     = sqlalchemy.Column('description', sqlalchemy.String)
    datetime_insert = sqlalchemy.Column('datetime_insert',
        sqlalchemy.String)
    datetime_recording  = sqlalchemy.Column('datetime_recording',
        sqlalchemy.String)
    path_length     = sqlalchemy.Column('path_length', sqlalchemy.Float)
    surface_length  = sqlalchemy.Column('surface_length',
        sqlalchemy.Float) 
    cylindric_volume    = sqlalchemy.Column('cylindric_volume',
        sqlalchemy.Float)
    frustum_volume  = sqlalchemy.Column('frustum_volume', sqlalchemy.Float)
    cylindric_lateral_area  = sqlalchemy.Column('cylindric_lateral_area',
        sqlalchemy.Float)
    frustum_lateral_area    = sqlalchemy.Column('frustum_lateral_area',
        sqlalchemy.Float)
    cylindric_surface_area  = sqlalchemy.Column('cylindric_surface_area',
        sqlalchemy.Float)
    frustum_surface_area    = sqlalchemy.Column('frustum_surface_area',
        sqlalchemy.Float)
    cylindric_mcse          = sqlalchemy.Column('cylindric_mcse',
        sqlalchemy.Float)
    frustum_mcse            = sqlalchemy.Column('frustum_mcse',
        sqlalchemy.Float)
    compartments            = sqlalchemy.Column('compartments',
        sqlalchemy.Integer)
    leafs                   = sqlalchemy.Column('leafs',
        sqlalchemy.Integer)
    branches                = sqlalchemy.Column('branches',
        sqlalchemy.Integer)
    age                     = sqlalchemy.Column('age',
        sqlalchemy.Integer)
    axon                    = sqlalchemy.Column('axon', sqlalchemy.String)


@morphjongleur.util.auto_string.auto_string
class MorphologyGroups(Identity):
    __tablename__ = PREFIX + 'morphology_groups'
    __mapper_args__ = {'polymorphic_identity': 'MorphologyGroups'}
    uuid = sa.Column(
        sqlalchemy.ForeignKey(Identity.uuid),
        primary_key=True
    )
    name        = sqlalchemy.Column('name', sqlalchemy.String)
    description = sqlalchemy.Column('description', sqlalchemy.String)
    age         = sqlalchemy.Column('age', sqlalchemy.Float)

    _morphology = sqlalchemy.orm.relation(Morphology,
        secondary=morphologies_morphology_groups_maps,
        backref='_groups')
