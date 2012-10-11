# -*- coding: utf-8 -*-
'''
@author: stransky
'''
import morphjongleur.model.morphology

def pca(self):
    import numpy
    import mdp
    m = morphjongleur.model.morphology.Morphology(
      name                  = self.name,
      file_origin           = self.file_origin,
      description           = self.description,
      datetime_recording    = self.datetime_recording
    )
    assert m.number_of_compartments == 0
    pca_cs  = mdp.pca( numpy.array([ [c.x, c.y, c.z] for c in self.compartments ] ) )
    assert self.number_of_compartments == len(pca_cs)
    for i in xrange( self.number_of_compartments ):
        m.add_compartment(
            morphjongleur.model.morphology.Compartment( 
                self._compartments[i].compartment_id, 
                self._compartments[i].compartment_parent_id, 
                self._compartments[i].radius, 
                x=pca_cs[i][0],
                y=pca_cs[i][1],
                z=pca_cs[i][2]
            ) 
        )
    assert self.number_of_compartments == m.number_of_compartments
    return m

def scale(self, scale_factor):
    m = morphjongleur.model.morphology.Morphology(
      name                  = self.name,
      file_origin           = self.file_origin,
      description           = self.description,
      datetime_recording    = self.datetime_recording
    )
    assert m.number_of_compartments == 0
    for i in xrange( self.number_of_compartments ):
        c   = morphjongleur.model.morphology.Compartment( 
                self._compartments[i].compartment_id, 
                self._compartments[i].compartment_parent_id, 
                self._compartments[i].radius
            )
        c._length   = scale_factor * self._compartments[i].length
        m.add_compartment( c )
    assert self.number_of_compartments == m.number_of_compartments
    return m

morphjongleur.model.morphology.Morphology.pca   = pca
morphjongleur.model.morphology.Morphology.scale = scale

if __name__ == '__main__':
    import morphjongleur.util.parser.swc
    import morphjongleur.util.metric_analysis

    m = morphjongleur.model.morphology.Morphology.swc_parse('../../data/test.swc')
    ma   = morphjongleur.util.metric_analysis.MetricAnalysis(m)
    (ks, vs)    = ma.variable_table(['name', 'compartments', #'datetime_recording', 
        'total_length', 'slant_length', 
        'number_of_branching_points', 'number_of_terminaltips', 
        'cylindric_volume', 'cylindric_surface_area', 'cylindric_compactness', 
        'frustum_volume', 'frustum_surface_area', 'frustum_compactness'
        ])

    print ks
    print vs

    m.pca()

    ms  = m.scale(0.5)
    mas = morphjongleur.util.metric_analysis.MetricAnalysis(ms)
    (ks, vs)    = mas.variable_table(['name', 'compartments', #'datetime_recording', 
        'total_length', 'slant_length', 
        'number_of_branching_points', 'number_of_terminaltips', 
        'cylindric_volume', 'cylindric_surface_area', 'cylindric_compactness', 
        'frustum_volume', 'frustum_surface_area', 'frustum_compactness'
        ])

    print vs
