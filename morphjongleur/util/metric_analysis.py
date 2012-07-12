#!env python
# -*- coding: utf-8 -*-
'''
frustrum v.s. cylinder

@author: stransky
'''
import math
import numpy
import mdp
from scipy.stats.stats import gmean, hmean#, mean, median, medianscore, mode
from chull import Vector, Hull
from morphjongleur.model.morphology import *

class MetricAnalysis(MorphologyInfo):
    '''
    classdocs
    total cell length = sum of all paths (based on center of each compartment)
    surface_length_frustum    total cell length running over the surface of an frustrum
    
    
    polyeder    #konvexen polyeder zur not emal in mathesoftware
    http://python.net/~gherman/convexhull.html
http://www.scipy.org/Cookbook
http://michelanders.blogspot.de/2012/02/3d-convex-hull-in-python.html
http://code.activestate.com/recipes/66527-finding-the-convex-hull-of-a-set-of-2d-points/

    cylindric_lateral_area
    without side circles of terminal tips


    self.frustum_lateral_area
    without side circles of terminal tips
    
    
    arithmetic_mean_branchpoint_radius
    mean cross section area

    geometric_mean_branchpoint_radius
    mean cross section area

    harmonic_mean_branchpoint_radius    

    arithmetic_mean_branchpoint_distance    

    geometric_mean_branchpoint_distance    

    harmonic_mean_branchpoint_distance    

    cylindric_arithmetic_mean_cross_section_area    

    frustum_arithmetic_mean_cross_section_area    
    
    '''

    @staticmethod
    def euclid_distance(x, y):
        if len(x) != len(y):
            raise Exception, 'different dimension';
        dist = 0;
        for i in range(len(x)):
            dist += (x[i] - y[i])**2;
        return math.sqrt(dist);
    

    @staticmethod
    def inverse_power_means(x,r):
        s = 0
        for x_i in x:
            s += x_i ** r
        return math.pow(s,1./r)
 
    @staticmethod
    def powermeans(X, W = None, r= 1):
        '''
        @see http://www.scipy.org/doc/api_docs/SciPy.stats.stats.html
        @see http://adorio-research.org/wordpress/?p=241
        Computes the mean using the formula
           [(sum x**r)/n]**(1/r).
        A general function, use specific functions
        for more accuracy, example, for r == 1,
        call the mean(X, method= ...). See previous blogs.
           X - original data
           W - optional weight array.
           r - power (usually integer valued.)
        '''
        n = math.len(X)
        if W is None:
            if r == -1: #harmonic mean
                return 1.0 / (sum([ (1.0/x) for x in X])/n)
            elif r == 0: # geometric mean
                return math.exp(sum([math.log(x) for x in X])/n)
            elif r == 1: # arithmetic mean. 
                return sum(X) / float(n)
            elif r == 2: # rms.
                return  math.sqrt(sum([x*x for x in X])/n)
            else:
                return (sum([x**r for x in X])/n)**(1.0/r)
        else:
            if r == -1: #harmonic mean
                return 1.0 / (sum([ w * (1.0/x) for x,w in zip(X,W)])/sum(W))
            elif r == 0: # geometric mean
                return math.exp(sum([w * math.log(x) for (x,w) in zip(X,W)])/sum(W))
            elif r == 1: # arithmetic mean. 
                return sum(w * x for (x,w) in zip(X,W)) / float(sum(W))
            elif r == 2: # rms.
                return  math.sqrt(sum([w * x*x for x,w in zip(X,W)])/sum(W))
            else:
                return (sum([w * x**r for x,w in zip(X,W)])/n)**(1.0/sum(W))
        
    def __init__(self, morphology):
        '''
        '''
        self.name               = morphology.name
        self.file_origin        = morphology.file_origin
        self.description        = morphology.description
        self.datetime_recording = morphology.datetime_recording
        self.compartments   = len(morphology.compartments)
        
        self.leafs  = morphology.leafs
        self.branching_points  = morphology.branching_points
        self.number_of_leafs = len(self.leafs)
        self.number_of_branching_points = len(self.branching_points)
#TODO:        if self.number_of_leafs == self.number_of_branching_points + 1:
#            print >> sys.stderr, "#leafs %i != #branches %i +1" % (self.number_of_leafs, self.number_of_branching_points + 1)
        #assert self.number_of_leafs == self.number_of_branching_points + 1

        self.total_cell_length    = 0.
        self.surface_length_frustrum = 0.
        
        self.cylindric_volume = 0.
        self.cylindric_lateral_area = 0.
        
        self.frustum_volume = 0.
        self.frustum_lateral_area   = 0.

        self.arithmetic_mean_cross_section_area  = 0.
        self.geometric_mean_cross_section_area   = 1.
        self.harmonic_mean_cross_section_area    = 0.
        
        self.arithmetic_mean_branchpoint_cross_section_area = 0.
        self.geometric_mean_branchpoint_cross_section_area  = 1.
        self.harmonic_mean_branchpoint_cross_section_area   = 0.

        self.arithmetic_mean_branchpoint_distance    = 0.
        self.geometric_mean_branchpoint_distance    = 1.
        self.harmonic_mean_branchpoint_distance    = 0.
        
        xs = []
        for compartment in morphology.getNonRootCompartments():
            xs.append([compartment.x, compartment.y, compartment.z])

            self.total_cell_length += compartment.parent_distance
            
            self.cylindric_volume += math.pi* compartment.radius**2 * compartment.parent_distance

            self.cylindric_lateral_area += 2 * math.pi* compartment.radius * compartment.parent_distance

            
            compartment.frustum_length = math.sqrt( (compartment.parent.radius - compartment.radius)**2 + compartment.parent_distance**2)
            self.surface_length_frustrum += compartment.frustum_length
            
            self.frustum_volume  += math.pi/3. * compartment.parent_distance * (compartment.parent.radius**2 + compartment.parent.radius * compartment.radius + compartment.radius**2)            
            
            self.frustum_lateral_area   += math.pi * compartment.frustum_length * (compartment.parent.radius + compartment.radius)
           
            self.arithmetic_mean_cross_section_area   += 2* math.pi * compartment.radius
            self.geometric_mean_cross_section_area    *= 2* math.pi * compartment.radius
            if compartment.radius > 0: 
                self.harmonic_mean_cross_section_area     += 1./(2* math.pi * compartment.radius)
            else:
#TODO:                print >> sys.stderr, "compartment %i in morphology %s has radius %f" % (compartment.compartment_id, morphology.name, compartment.radius)
                self.harmonic_mean_cross_section_area   = float('nan')

        for compartment in self.branching_points:
            cross_section_area  = 2* math.pi * compartment.radius
            self.arithmetic_mean_branchpoint_cross_section_area += cross_section_area
            self.geometric_mean_branchpoint_cross_section_area  *= cross_section_area
            if compartment.radius > 0: 
                self.harmonic_mean_branchpoint_cross_section_area   += 1./(cross_section_area)
            else:
                print >> sys.stderr, "compartment %i in morphology %s has radius %f" % (compartment.compartment_id, morphology.name, compartment.radius)
                self.harmonic_mean_branchpoint_cross_section_area   = float('nan')

            parent_brachpoint_distance = float('nan')
            #TODO: gewichtung nach volumen
            self.arithmetic_mean_branchpoint_distance   += parent_brachpoint_distance
            self.geometric_mean_branchpoint_distance    *= parent_brachpoint_distance
            if parent_brachpoint_distance > 0:
                self.harmonic_mean_branchpoint_distance     += 1./parent_brachpoint_distance
#            else:
#                print >> sys.stderr, "compartment %i in morphology %s has branchpoint_distance %f" % (compartment.compartment_id, morphology.name, parent_brachpoint_distance)
#                self.harmonic_mean_branchpoint_distance   = float('nan')
            
            
        terminaltips_cross_section_area = 0
        for compartment in self.leafs:
            terminaltip_cross_section_area  = 2* math.pi * compartment.radius
            terminaltips_cross_section_area += terminaltip_cross_section_area
        
        self.cylindric_surface_area = self.cylindric_lateral_area + terminaltips_cross_section_area
        self.frustum_surface_area   = self.frustum_lateral_area   + terminaltips_cross_section_area
        
        self.cylindric_mean_cross_section_area   = self.cylindric_volume / self.total_cell_length
        self.frustum_mean_cross_section_area     = self.frustum_volume / self.total_cell_length
        self.arithmetic_mean_cross_section_area  /= self.compartments
        self.geometric_mean_cross_section_area   = math.pow(self.geometric_mean_cross_section_area, 1./self.compartments)
        self.harmonic_mean_cross_section_area    = self.compartments/self.harmonic_mean_cross_section_area

        self.arithmetic_mean_branchpoint_cross_section_area /= self.number_of_branching_points
        self.geometric_mean_branchpoint_cross_section_area  = math.pow(self.geometric_mean_branchpoint_cross_section_area, 1./self.number_of_branching_points)
        self.harmonic_mean_branchpoint_cross_section_area   = self.number_of_branching_points / self.harmonic_mean_branchpoint_cross_section_area
        
        self.mean_branchpoint_distance   = self.total_cell_length / self.number_of_leafs   #not self.number_of_branching_points
        self.arithmetic_mean_branchpoint_distance    /= self.number_of_branching_points
        self.geometric_mean_branchpoint_distance    = math.pow(self.geometric_mean_branchpoint_distance, 1./self.number_of_branching_points)
        self.harmonic_mean_branchpoint_distance    = float('nan')#TODO: self.number_of_branching_points / self.harmonic_mean_branchpoint_distance

        mins   = [float('inf'), float('inf'), float('inf')]
        maxs   = [float('-inf'),float('-inf'),float('-inf')]
        for x in mdp.pca( numpy.array(xs) ):
            for d in range(3):
                if x[d] < mins[d]:
                    mins[d]    = x[d]
                if x[d] > maxs[d]:
                    maxs[d] = x[d]
        self.pca_lengths = (
            maxs[0] - mins[0],
            maxs[1] - mins[1],
            maxs[2] - mins[2]
        )
        self.pca_box= 2*( self.pca_lengths[0]*self.pca_lengths[1]
                        + self.pca_lengths[1]*self.pca_lengths[2]
                        + self.pca_lengths[2]*self.pca_lengths[0]
                        )
        
        self.pca_rhob2=float('nan')

        h    = Hull([Vector.fromArray(x) for x in xs])
        self.convex_enveloping_polyhedron_surface_area   = h.surface_area()
    

if __name__ == '__main__':
    '''
    Parameter: files, not directories

    cd data/amira/
    pathto/metric_analysis.py *.swc
    '''
    import sys
    import morphjongleur.util.parser.swc
    with_head   = True
    for swc in sys.argv[1:]:#['/tmp/mitsubachi/test.swc','/tmp/mitsubachi/H060602DB_10_2_zentai_.swc','/tmp/mitsubachi/H060602VB_10_2_zentai_.swc','/tmp/mitsubachi/H060607DB_10_2(zentai).swc','/tmp/mitsubachi/H060607VB_10_2(zentai).swc']:#['../../data/test.swc']:#
        
        

        
        
        m   = Morphology.swc_parse(swc, verbose=False)
        a   = MetricAnalysis(m)
        (ks, vs)    = a.variable_table(['name', 'datetime_recording', 
        'compartments', 'number_of_branching_points', 
        'total_cell_length', 'surface_length_frustrum', 
        'cylindric_volume', 'cylindric_surface_area', 
        'frustum_volume', 'frustum_surface_area', 
        'cylindric_mean_cross_section_area', 'frustum_mean_cross_section_area', 
        'mean_branchpoint_distance', 
        'pca_lengths', 'convex_enveloping_polyhedron_surface_area'
        ])
        
        if with_head:
            print ks
            with_head   = False
        print vs
        
        #print 80*'_'
        #print str(a)
        #print 80*'_'
        #print repr(a)
        #print 80*'_'
        #print m
        #print 80*'_'
    #a.convex_enveloping_polyeder_hull.Print()