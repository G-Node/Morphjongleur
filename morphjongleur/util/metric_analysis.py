'''
frustrum v.s. cylinder

@author: stransky
'''
import math
import numpy
import morphjongleur.util.auto_string
from morphjongleur.model.morphology import Morphology,Star

@morphjongleur.util.auto_string.auto_string
class MetricAnalysis(object):
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
            raise Exception, "different dimension";
        dist = 0;
        for i in range(len(x)):
            dist += (x[i] - y[i])**2;
        return math.sqrt(dist);
    
    def __init__(self, morphology):
        '''
        '''
        self.compartments   = len(morphology.compartments)
        
        #TODO: assert len(morphology.leafs) == morphology.branching_points + 1
        self.leafs  = morphology.leafs
        self.branching_points  = morphology.branching_points

        self.total_cell_length    = 0.
        self.surface_length_frustum = 0.
        self.cylindric_volume = 0.
        self.cylindric_lateral_area = 0.
        self.arithmetic_mean_branchpoint_radius = 0.
        self.geometric_mean_branchpoint_radius  = 1.
        xyz = []
        for compartment in morphology.getCompartments():
            xyz.append([compartment.x, compartment.y, compartment.z])
            
            self.total_cell_length += compartment.parent_distance
            
            self.cylindric_volume += math.pi* compartment.radius**2 * compartment.parent_distance

            self.cylindric_lateral_area += 2 * math.pi* compartment.radius * compartment.parent_distance

            
            compartment.frustum_length = math.sqrt( (compartment.parent.radius - compartment.radius)^2 + compartment.parent_distance^2)
            self.surface_length_frustrum += compartment.frustum_length
            
            self.frustum_volume  += math.pi/3. * compartment.parent_distance * (compartment.parent_radius^2 + compartment.parent_radius * compartment.radius + compartment.radius^2)            
            
            self.frustum_lateral_area   += math.pi * compartment.frustum_length * (compartment.parent_radius + compartment.radius)
            
            self.arithmetic_mean_branchpoint_radius += compartment.radius
            self.geometric_mean_branchpoint_radius  *= compartment.radius
            self.harmonic_mean_branchpoint_radius   += float('nan')
            
            
            self.cylindric_arithmetic_mean_cross_section_area   += float('nan')
            self.frustum_arithmetic_mean_cross_section_area     += float('nan')
        
        self.cylindric_surface_area += float('nan')            
        self.frustum_surface_area   += float('nan')
        
        self.arithmetic_mean_branchpoint_distance   = self.total_cell_length / self.branching_points
        self.geometric_mean_branchpoint_distance    += float('nan')
        self.harmonic_mean_branchpoint_distance += float('nan')
            

        import mdp
        x = None;
        for r in result:
            c = xyz#numpy.array([r['x'], r['y'], r['z']]);
            if(x == None):
                x = c;
            else:
                x = numpy.vstack((x,c));        
        p = mdp.pca(x);
        x_min   = float("inf");
        x_max   = float("-inf");
        for x,y,z in p:
            if x < x_min:
                x_min = x; 
            if x > x_max:
                x_max = x;
        self.pca_length = x_max - x_min

        self.polyeder   = float('nan')

if __name__ == '__main__':
    
    import morphjongleur.util.parser.swc
    m   = Morphology.swc_parse('../../data/test.swc')
    a   = MetricAnalysis(m)
    print a
    print m