#!/usr/bin/python
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
        self.terminaltips_distance   = 0
        for compartment in self.leafs:
            terminaltip_cross_section_area  = 2* math.pi * compartment.radius
            terminaltips_cross_section_area += terminaltip_cross_section_area
            self.terminaltips_distance += morphology.biggest.path_distance(compartment)
        self.mean_terminaltips_distance   = self.terminaltips_distance / len(self.leafs)
        
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
        self.pca_length_x  = maxs[0] - mins[0]
        self.pca_length_y  = maxs[1] - mins[1]
        self.pca_length_z  = maxs[2] - mins[2]
        
        self.pca_lengths = (
            self.pca_length_x,
            self.pca_length_y,
            self.pca_length_z
        )
        #self.pca_box_surface_area= 2*( self.pca_lengths[0]*self.pca_lengths[1]
        #                + self.pca_lengths[1]*self.pca_lengths[2]
        #                + self.pca_lengths[2]*self.pca_lengths[0]
        #                )
        '''
        2 * (
            1/2. * self.pca_lengths[0] * numpy.sqrt(numpy.square(self.pca_lengths[1]/2) + numpy.square(self.pca_lengths[2]/2)) 
            +
            1/2. * self.pca_lengths[0] * numpy.sqrt(numpy.square(self.pca_lengths[2]/2) + numpy.square(self.pca_lengths[1]/2)) 
            )
        '''
        #self.pca_rhombus =  self.pca_lengths[0] * numpy.sqrt(numpy.square(self.pca_lengths[2]) + numpy.square(self.pca_lengths[1]))

        h    = Hull([Vector.fromArray(x) for x in xs])
        self.convex_enveloping_polyhedron_surface_area, self.convex_enveloping_polyhedron_volume   = h.surface_area_and_volume()

    @staticmethod
    def plot_all_properties(morphologies=[], picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        MetricAnalysis.plot_property(morphologies, quantity="cylindric_volume",                                                        picture_file=picture_file, picture_formats=picture_formats)
        MetricAnalysis.plot_property(morphologies, quantity="frustum_volume",       yaxis_description=u'volume [µm³]',                 picture_file=picture_file, picture_formats=picture_formats)
        MetricAnalysis.plot_property(morphologies, quantity="cylindric_surface_area",                                                  picture_file=picture_file, picture_formats=picture_formats)
        MetricAnalysis.plot_property(morphologies, quantity="frustum_surface_area", yaxis_description=u'surface area [µm²]',           picture_file=picture_file, picture_formats=picture_formats)
        MetricAnalysis.plot_property(morphologies, quantity="branches",             yaxis_description='#branches',                     picture_file=picture_file, picture_formats=picture_formats)
        MetricAnalysis.plot_property(morphologies, quantity="path_length",          yaxis_description=u'total cell length [µm]',       picture_file=picture_file, picture_formats=picture_formats)
        MetricAnalysis.plot_property(morphologies, quantity="cylindric_mcsa",                                                          picture_file=picture_file, picture_formats=picture_formats)
        MetricAnalysis.plot_property(morphologies, quantity="frustum_mcsa",         yaxis_description=u"mean cross-section area [µm]", picture_file=picture_file, picture_formats=picture_formats)
        #TODO: MetricAnalysis.plot_property(morphologies, quantity="spatial_strech", yaxis_description='spatial_strech [$\mu$m]',               picture_file=picture_file, picture_formats=picture_formats)

    @staticmethod
    def plot_property(morphologies=[], quantity='volumes', yaxis_description=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import matplotlib.pyplot
        import numpy
        #matplotlib.rc('text', usetex=True): error with names

        values   = []
        xticks      = []
        for m in morphologies:
            values.append( m.info.__dict__[quantity] )
            xticks.append( m.name.replace('_',' ') )
        ind = numpy.arange(len(morphologies))    # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence
        
        if yaxis_description == None:
            matplotlib.pyplot.ylabel(quantity.replace('_',' '))
        else:
            matplotlib.pyplot.ylabel(yaxis_description)
        matplotlib.pyplot.title('change: DB: %f, VB: %f' %(float(morphologies[1].info.__dict__[quantity])/morphologies[0].info.__dict__[quantity], float(morphologies[3].info.__dict__[quantity])/morphologies[2].info.__dict__[quantity]))#TODO: verallgemeinern, verifizierren!
        print '%s change: DB: %f, VB: %f' %(quantity, float(morphologies[1].info.__dict__[quantity])/morphologies[0].info.__dict__[quantity], float(morphologies[3].info.__dict__[quantity])/morphologies[2].info.__dict__[quantity])
        matplotlib.pyplot.xticks(ind+width, xticks,rotation=9)#
        #matplotlib.pyplot.yticks(numpy.arange(0,81,10))
        matplotlib.pyplot.grid(True, color='lightgrey')
        
        matplotlib.pyplot.bar(ind, values, width, color='black')#

        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+quantity+'.'+picture_format,format=picture_format)
                except Exception, e:
                    import traceback
                    print picture_format 
                    print traceback.format_exc()
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

    @staticmethod
    def plot_endpoints_histogramm(morphology, xlim=None, ylim=None, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
        import matplotlib.pyplot
        import numpy
        #matplotlib.rc('text', usetex=True): error with names

        x   = []
        for c in morphology.leafs:
            #print "%i/%i\r" % (len(x),len(morphology.leafs)), 
            x.append(c.path_distance(morphology.biggest))
        mean   = numpy.mean(x)
        std    = numpy.std(x)

        #matplotlib.pyplot.title('Endpoints of %s' % (morphology.name.replace('_',' ')) )
        print 'Endpoints of %s : mean=%f, std=%f' % (morphology.name, mean, std)
        
        matplotlib.pyplot.axvline(x=mean, color='black')
        matplotlib.pyplot.grid(True, color='lightgrey')
        if xlim != None:#TODO: isnumber
            matplotlib.pyplot.hist(x, bins=range(xlim[0],xlim[1],(xlim[1]-xlim[0])/100), normed=0, color=color)#, label='my data'
        else:
            matplotlib.pyplot.hist(x, 20, normed=0, color=color)#, label='my data'

        matplotlib.pyplot.ylabel('#')#%
        if ylim != None:
            matplotlib.pyplot.ylim(ylim)
        matplotlib.pyplot.xlabel(u'distance [µm]')
        if xlim != None:
            matplotlib.pyplot.xlim(xlim)

        xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis()
        width    = std / (xmax-xmin)
        center    = (mean - xmin) / (xmax-xmin)
        matplotlib.pyplot.axhline(y=.5*(ymax-ymin), xmin=center-.5*width, xmax=center+.5*width, color='red')#TODO: höhe der Line = mittelwert der bins
        matplotlib.pyplot.legend( [ 'mean %f' % ( mean ), 'std %f' % ( std )  ] )

        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

if __name__ == '__main__':
    '''
    Parameter: files, not directories

    cd data/swc_files/
    path_to/metric_analysis.py *.swc
    '''
    import sys
    import morphjongleur.util.parser.swc
    with_head   = True
    for swc in ['/home/stransky/git/mitsubachi/data/mitsubachi/H060607VB_10_2(whole).swc']:#sys.argv[1:]:#['../../data/test.swc']:#
# H060602DB_10_2(whole).swc  H060607DB_10_2(whole).swc  H060602VB_10_2(whole).swc  H060607VB_10_2(whole).swc
# #00ff00 #008000 #00ff80 #008080
        morphology   = Morphology.swc_parse(swc, verbose=False)
        MetricAnalysis.plot_endpoints_histogramm(morphology, xlim=(0, 1000), ylim=(0, 50), color='#008080', picture_file='/tmp/endpointdistribution_'+str(morphology.name), picture_formats=['svg', 'png'])
        a   = MetricAnalysis(morphology)
        (ks, vs)    = a.variable_table(['name', #'datetime_recording', 
        'compartments', 'number_of_branching_points', 
        'total_cell_length', 'surface_length_frustrum', 
        'cylindric_volume', 'cylindric_surface_area', 
        'frustum_volume', 'frustum_surface_area', 
        'pca_length_x', 'pca_length_y', 'pca_length_z', 
        'convex_enveloping_polyhedron_volume', 'convex_enveloping_polyhedron_surface_area', 
        ])
        
        if with_head:
            print ks
            with_head   = False
        print vs
        
        morphology.plot(picture_file='/tmp/%s' % (morphology.name), picture_formats=['svg', 'png'])
        m_pca   = morphology.pca()
        m_pca.swc_write('/tmp/%s_pca.swc' % (m_pca.name) )
        m_pca.plot(picture_file='/tmp/%s_pca' % (m_pca.name), picture_formats=['svg', 'png'])
        #print 80*'_'
        #print str(a)
        #print 80*'_'
        #print repr(a)
        #print 80*'_'
        #print m
        #print 80*'_'
    #a.convex_enveloping_polyeder_hull.Print()
