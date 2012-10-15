#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
frustum v.s. cylinder

@author: stransky
'''
import math
import numpy
import mdp
import scipy.stats.stats #import gmean, hmean#, mean, median, medianscore, mode
import matplotlib.pyplot
import numpy
import morphjongleur.util.chull
from morphjongleur.model.morphology import MorphologyInfo,Morphology,Compartment
import traceback

class MetricAnalysis(MorphologyInfo):
    '''
    classdocs
    total cell length = sum of all paths (based on center of each compartment)
    slant_length    total cell length running over the surface of an frustum
    
    
    polyeder    #konvexen polyeder zur not emal in mathesoftware
    http://python.net/~gherman/convexhull.html
http://www.scipy.org/Cookbook
http://michelanders.blogspot.de/2012/02/3d-convex-hull-in-python.html
http://code.activestate.com/recipes/66527-finding-the-convex-hull-of-a-set-of-2d-points/

    cylindric_lateral_area
    without side circles of terminal tips


    self.frustum_lateral_area
    without side circles of terminal tips
    
    
    am_branch_point_radius
    mean cross section area

    gm_branch_point_radius
    mean cross section area

    hm_branch_point_radius    

    am_branch_point_distance    

    gm_branch_point_distance    

    hm_branch_point_distance    

    cylindric_am_cross_section_area    

    frustum_am_cross_section_area    

    http://openbook.galileocomputing.de/python/python_kapitel_13_009.htm
    TODO: glossary
    path_length         = %f, 
 surface_length         = %f, 
 cylindric_volume       = %f, 
   frustum_volume       = %f, 
 cylindric_lateral_area = %f, 
   frustum_lateral_area = %f, 
 cylindric_surface_area = %f, 
   frustum_surface_area = %f, 
 #branches              = %i
    ''' 

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
        self.morphology         = morphology

        #redundancy for db
        self.name               = morphology.name
        self.file_origin        = morphology.file_origin
        self.description        = morphology.description
        self.datetime_recording = morphology.datetime_recording
        self.compartments       = morphology.number_of_compartments
        
        self.number_of_terminal_tips     = morphology.number_of_terminal_tips
        self.number_of_branch_points = morphology.number_of_branch_points
#        if self.number_of_terminal_tips != self.number_of_branch_points + 1:
#            print >> sys.stderr, "#terminal_tips %i > #branching points %i +1, because of %s" % (self.number_of_terminal_tips, self.number_of_branch_points, ["%i:%i"%(pleb.compartment_id,len(pleb.children)) for pleb in morphology.plebs])

        terminal_tips_radii      = [compartment.radius for compartment in morphology.terminal_tips]
        terminal_tips_cross_section_area = 2 * math.pi * numpy.sum( terminal_tips_radii )

        self.total_length       = 0.
        self.slant_length       = 0.
        self.cylindric_volume       = 0.
        self.cylindric_lateral_area = 0.
        self.frustum_volume         = 0.
        self.frustum_lateral_area   = 0.
        morphology.root.volume      = 0
        morphology.root.lateral_area     = 0
        for compartment in morphology.non_root_compartments:
            self.total_length   += compartment.length
            compartment.frustum_length  = math.sqrt( (compartment.parent.radius - compartment.radius)**2 + compartment.length**2)
            self.slant_length   += compartment.frustum_length

            #mostly smaller radius  = compartment.radius
            radius  = (compartment.parent.radius + compartment.radius)/2.
            self.cylindric_volume   += radius**2 * compartment.length
            self.cylindric_lateral_area += radius * compartment.length

            self.frustum_volume         += compartment.length * (compartment.parent.radius**2 + compartment.parent.radius * compartment.radius + compartment.radius**2)            
            self.frustum_lateral_area   += compartment.frustum_length * (compartment.parent.radius + compartment.radius)
            
            compartment.volume          = compartment.length * (compartment.parent.radius**2 + compartment.parent.radius * compartment.radius + compartment.radius**2)          
            compartment.lateral_area    = compartment.frustum_length * (compartment.parent.radius + compartment.radius)
        self.cylindric_volume       *= math.pi
        self.cylindric_lateral_area *= 2 * math.pi
        self.cylindric_surface_area = self.cylindric_lateral_area + terminal_tips_cross_section_area
        self.frustum_volume         *= math.pi/3
        self.frustum_lateral_area   *= math.pi
        self.frustum_surface_area   = self.frustum_lateral_area   + terminal_tips_cross_section_area

        compartment_radii    = [compartment.radius for compartment in morphology.compartments]
        self.am_csa = 2 * math.pi * numpy.mean(compartment_radii)
        self.gm_csa  = 2 * math.pi * scipy.stats.stats.gmean(compartment_radii)
        try:
            self.hm_csa   = 2 * math.pi * scipy.stats.stats.hmean(compartment_radii)
        except Exception, e:
            self.hm_csa   = 0
        self.median_csa = 2 * math.pi * numpy.median(compartment_radii)

        branches_degrees    = {}
        branchingpoints_distances   = []
        for branch_point in morphology.branches.keys():
            branchingpoints_distances.append( branch_point.distance_path(branch_point.parent_node) )
            degree  = 0
            bp  = branch_point
            while vars(bp).has_key('parent_node') and bp.parent_node != None:
                bp = bp.parent_node
                degree += 1
            if not branches_degrees.has_key(degree):
                branches_degrees[degree]    = []
            branches_degrees[degree].append( morphology.branches[branch_point] )

        self.branches_degrees_mean_length    = {}
        self.branches_degrees_mean_volume    = {}
        self.branches_degrees_mean_lateral_area    = {}
        self.branches_degrees_mean_compactness    = {}
        for degree,branches in branches_degrees.iteritems():
            #===================================================================
            # lengths = []
            # volumes = []
            # lateral_areas = []
            # for branch in branches:
            #    lengths.append( numpy.sum([c.length for c in branch]) )
            #    volumes.append( numpy.sum([c.volume for c in branch]) )
            #    lateral_areas.append( numpy.sum([c.lateral_area for c in branch]) )
            #===================================================================
            self.branches_degrees_mean_length[degree]    = numpy.mean( [numpy.sum([c.length for c in branch]) for branch in branches]  )
            self.branches_degrees_mean_volume[degree]    = numpy.mean( [math.pi/3 * numpy.sum([compartment.length * (compartment.parent.radius**2 + compartment.parent.radius * compartment.radius + compartment.radius**2) for c in branch]) for branch in branches] )
            self.branches_degrees_mean_lateral_area[degree]    = numpy.mean( [math.pi * numpy.sum([compartment.frustum_length * (compartment.parent.radius + compartment.radius) for c in branch]) for branch in branches] )
            #self.branches_degrees_mean_compactness[degree]    = self.branches_degrees_mean_lateral_area[degree] / self.branches_degrees_mean_volume[degree]
            self.branches_degrees_mean_compactness[degree]     = numpy.sum([numpy.sum([c.lateral_area for c in branch]) for branch in branches]) / numpy.sum([numpy.sum([c.volume for c in branch]) for branch in branches]) 
        
        #=======================================================================
        # leafs   = [leaf for leaf in morphology.terminal_tips]
        # while len(leafs) > 0:
        #    new_leafs   = {}
        #    for compartment in leafs:
        #        if compartment.compartment_parent_id < 1:
        #            branchingpoints_distances.append(leaf.length/.2)
        #            continue
        #        distance = 0
        #        compartment   = compartment.parent
        #        while compartment.compartment_parent_id > 0 and len(compartment.children) == 1:
        #            distance += compartment.length
        #            compartment   = compartment.parent
        #        branchingpoints_distances.append(distance)
        #        if compartment.compartment_parent_id > 0:
        #            new_leafs[compartment]  = True
        #    leafs   = new_leafs.keys()
        #=======================================================================

        self.am_branch_point_distance= numpy.mean( branchingpoints_distances )
        self.gm_branch_point_distance = scipy.stats.stats.gmean( branchingpoints_distances )
        try:
            self.hm_branch_point_distance  = scipy.stats.stats.hmean( branchingpoints_distances )
        except Exception, e:
            self.hm_branch_point_distance  = 0
        self.median_branch_point_distance         = numpy.median( branchingpoints_distances )

        branchingpoints_radii   = [compartment.radius for compartment in morphology.branch_points]
        self.am_branch_point_csa = 2 * math.pi * numpy.mean(branchingpoints_radii)
        self.gm_branch_point_csa  = 2 * math.pi * scipy.stats.stats.gmean(branchingpoints_radii)
        try:
            self.hm_branch_point_csa   = 2 * math.pi * scipy.stats.stats.hmean(branchingpoints_radii)
        except Exception, e:
            self.hm_branch_point_csa   = 0
        self.median_branch_point_csa          = 2 * math.pi * numpy.median(branchingpoints_radii)
        
        terminal_tips_distances   = [morphology.root.distance_path(compartment) for compartment in morphology.terminal_tips]
        self.am_terminal_tip_distance = numpy.mean( terminal_tips_distances )
        self.gm_terminal_tip_distance   = scipy.stats.stats.gmean( terminal_tips_distances )
        try:
            self.hm_terminal_tip_distance    = scipy.stats.stats.hmean( terminal_tips_distances )
        except Exception, e:
            self.hm_terminal_tip_distance   = 0
        self.median_terminal_tip_distance = numpy.median( terminal_tips_distances )

        self.am_terminal_tip_csa = 2 * math.pi * numpy.mean(terminal_tips_radii)
        self.gm_terminal_tip_csa  = 2 * math.pi * scipy.stats.stats.gmean(terminal_tips_radii)
        try:
            self.hm_terminal_tip_csa   = 2 * math.pi * scipy.stats.stats.hmean(terminal_tips_radii)
        except Exception, e:
            self.hm_terminal_tip_csa   = 0
        self.median_terminal_tip_csa = 2 * math.pi * numpy.median(terminal_tips_radii)



    @property
    def cylindric_compactness(self):
        return self.cylindric_volume  / self.cylindric_lateral_area
    @property
    def frustum_compactness(self):
        return self.frustum_volume / self.cylindric_lateral_area
    @property
    def cylindric_mean_cross_section_area(self):
        return self.cylindric_volume / self.total_length
    @property
    def frustum_mean_cross_section_area(self):
        return self.frustum_volume / self.total_length
    @property
    def mcsa(self):
        return self.frustum_mean_cross_section_area
    @property
    def mean_branching_distance(self):
        return self.total_length / self.number_of_terminal_tips
    @property
    def mean_branch_point_distance(self):
        return self.total_length / self.number_of_branch_points

    @property
    def es_volume(self):
        '''
        volume of minimal sphere with equal surface area        [~m³]
        '''
        return 4./3 * math.pi * ( self.frustum_surface_area / 4. / math.pi) ** (3./2)
    @property
    def es_surface_area(self):
        '''
        surface area of minimal sphere with equal volume 
        '''
        return 4 * math.pi *  (3./4 * self.frustum_volume / math.pi) ** (2./3)
    @property
    def es_compactness(self):
        '''
        4/3 radius ?        [m]
        '''
        return self.es_volume / self.es_surface_area
    @property
    def esn_frustum_surface_area(self):
        '''
        surface area / equal volume minial shere surface area        [#]
        '''
        return self.frustum_surface_area / self.es_surface_area
    @property
    def esn_frustum_volume(self):
        '''
        volume / equal surface area sphere volume        [#]
        '''
        return self.frustum_volume / self.es_volume
    @property
    def esn_compactness(self):
        '''
        4/3 radius ?
        [m]
        TODO: == es ? 
        '''
        return self.esn_frustum_volume / self.esn_frustum_surface_area

    @property
    def volume(self):
        return self.frustum_volume
    @property
    def lateral_area(self):
        return self.frustum_lateral_area
    @property
    def surface_area(self):
        return self.frustum_surface_area
    @property
    def compactness(self):
        return self.volume / self.lateral_area

    def _pca(self):
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
        mins   = [float('inf'), float('inf'), float('inf')]
        maxs   = [float('-inf'),float('-inf'),float('-inf')]
        for x in mdp.pca( numpy.array([[compartment.x, compartment.y, compartment.z] for compartment in self.morphology.compartments]) ):
            for d in xrange(3):
                if x[d] < mins[d]:
                    mins[d]     = x[d]
                if x[d] > maxs[d]:
                    maxs[d]     = x[d]
        self._pca_length_x  = maxs[0] - mins[0]
        self._pca_length_y  = maxs[1] - mins[1]
        self._pca_length_z  = maxs[2] - mins[2]
        
        self._pca_lengths = (
            self._pca_length_x,
            self._pca_length_y,
            self._pca_length_z
        )
    @property
    def pca_length_x(self):
        if not vars(self).has_key('_pca_length_x') or self._pca_length_x == None:
            self._pca()
        return self._pca_length_x
    @property
    def pca_length_y(self):
        if not vars(self).has_key('_pca_length_y') or self._pca_length_y == None:
            self._pca()
        return self._pca_length_y
    @property
    def pca_length_z(self):
        if not vars(self).has_key('_pca_length_z') or self._pca_length_z == None:
            self._pca()
        return self._pca_length_z
    @property
    def pca_lengths(self):
        if not vars(self).has_key('_pca_lengths') or self._pca_lengths == None:
            self._pca()
        return self._pca_lengths

    def _cep(self):
        '''
        convex_enveloping_polyhedron_volume / convex_enveloping_polyhedron_lateral_area        [m]
        frustum volume / polygon volume        [#]
        frustum surface area / polygon surface area        [#]
        frustum / polygon (volume / surface area)        [#]
        '''
        try:
            self.cep_hull    = morphjongleur.util.chull.Hull([morphjongleur.util.chull.Vector.fromArray([compartment.x, compartment.y, compartment.z]) for compartment in self.morphology.compartments])
            convex_enveloping_polyhedron_surface_area, convex_enveloping_polyhedron_volume  = self.cep_hull.surface_area_and_volume()
        except Exception, e:
            convex_enveloping_polyhedron_surface_area, convex_enveloping_polyhedron_volume  = float('nan'),float('nan')
            print traceback.format_exc()
        self._cep_surface_area  = convex_enveloping_polyhedron_surface_area
        self._cep_volume        = convex_enveloping_polyhedron_volume
    @property
    def cep_surface_area(self):
        if not vars(self).has_key('_cep_surface_area') or self._cep_surface_area == None:
            self._cep()
        return self._cep_surface_area
    @property
    def cep_volume(self):
        if not vars(self).has_key('_cep_volume') or self._cep_volume == None:
            self._cep()
        return self._cep_volume
    @property
    def cep_compactness(self):
        return self.cep_volume / self.cep_surface_area
    @property
    def cepn_volume(self):
        return self.volume / self.cep_volume
    @property
    def cepn_surface_area(self):
        return self.surface_area / self.cep_surface_area
    @property
    def cepn_compactness(self):
        '''
        cepn_volume / cepn_surface_area = compactness / cep_compactness
        '''
        return self.cepn_volume / self.cepn_surface_area
    @property
    def esn_cep_volume(self):
        return self.cep_volume / self.es_volume
    @property
    def esn_cep_surface_area(self):
        return self.cep_surface_area / self.es_surface_area
    @property
    def esn_cep_compactness(self):
        return self.esn_cep_volume / self.esn_cep_surface_area

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
        #matplotlib.pyplot.title('change: DB: %f, VB: %f' %(float(morphologies[1].info.__dict__[quantity])/morphologies[0].info.__dict__[quantity], float(morphologies[3].info.__dict__[quantity])/morphologies[2].info.__dict__[quantity]))#TODO: verallgemeinern, verifizierren!
        print '%s change: DB: %f, VB: %f' %(quantity, float(morphologies[1].info.__dict__[quantity])/morphologies[0].info.__dict__[quantity], float(morphologies[3].info.__dict__[quantity])/morphologies[2].info.__dict__[quantity])
        matplotlib.pyplot.xticks(ind+width, xticks,rotation=9)#
        #matplotlib.pyplot.yticks(numpy.arange(0,81,10))
        matplotlib.pyplot.grid(True, color='lightgrey')
        
        matplotlib.pyplot.bar(ind, values, width, color='black')#

        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+quantity+'.'+picture_format, format=picture_format, transparent=True)
                except Exception, e:
                    print picture_format 
                    print traceback.format_exc()
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

    @staticmethod
    def bars_plot(v, bars, xs, colors, y_label='', rotation=0, horizontal=False, tex=False, ratio=None, xlimit=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        matplotlib.rc('text', usetex=tex)
        v_means = {}
        v_stds  = {}
        for b in bars:
            v_means[b]  = [ numpy.mean( v[b][x] ) for x in xs ]
            v_stds[b]   = [ numpy.std( v[b][x] )  for x in xs ]

        ind = numpy.arange( len(xs) )  # the x locations for the groups*width
        width = 1.0/(len(bars)+1)      # the width of the bars

        matplotlib.pyplot.subplot(111)
        rects   = []
        for i in xrange(len(bars)):
            if horizontal:
                rect    = matplotlib.pyplot.barh(ind+i*width+width/2., 
                                v_means[ bars[i] ],#TODO: (numpy.mean( v[b][x] ) for x in xs) 
                                width,
                                color=colors[i],
                                #yerr=v_stds[ bars[i] ],#(numpy.std( v[b][x] ) for x in xs)
                                #error_kw=dict(elinewidth=6, ecolor='red')
                                )
            else:
                rect    = matplotlib.pyplot.bar(ind+i*width+width/2., 
                            v_means[ bars[i] ],#TODO: (numpy.mean( v[b][x] ) for x in xs) 
                            width,
                            color=colors[i],
                            #yerr=v_stds[ bars[i] ],#(numpy.std( v[b][x] ) for x in xs)
                            #error_kw=dict(elinewidth=6, ecolor='red')
                            )
            rects.append(rect[0])

        xs  = [str(x).replace('_',' ') for x in xs]
        if horizontal:
            matplotlib.pyplot.xlabel(y_label)
            #matplotlib.pyplot.ylim(ymin=-1)
            matplotlib.pyplot.yticks(ind + (len(bars)-1)*width+width/2., xs , rotation=rotation )
            matplotlib.pyplot.axvline(x=0, color='black')
        else:
            matplotlib.pyplot.ylabel(y_label)
            #matplotlib.pyplot.ylim(ymin=-1)
            matplotlib.pyplot.xticks(ind + (len(bars)-1)*width+width/2., xs, rotation=rotation )
            matplotlib.pyplot.axhline(y=0, color='black')

        matplotlib.pyplot.grid(True, axis='y', color='lightgrey')
        
        if ratio != None:#matplotlib.figure.figaspect(arg)
            fig = matplotlib.pyplot.gcf()
            fig.set_size_inches(ratio[0],ratio[1])
        
        if xlimit != None:#matplotlib.figure.figaspect(arg)
            matplotlib.pyplot.xlim(xlimit)

        matplotlib.pyplot.legend( tuple(rects), tuple(map(str, bars)), loc='best')
        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig('%s.%s' % (str(picture_file), str(picture_format)), format=picture_format, transparent=True)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()
        
    @staticmethod
    def farther_away(compartment_iterable, center, upper_threshold=450, lower_threshold=0):
        for c in compartment_iterable:
            distance    = center.distance_path(c)
            if lower_threshold > distance or distance > upper_threshold:
                yield(c)


if __name__ == '__main__':
    '''
    Parameter: files, not directories

    cd data/swc_files/
    path_to/metric_analysis.py *.swc
    '''
    import sys
    import morphjongleur.util.parser.swc
    import morphjongleur.util.transformations
    picture_formats = ['png','svg', 'pdf']#
    colors = ['#000000', '#00ff00','#008000', '#0000ff','#000080']# H060602DB_10_2(whole).swc H060607DB_10_2(whole).swc H060602VB_10_2(whole).swc H060607VB_10_2(whole).swc

    with_head   = True
    i = 0
    for swc in sys.argv[1:]:#['../../data/test.swc']:#
        color = colors[i % len(colors)]
        i += 1

        morphology   = Morphology.swc_parse(swc, verbose=False)
        
        m_pca   = morphology.pca()
        Compartment.write_svg('/tmp/%s_pca_faraway.svg' % (m_pca.name), 
            [m_pca.compartments,  MetricAnalysis.farther_away(m_pca.compartments, m_pca.root, 650)], 
            [colors[1], 'black']
        )

        morphology.root.plot_distance(morphology.compartments,       morphology.name,            xlim=900, ylim=8,   color=color, ratio=(4,4.5), picture_file='/tmp/distance_compartments_'+str(morphology.name),               picture_formats=picture_formats)
        morphology.root.plot_distance(morphology.branch_points,   morphology.name,            xlim=900, ylim=8,   color=color, ratio=(4,4.5), picture_file='/tmp/distance_branch_points_'+str(morphology.name),               picture_formats=picture_formats)
        morphology.root.plot_distance(morphology.terminal_tips,       morphology.name,            xlim=900, ylim=8,   color=color, ratio=(4,4.5), picture_file='/tmp/distance_terminal_tips_'+str(morphology.name),               picture_formats=picture_formats)

        morphology.plot_distance_distribution(morphology.root,                                    morphology.name,   color=color,    bins=20,   xlim=900, ylim=900, ratio=(4,4.5), picture_file='/tmp/distance_distribution_compartments_'+str(morphology.name),  picture_formats=picture_formats)
        morphology.plot_distance_distributions([morphology.branch_points], [morphology.root],  [morphology.name], colors=[color], bins=20,   xlim=900, ylim=65, ratio=(4,4.5),  picture_file='/tmp/distance_distribution_branch_points_'+str(morphology.name),  picture_formats=picture_formats)
        morphology.plot_distance_distributions([morphology.terminal_tips],     [morphology.root],  [morphology.name], colors=[color], bins=20,   xlim=900, ylim=70, ratio=(4,4.5),  picture_file='/tmp/distance_distribution_terminal_tips_'+str(morphology.name),  picture_formats=picture_formats)#900, 43


        print morphology.name
        continue

        a   = MetricAnalysis(morphology)
    
        (ks, vs)    = a.variable_table(['name', 'total_length', 'slant_length', 'volume', 'lateral_area', 'mcsa', 'compactness', 'surface_area', 
            'pca_length_x', 'pca_length_y', 'pca_length_z', 
            'number_of_terminal_tips', 'number_of_branch_points', 
            'mean_branching_distance', 
            'es_volume', 'es_surface_area', 'es_compactness', 
            'esn_frustum_volume', 'esn_frustum_surface_area', 'esn_compactness',
            'cep_volume', 'cep_surface_area', 'cep_compactness', 
            'cepn_volume','cepn_surface_area', 'cepn_compactness',
            'esn_cep_volume', 'esn_cep_surface_area', 'esn_cep_compactness',
            'am_branch_point_distance', 'gm_branch_point_distance', 'hm_branch_point_distance', 'median_branch_point_distance',
            'am_branch_point_csa', 'gm_branch_point_csa', 'hm_branch_point_csa', 'median_branch_point_csa',
            'am_terminal_tip_distance', 'gm_terminal_tip_distance', 'hm_terminal_tip_distance', 'median_terminal_tip_distance',
            'am_terminal_tip_csa', 'gm_terminal_tip_csa', 'hm_terminal_tip_csa', 'median_terminal_tip_csa',
            ])#,g=6,texify=True)

        if vars(a).has_key('cep_hull'):
            a.cep_hull.write('/tmp/%s_pca' % (morphology.name))
        
        if with_head:
            print ks
            with_head   = False
        print vs

        #print morphology.name
        morphology.write_svg(svg_file='/tmp/%s.svg' % (morphology.name), color=color)
        Compartment.write_svg('/tmp/%s_color.svg' % (morphology.name), 
            [morphology.compartments,  morphology.terminal_tips, [morphology.root]], 
            [color, 'yellow', 'red']
        )

        #=======================================================================
        # morphology.plot(color=color, picture_file='/tmp/%s' % (morphology.name), picture_formats=picture_formats)
        # Compartment.plot(
        #    [morphology.compartments,  morphology.terminal_tips, [morphology.root]], 
        #    [color, 'yellow', 'red'], 
        #    picture_file='/tmp/%s_color' % (morphology.name), picture_formats=picture_formats#
        # )
        #=======================================================================
        Compartment.plot3d(morphology.compartments, ratio=(8,6), picture_file='/tmp/%s_3d' % (morphology.name), picture_formats=picture_formats)

        continue
        m_pca   = morphology.pca()
        m_pca.swc_write('/tmp/%s_pca.swc' % (m_pca.name) )
        m_pca.plot(color=color, picture_file='/tmp/%s_pca' % (m_pca.name), picture_formats=picture_formats)
        m_pca.write_svg(svg_file='/tmp/%s_pca.svg' % (m_pca.name), color=color)
        Compartment.plot(
            [m_pca.compartments,  m_pca.terminal_tips, [m_pca.root]], 
            [color, 'yellow', 'red'],
            picture_file='/tmp/%s_pca_color' % (m_pca.name), picture_formats=picture_formats#
        )
        Compartment.write_svg('/tmp/%s_pca_color.svg' % (morphology.name), 
            [morphology.compartments,  morphology.terminal_tips, [morphology.root]], 
            [color, 'yellow', 'red']
        )

        #continue
        #continue
        #------------------------------------------------------------------ try:
            #-------------------------------------------------------------- pass
            #-------------------------------------- #a   = MetricAnalysis(m_pca)
            #------------------- #a.cep_hull.write('/tmp/%s_pca' % (m_pca.name))
        #-------------------------------------------------- except Exception, e:
            #-------------------------------------- print traceback.format_exc()




    #sys.exit(0)

    morphologies    = [Morphology.swc_parse(swc, verbose=False) for swc in sys.argv[1:6]]

    a   = [MetricAnalysis(ma).variable_map()[0] for ma in morphologies[1:5]]
    v   = {'dorsal branch':{}, 'ventral branch':{}}
    bars= ['dorsal branch','ventral branch']
    xs  = ['total_length', 'slant_length', 'volume', 'lateral_area', 'mcsa', 'compactness', 'surface_area', 
    'pca_length_x', 'pca_length_y', 'pca_length_z', 
    'number_of_terminal_tips', 'number_of_branch_points', 
    'mean_branching_distance', 
    'es_volume', 'es_surface_area', 'es_compactness', 
    'esn_frustum_volume', 'esn_frustum_surface_area', 'esn_compactness',
    'cep_volume', 'cep_surface_area', 'cep_compactness', 
    'cepn_volume','cepn_surface_area', 'cepn_compactness',
    'esn_cep_volume', 'esn_cep_surface_area', 'esn_cep_compactness',
    ]
    for key in xs:#TODO: more time efficient with properties and only needed
        if a[0][key] != 0 and a[2][key] != 0:
            v['dorsal branch'][key]     = float(a[1][key]) / a[0][key] - 1
            v['ventral branch'][key]    = float(a[3][key]) / a[2][key] - 1
        else: 
            v['dorsal branch'][key]     = float('nan')
            v['ventral branch'][key]    = float('nan')
    xs.reverse()
    MetricAnalysis.bars_plot(
        v=v, 
        bars=bars, 
        xs=xs, 
        colors=['#00ff00','#0000ff'], 
        rotation=0, 
        horizontal=True,
        ratio=(8,11), #8.268,11.6
        y_label='change: forager / nurse - 1', 
        picture_file='/tmp/change_morphology', picture_formats=picture_formats
    )


    morphologies[0].name    = 'test'
    morphologies[1].name    = 'nurse'
    morphologies[2].name    = 'forager'
    morphologies[3].name    = 'nurse'
    morphologies[4].name    = 'forager'

    m_pca   = morphology.pca()
    Compartment.write_svg('/tmp/%s_pca_faraway.svg' % (m_pca.name), 
        [m_pca.compartments,  MetricAnalysis.farther_away(m_pca.compartments, m_pca.root, 650)], 
        [colors[1], 'black']
    )
    sys.exit(0)

    Compartment.plot_distances([m.compartments for m in morphologies[1:3]],
                               [m.root for m in morphologies[1:3]], 
                               [m.name for m in morphologies[1:3]],            
                               xlim=900, ylim=8,   
                               colors=colors[1:3], 
                               picture_file='/tmp/distance_compartments_db',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.branch_points for m in morphologies[1:3]],
                               [m.root for m in morphologies[1:3]], 
                               [m.name for m in morphologies[1:3]],            
                               xlim=900, ylim=8,   
                               colors=colors[1:3], 
                               picture_file='/tmp/distance_branch_points_db',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.terminal_tips for m in morphologies[1:3]],
                               [m.root for m in morphologies[1:3]], 
                               [m.name for m in morphologies[1:3]],            
                               xlim=900, ylim=8,   
                               colors=colors[1:3], 
                               picture_file='/tmp/distance_terminal_tips_db',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.compartments for m in morphologies[3:5]],
                               [m.root for m in morphologies[3:5]], 
                               [m.name for m in morphologies[3:5]],            
                               xlim=900, ylim=8,   
                               colors=colors[3:5], 
                               picture_file='/tmp/distance_compartments_vb',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.branch_points for m in morphologies[3:5]],
                               [m.root for m in morphologies[3:5]], 
                               [m.name for m in morphologies[3:5]],            
                               xlim=900, ylim=8,   
                               colors=colors[3:5], 
                               picture_file='/tmp/distance_branch_points_vb',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.terminal_tips for m in morphologies[3:5]],
                               [m.root for m in morphologies[3:5]], 
                               [m.name for m in morphologies[3:5]],            
                               xlim=900, ylim=8,   
                               colors=colors[3:5], 
                               picture_file='/tmp/distance_terminal_tips_vb',   
                               picture_formats=['png','svg']
    )

    Morphology.plot_distance_distributions(
        [m.terminal_tips for m in morphologies[1:3]],
        [m.root for m in morphologies[1:3]], 
        [m.name for m in morphologies[1:3]], 
        colors=colors[1:3], 
        bins=20, xlim=900, ylim=70, 
        picture_file='/tmp/distance_distribution_terminal_tips_db',               
        picture_formats=['png','svg']
    )
    Morphology.plot_distance_distributions(
        [m.terminal_tips for m in morphologies[3:5]],
        [m.root for m in morphologies[3:5]], 
        [m.name for m in morphologies[3:5]], 
        colors=colors[3:5], 
        bins=20, xlim=900, ylim=70, #xlim=600, ylim=40, 
        picture_file='/tmp/distance_distribution_terminal_tips_vb',              
        picture_formats=['png','svg']
    )

    Morphology.plot_distance_distributions(
        [m.branch_points for m in morphologies[1:3]],
        [m.root for m in morphologies[1:3]], 
        [m.name for m in morphologies[1:3]], 
        colors=colors[1:3], 
        bins=20, xlim=900, ylim=65, 
        picture_file='/tmp/distance_distribution_branch_points_db',               
        picture_formats=['png','svg']
    )
    Morphology.plot_distance_distributions(
        [m.branch_points for m in morphologies[3:5]],
        [m.root for m in morphologies[3:5]], 
        [m.name for m in morphologies[3:5]], 
        colors=colors[3:5], 
        bins=20, xlim=900, ylim=65, #xlim=600, ylim=40, 
        picture_file='/tmp/distance_distribution_branch_points_vb',               
        picture_formats=['png','svg']
    )

    Morphology.plot_distance_distributions(
        [m.compartments for m in morphologies[1:3]],
        [m.root for m in morphologies[1:3]], 
        [m.name for m in morphologies[1:3]], 
        colors=colors[1:3], 
        bins=20, xlim=900, ylim=900, 
        picture_file='/tmp/distance_distribution_compartments_db',               
        picture_formats=['png','svg']
    )
    Morphology.plot_distance_distributions(
        [m.compartments for m in morphologies[3:5]],
        [m.root for m in morphologies[3:5]], 
        [m.name for m in morphologies[3:5]], 
        colors=colors[3:5], 
        bins=20, xlim=900, ylim=900, #xlim=600, ylim=700, 
        picture_file='/tmp/distance_distribution_compartments_vb',               
        picture_formats=['png','svg']
    )
