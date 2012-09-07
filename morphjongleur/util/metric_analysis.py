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
    surface_length_frustum    total cell length running over the surface of an frustum
    
    
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
        self.compartments       = morphology.number_of_compartments
        
        self.number_of_terminaltips     = morphology.number_of_terminaltips
        self.number_of_branching_points = morphology.number_of_branching_points
#        if self.number_of_terminaltips != self.number_of_branching_points + 1:
#            print >> sys.stderr, "#terminaltips %i > #branching points %i +1, because of %s" % (self.number_of_terminaltips, self.number_of_branching_points, ["%i:%i"%(pleb.compartment_id,len(pleb.children)) for pleb in morphology.plebs])

        terminaltips_radii      = [compartment.radius for compartment in morphology.terminaltips]
        terminaltips_cross_section_area = 2 * math.pi * numpy.sum( terminaltips_radii )

        self.total_cell_length      = 0.
        self.surface_length_frustum= 0.
        self.cylindric_volume       = 0.
        self.cylindric_lateral_area = 0.
        self.frustum_volume         = 0.
        self.frustum_lateral_area   = 0.
        for compartment in morphology.non_root_compartments():
            self.total_cell_length  += compartment.length
            compartment.frustum_length  = math.sqrt( (compartment.parent.radius - compartment.radius)**2 + compartment.length**2)
            self.surface_length_frustum+= compartment.frustum_length

            self.cylindric_volume   += compartment.radius**2 * compartment.length
            self.cylindric_lateral_area += compartment.radius * compartment.length

            self.frustum_volume         += compartment.length * (compartment.parent.radius**2 + compartment.parent.radius * compartment.radius + compartment.radius**2)            
            self.frustum_lateral_area   += compartment.frustum_length * (compartment.parent.radius + compartment.radius)
        self.cylindric_volume       *= math.pi
        self.cylindric_lateral_area *= 2 * math.pi
        self.cylindric_surface_area = self.cylindric_lateral_area + terminaltips_cross_section_area
        self.cylindric_sparsity     = self.cylindric_surface_area / self.cylindric_volume 
        self.frustum_volume         *= math.pi/3
        self.frustum_lateral_area   *= math.pi
        self.frustum_surface_area   = self.frustum_lateral_area   + terminaltips_cross_section_area
        self.frustum_sparsity       = self.cylindric_surface_area / self.frustum_volume

        self.cylindric_mean_cross_section_area  = self.cylindric_volume / self.total_cell_length
        self.frustum_mean_cross_section_area    = self.frustum_volume / self.total_cell_length
        compartment_radii    = [compartment.radius for compartment in morphology.compartments]
        self.arithmetic_mean_cross_section_area = 2 * math.pi * numpy.mean(compartment_radii)
        self.geometric_mean_cross_section_area  = 2 * math.pi * scipy.stats.stats.gmean(compartment_radii)
        self.harmonic_mean_cross_section_area   = 2 * math.pi * scipy.stats.stats.hmean(compartment_radii)
        self.median_cross_section_area = 2 * math.pi * numpy.median(compartment_radii)

        self.mean_branching_distance    = self.total_cell_length / self.number_of_terminaltips
        self.mean_branchpoint_distance  = self.total_cell_length / self.number_of_branching_points
        branchingpoints_distances   = []
        leafs   = [leaf for leaf in morphology.terminaltips]
        while len(leafs) > 0:
            new_leafs   = {}
            for leaf in leafs:
                compartment = leaf.parent
                distance = leaf.length
                while compartment.compartment_parent_id > 0 and len(compartment.children) == 1:
                    distance += compartment.length
                    compartment   = compartment.parent
                branchingpoints_distances.append(distance)
                if compartment.compartment_parent_id > 0:
                    new_leafs[compartment]  = True
            leafs   = new_leafs.keys()
        self.arithmetic_mean_branchpoint_distance= numpy.mean( branchingpoints_distances )
        self.geometric_mean_branchpoint_distance = scipy.stats.stats.gmean( branchingpoints_distances )
        self.harmonic_mean_branchpoint_distance  = scipy.stats.stats.hmean( branchingpoints_distances )
        self.median_branchpoint_distance         = numpy.median( branchingpoints_distances )

        branchingpoints_radii   = [compartment.radius for compartment in morphology.branching_points]
        self.arithmetic_mean_branchpoint_cross_section_area = 2 * math.pi * numpy.mean(branchingpoints_radii)
        self.geometric_mean_branchpoint_cross_section_area  = 2 * math.pi * scipy.stats.stats.gmean(branchingpoints_radii)
        self.harmonic_mean_branchpoint_cross_section_area   = 2 * math.pi * scipy.stats.stats.hmean(branchingpoints_radii)
        self.median_branchpoint_cross_section_area          = 2 * math.pi * numpy.median(branchingpoints_radii)
        
        terminaltips_distances   = [morphology.root_biggest_child.path_distance(compartment) for compartment in morphology.terminaltips]
        self.arithmetic_mean_terminaltip_distance = numpy.mean( terminaltips_distances )
        self.geometric_mean_terminaltip_distance   = scipy.stats.stats.gmean( terminaltips_distances )
        self.harmonic_mean_terminaltip_distance    = scipy.stats.stats.hmean( terminaltips_distances )
        self.median_terminaltip_distance = numpy.median( terminaltips_distances )

        self.arithmetic_mean_terminaltip_cross_section_area = 2 * math.pi * numpy.mean(terminaltips_radii)
        self.geometric_mean_terminaltip_cross_section_area  = 2 * math.pi * scipy.stats.stats.gmean(terminaltips_radii)
        self.harmonic_mean_terminaltip_cross_section_area   = 2 * math.pi * scipy.stats.stats.hmean(terminaltips_radii)
        self.median_terminaltip_cross_section_area = 2 * math.pi * numpy.median(terminaltips_radii)


        mins   = [float('inf'), float('inf'), float('inf')]
        maxs   = [float('-inf'),float('-inf'),float('-inf')]
        for x in mdp.pca( numpy.array([[compartment.x, compartment.y, compartment.z] for compartment in morphology.compartments]) ):
            for d in range(3):
                if x[d] < mins[d]:
                    mins[d]     = x[d]
                if x[d] > maxs[d]:
                    maxs[d]     = x[d]
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

        self.cep_hull    = morphjongleur.util.chull.Hull([morphjongleur.util.chull.Vector.fromArray([compartment.x, compartment.y, compartment.z]) for compartment in morphology.compartments])
        convex_enveloping_polyhedron_surface_area, convex_enveloping_polyhedron_volume   = self.cep_hull.surface_area_and_volume()
        self.cep_surface_area  = convex_enveloping_polyhedron_surface_area
        self.cep_volume        = convex_enveloping_polyhedron_volume

        '''
        volume of minial shere with equal surface area        [~m³]
        surface area of minial shere with equal volume 
        4/3 radius ?        [m]
        '''
        self.es_volume          = 4./3 * math.pi * ( self.frustum_surface_area / 4. / math.pi) ** (3./2)
        self.es_surface_area    = 4 * math.pi *  (3./4 * self.frustum_volume / math.pi) ** (2./3)
        self.es_sparsity        = self.es_surface_area / self.es_volume 

        '''
        surface area / equal volume minial shere surface area        [#]
        volume / equal surface area sphere volume        [#]
        4/3 radius ?
        [m]
        TODO: == es ? 
        '''
        self.esn_frustum_surface_area   = self.frustum_surface_area / self.es_surface_area
        self.esn_frustum_volume         = self.frustum_volume / self.es_volume
        self.esn_sparsity               = self.esn_frustum_surface_area / self.esn_frustum_volume
    
        '''
        convex_enveloping_polyhedron_volume / convex_enveloping_polyhedron_surface_area        [m]
        frustum volume / polygon volume        [#]
        frustum surface area / polygon surface area        [#]
        frustum / polygon (volume / surface area)        [#]
        '''
        self.cep_volume_div_cep_surface_area    = self.cep_volume / self.cep_surface_area
        self.cepn_volume    = self.frustum_volume / self.cep_volume
        self.cepn_surface_area     = self.frustum_surface_area / self.cep_surface_area
        self.cepn_sparsity   = self.cepn_surface_area / self.cepn_volume
    
        self.esn_cep_volume         = self.cep_volume / self.es_volume
        self.esn_cep_surface_area   = self.cep_surface_area / self.es_surface_area
        self.esn_cep_sparsity       = self.esn_cep_volume / self.esn_cep_surface_area
    
    @property
    def pesn_cep_sparsity(self):
        if not vars(self).has_key('_esn_cep_sparsity') or self._esn_cep_sparsity == None:
            self._esn_cep_sparsity   = self.esn_cep_volume / self.esn_cep_surface_area
        return self._esn_cep_sparsity

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
        matplotlib.pyplot.title('change: DB: %f, VB: %f' %(float(morphologies[1].info.__dict__[quantity])/morphologies[0].info.__dict__[quantity], float(morphologies[3].info.__dict__[quantity])/morphologies[2].info.__dict__[quantity]))#TODO: verallgemeinern, verifizierren!
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
    def bars_plot(v, bars, xs, colors, y_label='', rotation=0, horizontal=False, tex=False, ratio=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
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
        for i in range(len(bars)):
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

        #xs  = ['' for x in xs]
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
        
        if ratio != None:
            fig = matplotlib.pyplot.gcf()
            fig.set_size_inches(ratio[0],ratio[1])

        matplotlib.pyplot.legend( tuple(rects), tuple(map(str, bars)), loc='best')
        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=True)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()
        
    @staticmethod
    def farther_away(compartment_iterable, center, upper_threshold=450, lower_threshold=0):
        for c in compartment_iterable:
            distance    = center.path_distance(c)
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
    colors = ['#000000', '#00ff00','#008000', '#0000ff','#000080']# H060602DB_10_2(whole).swc H060607DB_10_2(whole).swc H060602VB_10_2(whole).swc H060607VB_10_2(whole).swc
    with_head   = True
    i = 0
    for swc in sys.argv[1:]:#['../../data/test.swc']:#
        color = colors[i % len(colors)]
        i += 1

        morphology   = Morphology.swc_parse(swc, verbose=False)
        print morphology.name
        #morphology.plot(color=color, picture_file='/tmp/%s' % (morphology.name), picture_formats=['png', 'svg'])
        Compartment.plot(
            [morphology.compartments,  morphology.terminaltips, [morphology.root_biggest_child]], 
            [color, 'yellow', 'red'],
            picture_file='/tmp/%s_color' % (morphology.name), picture_formats=['png', 'svg']#
        )

        continue

        Compartment.plot(
            [morphology.compartments,  MetricAnalysis.farther_away(morphology.compartments, morphology.root_biggest_child, 450)], 
            [color, 'black'],
            picture_file='/tmp/%s_faraway' % (morphology.name), picture_formats=['png', 'svg']#
        )
        m_pca   = morphology.pca()
        #m_pca.swc_write('/tmp/%s_pca.swc' % (m_pca.name) )
        #m_pca.plot(color=color, picture_file='/tmp/%s_pca' % (m_pca.name), picture_formats=['png', 'svg'])
        Compartment.plot(
            [m_pca.compartments, MetricAnalysis.farther_away(m_pca.compartments, m_pca.root_biggest_child, 450)], 
            [color, 'black'],
            picture_file='/tmp/%s_color_pca' % (m_pca.name), picture_formats=['png', 'svg']#
        )
        Compartment.plot(
            [m_pca.compartments,  MetricAnalysis.farther_away(m_pca.compartments, m_pca.root_biggest_child, 450)], 
            [color, 'black'],
            picture_file='/tmp/%s_faraway_pca' % (m_pca.name), picture_formats=['png', 'svg']#
        )
        try:
            a   = MetricAnalysis(m_pca)
            a.cep_hull.write('/tmp/%s_pca' % (m_pca.name))
        except Exception, e:
            print traceback.format_exc()

        continue

        morphology.root_biggest_child.plot_distance(morphology.compartments,       morphology.name,            xlim=900, ylim=8,   color=color, picture_file='/tmp/distance_compartments_'+str(morphology.name),               picture_formats=['png'])
        morphology.root_biggest_child.plot_distance(morphology.branching_points,   morphology.name,            xlim=900, ylim=8,   color=color, picture_file='/tmp/distance_branchpoints_'+str(morphology.name),               picture_formats=['png'])
        morphology.root_biggest_child.plot_distance(morphology.terminaltips,       morphology.name,            xlim=900, ylim=8,   color=color, picture_file='/tmp/distance_terminaltips_'+str(morphology.name),               picture_formats=['png'])

        morphology.plot_distance_distribution(morphology.root_biggest_child,                                    morphology.name,   color=color,    bins=20,   xlim=900, ylim=900, picture_file='/tmp/distance_distribution_compartments_'+str(morphology.name),  picture_formats=['png'])
        morphology.plot_distance_distributions([morphology.branching_points], [morphology.root_biggest_child],  [morphology.name], colors=[color], bins=20,   xlim=900, ylim=65,  picture_file='/tmp/distance_distribution_branchpoints_'+str(morphology.name),  picture_formats=['png'])
        morphology.plot_distance_distributions([morphology.terminaltips],     [morphology.root_biggest_child],  [morphology.name], colors=[color], bins=20,   xlim=900, ylim=70,  picture_file='/tmp/distance_distribution_terminaltips_'+str(morphology.name),  picture_formats=['png'])#900, 43

        a   = MetricAnalysis(morphology)
        a.cep_hull.write('/tmp/%s_pca' % (morphology.name))

        (ks, vs)    = a.variable_table(['name', 'compartments', #'datetime_recording', 
        'total_cell_length', 'surface_length_frustum', 
        'number_of_branching_points', 'number_of_terminaltips', 

        'cylindric_volume', 'cylindric_surface_area', 'cylindric_sparsity', 
        'frustum_volume', 'frustum_surface_area', 'frustum_sparsity', 
        
        'cylindric_mean_cross_section_area', 'frustum_mean_cross_section_area', 
        'arithmetic_mean_cross_section_area', 'geometric_mean_cross_section_area', 'harmonic_mean_cross_section_area', 'median_cross_section_area',
        'mean_branchpoint_distance', 'mean_branching_distance', 
        'arithmetic_mean_branchpoint_distance', 'geometric_mean_branchpoint_distance', 'harmonic_mean_branchpoint_distance', 'median_branchpoint_distance',
        'arithmetic_mean_branchpoint_cross_section_area', 'geometric_mean_branchpoint_cross_section_area', 'harmonic_mean_branchpoint_cross_section_area', 'median_branchpoint_cross_section_area',
        'arithmetic_mean_terminaltip_distance', 'geometric_mean_terminaltip_distance', 'harmonic_mean_terminaltip_distance', 'median_terminaltip_distance',
        'arithmetic_mean_terminaltip_cross_section_area', 'geometric_mean_terminaltip_cross_section_area', 'harmonic_mean_terminaltip_cross_section_area', 'median_terminaltip_cross_section_area',

        'pca_length_x', 'pca_length_y', 'pca_length_z', 
        'es_volume', 'es_surface_area', 'es_sparsity', 
        'esn_frustum_surface_area', 'esn_frustum_volume', 'es_sparsity',
        'cep_volume', 'cep_surface_area', 'cep_volume_div_cep_surface_area', 
        'cepn_volume','cepn_surface_area', 'cepn_sparsity',
        'esn_cep_volume', 'esn_cep_surface_area', 'esn_cep_sparsity'
        ])
        
        if with_head:
            print ks
            with_head   = False
        print vs

    sys.exit(0)

    morphologies    = [Morphology.swc_parse(swc, verbose=False) for swc in sys.argv[1:6]]
    morphologies[0].name    = 'test'
    morphologies[1].name    = 'nurse'
    morphologies[2].name    = 'forager'
    morphologies[3].name    = 'nurse'
    morphologies[4].name    = 'forager'

    Compartment.plot_distances([m.compartments for m in morphologies[1:3]],
                               [m.root_biggest_child for m in morphologies[1:3]], 
                               [m.name for m in morphologies[1:3]],            
                               xlim=900, ylim=8,   
                               colors=colors[1:3], 
                               picture_file='/tmp/distance_compartments_db',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.branching_points for m in morphologies[1:3]],
                               [m.root_biggest_child for m in morphologies[1:3]], 
                               [m.name for m in morphologies[1:3]],            
                               xlim=900, ylim=8,   
                               colors=colors[1:3], 
                               picture_file='/tmp/distance_branchpoints_db',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.terminaltips for m in morphologies[1:3]],
                               [m.root_biggest_child for m in morphologies[1:3]], 
                               [m.name for m in morphologies[1:3]],            
                               xlim=900, ylim=8,   
                               colors=colors[1:3], 
                               picture_file='/tmp/distance_terminaltips_db',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.compartments for m in morphologies[3:5]],
                               [m.root_biggest_child for m in morphologies[3:5]], 
                               [m.name for m in morphologies[3:5]],            
                               xlim=900, ylim=8,   
                               colors=colors[3:5], 
                               picture_file='/tmp/distance_compartments_vb',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.branching_points for m in morphologies[3:5]],
                               [m.root_biggest_child for m in morphologies[3:5]], 
                               [m.name for m in morphologies[3:5]],            
                               xlim=900, ylim=8,   
                               colors=colors[3:5], 
                               picture_file='/tmp/distance_branchpoints_vb',   
                               picture_formats=['png','svg']
    )
    Compartment.plot_distances([m.terminaltips for m in morphologies[3:5]],
                               [m.root_biggest_child for m in morphologies[3:5]], 
                               [m.name for m in morphologies[3:5]],            
                               xlim=900, ylim=8,   
                               colors=colors[3:5], 
                               picture_file='/tmp/distance_terminaltips_vb',   
                               picture_formats=['png','svg']
    )

    Morphology.plot_distance_distributions(
        [m.terminaltips for m in morphologies[1:3]],
        [m.root_biggest_child for m in morphologies[1:3]], 
        [m.name for m in morphologies[1:3]], 
        colors=colors[1:3], 
        bins=20, xlim=900, ylim=70, 
        picture_file='/tmp/distance_distribution_terminaltips_db',               
        picture_formats=['png','svg']
    )
    Morphology.plot_distance_distributions(
        [m.terminaltips for m in morphologies[3:5]],
        [m.root_biggest_child for m in morphologies[3:5]], 
        [m.name for m in morphologies[3:5]], 
        colors=colors[3:5], 
        bins=20, xlim=900, ylim=70, #xlim=600, ylim=40, 
        picture_file='/tmp/distance_distribution_terminaltips_vb',              
        picture_formats=['png','svg']
    )

    Morphology.plot_distance_distributions(
        [m.branching_points for m in morphologies[1:3]],
        [m.root_biggest_child for m in morphologies[1:3]], 
        [m.name for m in morphologies[1:3]], 
        colors=colors[1:3], 
        bins=20, xlim=900, ylim=65, 
        picture_file='/tmp/distance_distribution_branchpoints_db',               
        picture_formats=['png','svg']
    )
    Morphology.plot_distance_distributions(
        [m.branching_points for m in morphologies[3:5]],
        [m.root_biggest_child for m in morphologies[3:5]], 
        [m.name for m in morphologies[3:5]], 
        colors=colors[3:5], 
        bins=20, xlim=900, ylim=65, #xlim=600, ylim=40, 
        picture_file='/tmp/distance_distribution_branchpoints_vb',               
        picture_formats=['png','svg']
    )

    Morphology.plot_distance_distributions(
        [m.compartments for m in morphologies[1:3]],
        [m.root_biggest_child for m in morphologies[1:3]], 
        [m.name for m in morphologies[1:3]], 
        colors=colors[1:3], 
        bins=20, xlim=900, ylim=900, 
        picture_file='/tmp/distance_distribution_compartments_db',               
        picture_formats=['png','svg']
    )
    Morphology.plot_distance_distributions(
        [m.compartments for m in morphologies[3:5]],
        [m.root_biggest_child for m in morphologies[3:5]], 
        [m.name for m in morphologies[3:5]], 
        colors=colors[3:5], 
        bins=20, xlim=900, ylim=900, #xlim=600, ylim=700, 
        picture_file='/tmp/distance_distribution_compartments_vb',               
        picture_formats=['png','svg']
    )

    a   = [MetricAnalysis(ma).variable_map()[0] for ma in morphologies[1:5]]
    v   = {'dorsal branch':{}, 'ventral branch':{}}
    bars= ['dorsal branch','ventral branch']
    xs=['number_of_terminaltips','total_cell_length','frustum_volume','frustum_surface_area','frustum_sparsity', 'esn_sparsity']#, 'cepn_volume','cepn_surface_area','esn_frustum_volume','esn_frustum_surface_area','esn_sparsity', 'esn_cep_volume','esn_cep_surface_area','esn_cep_sparsity'
    for key in xs:#TODO: more time efficient with properties and only needed
        v['dorsal branch'][key]     = float(a[1][key]) / a[0][key] - 1
        v['ventral branch'][key]    = float(a[3][key]) / a[2][key] - 1
    xs.reverse()
    MetricAnalysis.bars_plot(
        v=v, 
        bars=bars, 
        xs=xs, 
        colors=['#00ff00','#0000ff'], 
        rotation=0, horizontal=True,#y_label='change: forager / nurse - 1', 
        picture_file='/tmp/change_morphology', picture_formats=['png','svg']
    )
