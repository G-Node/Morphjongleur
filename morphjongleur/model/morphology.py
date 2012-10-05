# -*- coding: utf-8 -*-
'''
:Author: stransky
'''
import morphjongleur.util.auto_string
import numpy

class Compartment(object):
    '''
    classdocs
    @see http://web.mit.edu/neuron_v7.1/doc/help/neuron/neuron/classes/python.html#Section
    '''
    def __init__(self, compartment_id, compartment_parent_id, radius=1.0, x=float('nan'),y=float('nan'),z=float('nan'), morphology=None):#, compartment_key=None
        '''
        classdocs
        '''
        self.compartment_id = int(compartment_id)
        self.compartment_parent_id  = int(compartment_parent_id)
        self.radius         = float(radius)
        self.xyz            = numpy.array([float(x), float(y), float(z)])

        self.parent         = None  #parent automatisch mappen lassen
        self.children       = [] #necessary for Morphology._crate_tree()

        self._morphology    = morphology
        self._info          = [None]
        self._groups        = []

        self.synapse        = None

    def neuron_create(self, parent, parent_location=1, self_location=0 ):
        '''
        @see http://web.mit.edu/neuron_v7.1/doc/help/neuron/neuron/geometry.html
        '''
        import neuron
        self.neuron_h_Section = neuron.h.Section()
        #??? must be defined BEOFRE to many connections !!!
        if self.length == 0:
            self.length = numpy.finfo(self.length).tiny
            import sys
            print >> sys.stderr, "distance from %s to its parent %s = 0" %(self.__repr__(), self.parent.__repr__())#oder neu einhängen
        self.neuron_h_Section.L     = self.length
        self.neuron_h_Section.diam  = self.radius + parent.radius   #2 * self.radius
#        self.neuron_h_Section.Ra    = 
#        self.neuron_h_Section.ri    =
#TODO:        if not ( numpy.isnan(self.x) and numpy.isnan(self.y) and numpy.isnan(self.z) ) :
#            self.neuron_h_Section.x3d = self.x
#            self.neuron_h_Section.y3d = self.y
#            self.neuron_h_Section.z3d = self.z
        self.neuron_h_Section.connect( parent.neuron_h_Section, parent_location, self_location) #connect c 0 with parent(1)

    @property
    def info(self):
        if self._info[0] == None:
            self._info[0]   = Compartment_info()
        return self._info[0]

    @property
    def x(self):
        return self.xyz[0]
    @x.setter
    def x(self, x):
        self.xyz[0]  = x
    @property
    def y(self):
        return self.xyz[1]
    @y.setter
    def y(self, y):
        self.xyz[1]  = y
    @property
    def z(self):
        return self.xyz[2]
    @z.setter
    def z(self, z):
        self.xyz[2]  = z

    @staticmethod
    def huge():
        return Compartment(-7, -7, radius=float('+inf'))

    @staticmethod
    def tiny():
        return Compartment(-7, -7, radius=float('-inf'))

    def distance_euclid(self, compartment):
        return (      ((compartment.x - self.x) ** 2)
                 +    ((compartment.y - self.y) ** 2)
                 +    ((compartment.z - self.z) ** 2)
                 )**0.5

    @property
    def length(self):
        """
        parent_distance
        """
        if not vars(self).has_key('_length') or self._length == None:
            if self.parent == None:
                self._length    = 0
            else:
                self._length    =  self.distance_euclid(self.parent)
        return self._length

    def lca(self, compartment):
        '''
        lowest common ancestor
        '''
        left    = self
        right   = compartment
        leftp   = {left.compartment_id : True}
        rightp  = {right.compartment_id : True}
        while left != right:
            if(left.compartment_parent_id > 0):
                left    = left.parent
                leftp[left.compartment_id] = True
                if(rightp.get(left.compartment_id) != None):
                    return left;
            if(right.compartment_parent_id > 0):
                right   = right.parent
                rightp[right.compartment_id] = True;
                if(leftp.get(right.compartment_id) != None):
                    return right;
        return left

    def distance_path(self, compartment):
        a = self.lca(compartment);
        left    = self
        right   = compartment
        dist = 0;
        while(left  != a):
            dist    +=left.length
            left    = left.parent
        while(right != a):
            dist    +=right.length
            right   = right.parent
        return dist

    @staticmethod
    def write_svg(svg_file, compartment_iterables, colors=['#000000'], x=0, y=1):
        '''
        write Scalable Vector Graphics file
        '''
        import itertools
        import math
        compartment_iterables1  = []
        compartment_iterables2  = []
        for compartment_iterable in compartment_iterables:
            compartment_iterable1,compartment_iterable2    = itertools.tee(compartment_iterable)
            compartment_iterables1.append(compartment_iterable1) 
            compartment_iterables2.append(compartment_iterable2)
        compartment_iterables   = compartment_iterables1

        stream  = open(svg_file, "w")
        stream.write('<?xml version="1.0"?>\n')
        stream.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
        xyz_min = [float("+inf"),float("+inf"),float("+inf")]
        xyz_max = [float("-inf"),float("-inf"),float("-inf")]
        radius_max = float("-inf")
        for compartment_iterable in compartment_iterables2:
            for c in compartment_iterable:
                for d in [x,y]:#range(3):
                    if xyz_min[d] > c.xyz[d]:
                        xyz_min[d] = c.xyz[d]
                    if xyz_max[d] < c.xyz[d]:
                        xyz_max[d] = c.xyz[d]
                if radius_max < c.radius:
                    radius_max = c.radius
        xyz_shift   = [0,0,0]
        for d in [x,y]:#range(3):
            xyz_shift[d]   = -xyz_min[d]+radius_max
        width   = int( xyz_max[x]-xyz_min[x]+2*radius_max + 0.5)
        height  = int( xyz_max[y]-xyz_min[y]+2*radius_max+1 + 0.5)
        stream.write('<svg xmlns="http://www.w3.org/2000/svg" width="%i" height="%i">\n' % (width,height) )
        for i in xrange(len(compartment_iterables)):
            color   = colors[i%len(colors)]
            for c in compartment_iterables[i]:#id="%i" c.compartment_id, 
                stream.write('<circle cx="%f" cy="%f" r="%f" fill="%s" />\n' % (c.xyz[x]+xyz_shift[x], c.xyz[y]+xyz_shift[y], c.radius, color))
        bar_size    = 10**int(math.log10(width-0.1))
        stream.write('<line x1="%i" y1="%i" x2="%i" y2="%i" style="stroke:black;stroke-width:1"/>\n' % (0, height-3, bar_size, height-3))
        stream.write('<text x="%i"  y="%i">%i µm</text>\n' % (bar_size + 3,height-3,bar_size))
        stream.write('</svg>\n')
        stream.close()

    @staticmethod
    def plot(compartment_iterables, colors, x=0, y=1, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import matplotlib.pyplot    # 00ff00 00800 00ff80 008080
        import numpy
        
        for i in xrange(len(compartment_iterables)):
            xs  = [] 
            ys  = []
            ss  = []
            for c in compartment_iterables[i]:
                xyz = (c.x, c.y, c.z)
                xs.append(xyz[x])
                ys.append(xyz[y])
                ss.append(numpy.pi * c.radius**2)
                #matplotlib.pyplot.axes().add_artist(
                #    matplotlib.patches.Circle(xy=(c.x, c.y),
                #      color=colors[i%len(colors)],
                #      radius=c.radius
                #    )
                #)
            ss = numpy.diff(matplotlib.pyplot.axes().transData.transform(zip([0]*len(ss), ss))) 
            matplotlib.pyplot.scatter(xs, ys, s=ss, c=colors[i], marker='.', edgecolors=colors[i])#'. o
        matplotlib.pyplot.axes().set_aspect('equal', 'datalim')
        #fig = matplotlib.pyplot.gcf()
        #fig.set_size_inches(7,7)

        #matplotlib.pyplot.axes().yaxis.set_ticks([])
        #matplotlib.pyplot.axes().xaxis.set_ticks([])
        #for spine in matplotlib.pyplot.axes().spines.values():
        #    spine.set_visible(spines_visible)
        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, dpi=600, transparent=True)
                except Exception, e:
                    import traceback
                    print picture_format 
                    print traceback.format_exc()
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

    @staticmethod
    def plot3d(compartment_iterable, color='black', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import mpl_toolkits.mplot3d
        import matplotlib.pyplot
        import numpy
        
        fig = matplotlib.pyplot.figure()
        #for c in compartment_iterable:
            #draw sphere
            #u, v = numpy.mgrid[0:2*numpy.pi:20j, 0:numpy.pi:10j]
            #x=numpy.cos(u)*numpy.sin(v)
            #y=numpy.sin(u)*numpy.sin(v)
            #z=numpy.cos(v)
            #ax.plot_wireframe(x, y, z, color="r")
        ax = fig.add_subplot(111, projection='3d')
        for ci in [compartment_iterable]:
            xs  = [] 
            ys  = []
            zs  = []
            cs  = []
            ss  = []
            colorbar    = True
            for c in compartment_iterable:
                xs.append(c.x)
                ys.append(c.y)
                zs.append(c.z)
                if vars(c).has_key('color'):
                    cs.append(float(c.color))
                else:
                    colorbar    = False
                    cs.append(color)
                ss.append(numpy.pi * c.radius**2)
            ss = numpy.diff(ax.transData.transform(zip([0]*len(ss), ss)))
            p = ax.scatter(xs, ys, zs, c=cs, marker='.', s=ss)
        if False and colorbar:
            matplotlib.pyplot.figure().colorbar(p)
        ax.set_aspect('equal')
        #fig.set_size_inches(7,7)

        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, dpi=600, transparent=True)
                except Exception, e:
                    import traceback
                    print picture_format 
                    print traceback.format_exc()
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

    @staticmethod
    def plot_color(compartment_iterables, x=0, y=1, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        '''
        requires c.color
        '''
        import matplotlib.pyplot    # 00ff00 00800 00ff80 008080
        fig = matplotlib.pyplot.figure()
        for i in xrange(len(compartment_iterables)):
            xs  = [] 
            ys  = []
            ss  = []
            cs  = []
            for c in compartment_iterables[i]:
                xyz = (c.x, c.y, c.z)
                xs.append(xyz[x])
                ys.append(xyz[y])
                ss.append(c.radius)
                cs.append(float(c.color))
                matplotlib.pyplot.axes().add_artist(
                    matplotlib.patches.Circle(
                        xy=(c.x, c.y),
                        radius=c.radius,
                        color=c.color
                    )
                )
            #ss = numpy.diff(matplotlib.pyplot.axes().transData.transform(zip([0]*len(ss), ss))) 
            p = matplotlib.pyplot.scatter(xs, ys, s=ss, c=cs, marker='.', edgecolors=cs)#'. o
        #fig.colorbar(p)
        matplotlib.pyplot.axes().set_aspect('equal', 'datalim')
        #fig = matplotlib.pyplot.gcf()
        #fig.set_size_inches(7,7)

        #matplotlib.pyplot.axes().yaxis.set_ticks([])
        #matplotlib.pyplot.axes().xaxis.set_ticks([])
        #for spine in matplotlib.pyplot.axes().spines.values():
        #    spine.set_visible(spines_visible)
        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=True)#, dpi=600
                except Exception, e:
                    import traceback
                    print picture_format 
                    print traceback.format_exc()
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

    def plot_distance(self, compartment_iterable, name='', xlim=None, ylim=None, color='#000000', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        Compartment.plot_distances([compartment_iterable], [self], [name], xlim, ylim, [color], picture_file, picture_formats)

    @staticmethod
    def plot_distances(compartment_iterables, centers, names, xlim=None, ylim=None, colors=['#000000'],  picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import scipy.stats.stats
        import matplotlib.pyplot
        legends = []
        for (compartment_iterable,center,color,name) in zip(compartment_iterables,centers,colors,names):
            import itertools
            compartment_iterable, compartment_iterable2 = itertools.tee(compartment_iterable)
            radii   = [c.radius for c in compartment_iterable]
            distances   = [c.distance_path(center) for c in compartment_iterable2]
            matplotlib.pyplot.scatter(distances, radii, c=color, marker='.', edgecolors=color)
            legends.append("%s" % (name))

            if False:
                r_mean  = numpy.mean(radii)
                r_std   = numpy.std(radii)
                r_gmean = scipy.stats.stats.gmean(radii)
                r_hmean = scipy.stats.stats.hmean(radii)
                r_median= numpy.median(radii)
                
                matplotlib.pyplot.axhline(y=r_mean,     color='blue')
                xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis()
                width    = r_std / (ymax-ymin)
                center    = (r_mean - ymin) / (ymax-ymin)
                matplotlib.pyplot.axvline(x=.5*(xmax-xmin), ymin=center-.5*width, ymax=center+.5*width, color='red')
                matplotlib.pyplot.axhline(y=r_gmean,    color='violet')
                matplotlib.pyplot.axhline(y=r_hmean,    color='orange')
                matplotlib.pyplot.axhline(y=r_median,   color='black')
            
                d_mean  = numpy.mean(distances)
                d_std   = numpy.std(distances)
                if numpy.min(distances) <= 0:
                    d_gmean = 0
                    d_hmean = 0
                else:
                    d_gmean = scipy.stats.stats.gmean(distances)
                    d_hmean = scipy.stats.stats.hmean(distances)
                d_median= numpy.median(distances)
                
                matplotlib.pyplot.axvline(x=d_mean,     color='blue')
                xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis()
                width    = d_std / (xmax-xmin)
                center    = (d_mean - xmin) / (xmax-xmin)
                matplotlib.pyplot.axhline(y=.5*(xmax-xmin), xmin=center-.5*width, xmax=center+.5*width, color='red')
                matplotlib.pyplot.axvline(x=d_gmean,    color='violet')
                matplotlib.pyplot.axvline(x=d_hmean,    color='orange')
                matplotlib.pyplot.axvline(x=d_median,   color='black')
    
                matplotlib.pyplot.legend( [ 
                    'r mean %f' % ( r_mean ), 
                    'r std %f' % ( r_std ), 
                    'r gmean %f' % ( r_gmean ),
                    'r hmean %f' % ( r_hmean ),
                    'r median %f' % ( r_median ),
                    'd mean %f' % ( d_mean ), 
                    'd std %f' % ( d_std ), 
                    'd gmean %f' % ( d_gmean ),
                    'd hmean %f' % ( d_hmean ),
                    'd median %f' % ( d_median ),
                    'compartments'
                ] )

        matplotlib.pyplot.ylabel(u'radius [µm]')
        if ylim != None:
            matplotlib.pyplot.ylim((0,ylim))
        matplotlib.pyplot.xlabel(u'distance [µm]')
        if xlim != None:
            matplotlib.pyplot.xlim((0,xlim))
        matplotlib.pyplot.legend( legends )
        
        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=True)
                except Exception, e:
                    print picture_format 
                    print traceback.format_exc()
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()
    def __repr__(self):
        return "Compartment(%i, %i, %f, %f,%f,%f)" % (
            self.compartment_id, self.compartment_parent_id, self.radius, self.x, self.y, self.z)
    def __str__(self):
        return """<Compartment(
 compartment_key        = '%s' 
 compartment_id         = %i
 compartment_parent_id  = %i
 radius                 = %f
%s groups = [
%s          ]
)>""" % (
            str(self.compartment_key if vars(self).has_key('compartment_key') else ''), 
            int(self.compartment_id), 
            int(self.compartment_parent_id), 
            float(self.radius),
            str( self.info if self.info != None else ''),
            "           ,\n".join(map(str, self._groups)),
            )


class Morphology(object):
    '''
    Model of a simulatable Neuron
    
    To check the morphology with NEURON gui:
    >>> from neuron import gui
    '''

    _type_compartment   = Compartment

    def __init__(self, name, file_origin, description, datetime_recording, compartments=[]):
        '''
        Constructor
        '''
        #assert len(compartments) == 0
        #self.morphology_key     = morphology_key
        self.name               = name
        self.file_origin        = file_origin
        self.description        = description
        self.datetime_recording = datetime_recording
        if compartments == []:#or list(compartments) to prevent static list
            self._compartments  = []
        else:
            self._compartments  = compartments
        self._info              = [None]
        self._groups            = []

        #TODO: map orm
        self._root              = None
        self._biggest           = None  # biggest compartment (probable part of soma)
        self._leafs             = []
        self._branching_points  = []
        #not orm mapped
        self._compartments_map  = {}

    @property
    def info(self):
        if self._info[0] == None:
            self._info[0]   = Morphology_info()
        return self._info[0]

    def _create_tree(self):
        if vars(self).has_key('_compartments_map') and self._compartments_map != {}:
            return
        self._compartments_map  = {}
        self._biggest   = self._compartments[0] if self._compartments[0].compartment_parent_id != -1 else self._compartments[1]

        for c in self.compartments:
            self._compartments_map[ c.compartment_id ] = c
            c.children  = []
            if(c.radius > self._biggest.radius and c.compartment_parent_id != -1):
                self._biggest   = c
        
        for c in self.compartments:
            assert c.compartment_parent_id != "-1"
            if c.compartment_parent_id != -1:
                c.parent  = self._compartments_map[ c.compartment_parent_id ]
                c.parent.children.append( c )
            else:
                if vars(self).has_key('_root') and self._root != None:
                    import sys
                    print >> sys.stderr, "multiple roots: self.root %s != %s" % (self.root , c)
                self._root = c
                
    def neuron_create(self):
        import neuron
        if vars(self.root.children[0]).has_key('neuron_h_Section'):
            return

        todo_stack    = []
        root    = self.root
        first_child = root.children[0]
        root.neuron_h_Section = neuron.h.Section()
        first_child.neuron_h_Section = root.neuron_h_Section
        root.neuron_h_Section.L     = first_child.length
        root.neuron_h_Section.diam  = root.radius + first_child.radius  # 2 * first_child.radius
        todo_stack.extend( first_child.children )
        for c in root.children[1:]:
            c.neuron_create( first_child, parent_location=0, self_location=0 ) #connect c 0 with parent(0)
            todo_stack.extend( c.children )
        
        #deepth-first
        while len( todo_stack ) > 0:
            c   = todo_stack.pop()
            c.neuron_create( c.parent, parent_location=1, self_location=0 ) #connect c 0 with parent(1)
            todo_stack.extend( c.children )

    @property
    def compartments_map(self):
        if not vars(self).has_key('_compartments_map') or self._compartments_map == {}:
            self._create_tree()
        return self._compartments_map

    @property
    def root(self):
        if not vars(self).has_key('_root') or self._root == None:
            self._create_tree()
        return self._root

    @property
    def biggest(self):
        if not vars(self).has_key('_biggest') or self._biggest == None:
            self._create_tree()
        return self._biggest

    @property
    def root_biggest_child(self):
        if not vars(self).has_key('_root_biggest_child') or self._root_biggest_child == None:
            self._create_tree()
            self._root_biggest_child    = self._root.children[0]
            for c in self._root.children:
                if(c.radius > self._root_biggest_child):
                    self._root_biggest_child   = c
        return self._root_biggest_child

    @property
    def terminaltips(self):
        if not vars(self).has_key('_leafs'):
            self._leafs = []
            #self._branching_points = []
        if self._leafs == []:
            self._create_tree()
            for c in self.compartments:
                if len(c.children) == 0 or len(c.children) == 1 and c.parent == None:
                    self._leafs.append(c)
#                elif len(c.children) > 1:
#                    self.branching_points.append(c)
        for leaf in self._leafs:
            yield(leaf)

    @property
    def number_of_terminaltips(self):
        if not vars(self).has_key('_leafs') or self._leafs == []:
            for l in self.terminaltips:
                pass
        return len(self._leafs)

    @property
    def terminaltips_biggest(self):
        tb  = Compartment.tiny()
        for t in self.terminaltips:
            if tb.radius < t.radius:
                tb = t
        return tb

    @property
    def branching_points(self):
        if not vars(self).has_key('_branching_points'):
            self._branching_points = []
            #self._leafs = []
        if self._branching_points == []:
            self._create_tree()
            for c in self.compartments:
#                if len(c.children) == 0:
#                    self._leafs.append(c)
                if len(c.children) > 2 or len(c.children) > 1 and c.parent != None:
                    self._branching_points.append(c)

        for branching_point in self._branching_points:
            yield(branching_point)

    @property
    def number_of_branching_points(self):
        if not vars(self).has_key('_branching_points') or self._branching_points == []:
            for b in self.branching_points:
                pass
        return len(self._branching_points)

    @property
    def plebs(self):
        for branching_point in self._branching_points:
            if len(branching_point.children) > 2:
                yield(branching_point)

    @property
    def compartments(self):
        '''
        generator over compartments
        '''
        for compartment in self._compartments:
            yield(compartment)

    @property
    def number_of_compartments(self):
        return len(self._compartments)

    @property
    def non_root_compartments(self):
        '''
        generator over compartments
        '''
        if not vars(self).has_key('_root') or self._root == None:
            self._create_tree()
        for compartment in self.compartments:
            if compartment.parent != None:
                yield(compartment)

    def get_compartment(self, compartment_id):
        '''
        in order to be conform with Model definition, it is possible to get an compartment by i.
        '''
        return self.compartments_map[ compartment_id ];

    def add_compartment(self, compartment):
        """
        doc-string to be added.
        """
        self._compartments.append(compartment)
        #compartment = Compartment(parent, id, length)
        #parent.children[id] = compartment

    def subtree(self, compartment):
        if not vars(self).has_key('_compartments_map') or self._compartments_map == {}:
            self._create_tree()
        if type(compartment) == type(1):
            compartment = self.get_compartment(compartment)

        todo_stack    = []
        yield( compartment )
        todo_stack.append( compartment )
        while len( todo_stack ) > 0:
            c   = todo_stack.pop()
            todo_stack.extend( c.children )
            for cc in c.children:
                yield( cc )
    
    def pca(self):
        import mdp
        pca_cs  = mdp.pca( numpy.array([ [c.x, c.y, c.z] for c in self.compartments ] ) )
        assert self.number_of_compartments == len(pca_cs)
        compartments    = []
        for i in xrange( self.number_of_compartments ):
            compartments.append(    # m.add_compartment(
                Compartment( 
                    self._compartments[i].compartment_id, 
                    self._compartments[i].compartment_parent_id, 
                    self._compartments[i].radius, 
                    x=pca_cs[i][0],
                    y=pca_cs[i][1],
                    z=pca_cs[i][2]
                ) 
            )
        return Morphology(self.name, self.file_origin, self.description, self.datetime_recording, compartments=compartments)

    def plot_distance_distribution(self, center, name='', xlim=None, ylim=None, color='#000000', bins=20, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
        Morphology.plot_distance_distributions([self.compartments], [center], [name], [color], bins, xlim, ylim, picture_file, picture_formats)

    @staticmethod
    def plot_distance_distributions(compartment_iterables, centers, names=[''], colors=['#000000'], bins=20, xlim=None, ylim=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
        import matplotlib
        #matplotlib.rc('text', usetex=True): error with names
        legends = []
        for (compartment_iterable,center,color,name) in zip(compartment_iterables,centers,colors,names):
            x   = [c.distance_path(center) for c in compartment_iterable]
            if len(x) == 0:
                import sys
                print >> sys.stderr, "iterable list has 0 elements"
                continue
            mean   = numpy.mean(x)
            legends.append(u"%7s: %i µm" % (name, round(mean)))
            std    = numpy.std(x)
            matplotlib.pyplot.axvline(x=mean, color=color, label='mean'+name)
            if xlim != None:#TODO: isnumber
                matplotlib.pyplot.hist(x, bins=range(0,xlim,bins), normed=0, color=color, edgecolor=color, alpha=0.6)
            else:
                matplotlib.pyplot.hist(x, bins, normed=0, color=color, edgecolor=color, alpha=0.6)

        #matplotlib.pyplot.title('Endpoints of %s' % (name.replace('_',' ')) )
        #print 'distribution of %s : mean=%f, std=%f' % (name, mean, std)


        matplotlib.pyplot.grid(True, color='lightgrey')
        matplotlib.pyplot.ylabel('#')#%
        if ylim != None:
            matplotlib.pyplot.ylim((0,ylim))
        matplotlib.pyplot.xlabel(u'distance [µm]')
        if xlim != None:
            matplotlib.pyplot.xlim((0,xlim))
        matplotlib.pyplot.legend( legends )

        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format, format=picture_format, transparent=True)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close('all')

    def plot(self, x=0, y=1, color='#000000', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        Compartment.plot([self.compartments], colors=[color], x=x, y=y, picture_file=picture_file, picture_formats=picture_formats)

    def write_svg(self, svg_file, color='#000000', x=0, y=1):
        '''
        write Scalable Vector Graphics file
        '''
        Compartment.write_svg(svg_file=svg_file, compartment_iterables=[self.compartments], colors=[color], x=x, y=y)

    def __repr__(self):
        return "Morphology('%s', '%s', '%s', '%s')" % ( 
            str(self.name), str(self.file_origin), str(self.description), str(self.datetime_recording) )

    def __str__(self):
        return """<Morphology(
 morphology_key         = '%s' 
 name                   = '%s' 
 file_origin            = '%s' 
 description            = '%s' 
 datetime_recording     = '%s'
%s groups = [
%s          ]
%s
)>""" % (
            str(self.morphology_key if vars(self).has_key('morphology_key') else ''), 
            str(self.name), 
            str(self.file_origin), 
            str(self.description), 
            str(self.datetime_recording),
            str( self.info if self.info != None else ''),
#            "           ,\n".join(map(str, self._groups)) if self.__dict__.has_key('_groups') else '',
            "           ,\n".join(map(str, self._groups)),
            str( self.analysis if self.__dict__.has_key('analysis') else '')
        )

class Star(Morphology):
    '''
    n-dendrites = n+1-Compartments Model
     
    
    This class will produce Neuron with simple star formation 
    with a standard soma (L=40 um, diam=20 um) with identical dendrites connected on opposite sites of the
    soma. 
    
    For the dendrites the following parameters can be changed:
    * dendrite_Ra:   
    * dendrite_length:   length of each dendrite
    * dendrite_diameter: diameter of each dendrite
    * dendrite_nseg: 
    '''
 
    def __init__(self, medials=1, laterals=1, use_passive_channes=True, gp=0.004, E=-60, Ra=200,
        soma_Ra=1, soma_length=40, soma_diameter=20, soma_nseg=10, 
        dendrite_Ra=1, dendrite_length=150, dendrite_diameter=3, dendrite_nseg=int(150/10)):
        '''
        consisting of `medials` medial dendrite and `laterals` lateral dendrite
        
        soma_Ra: Axial Resistivity (Ra): 200 [Ohm * cm]
        soma_L in µm; stored as a float number
        soma_diam diameter (soma.diam): [µm]
        soma_nseg: stored as an integer
        '''
        self.name   = "Star_%i_%i" % (medials, laterals)
        
        import neuron
        #Morphology.__init__(self);
        # soma (compartment: neuron.h.Section() )
        self.soma = Compartment(1,-1, radius=soma_diameter/2.0);
        self.medial_dendrites = [];
        self.lateral_dendrites = [];
        for k in xrange(int(medials)): # medial dendrites
            self.medial_dendrites.append( Compartment(2*k+2, 1) )
        assert len(self.medial_dendrites) == medials
        for k in xrange(int(laterals)): # lateral dendrites
            self.lateral_dendrites.append( Compartment(2*k+3, 1) )
        assert len(self.lateral_dendrites) == laterals

        #http://code.activestate.com/recipes/52235-calling-a-superclasss-implementation-of-a-method/
        
        self._compartments  = []
        self._compartments.extend(self.medial_dendrites)
        self._compartments.extend(self.lateral_dendrites)
        self.soma.neuron_h_Section = neuron.h.Section()

        for c in self.compartments:
            c.neuron_h_Section = neuron.h.Section()
            c.neuron_h_Section.L = dendrite_length;
            c.neuron_h_Section.diam = dendrite_diameter
            c.neuron_h_Section.nseg = dendrite_nseg
            c.neuron_h_Section.Ra = dendrite_Ra

        self._compartments.append(self.soma)
            
        for medial_dendrite in self.medial_dendrites:
            medial_dendrite.neuron_h_Section.connect(self.soma.neuron_h_Section, 1, 0) # connect soma(1) with medial_dendrite(0)
        for lateral_dendrite in self.lateral_dendrites:
            lateral_dendrite.neuron_h_Section.connect(self.soma.neuron_h_Section, 0, 0) # connect soma(0) with lateral_dendrites(0)


    def get_compartment(self, compartment_id):
        '''
        in order to be conform with Model definition, it is possible to get an compartment by i.
        compartment_id = 0 returns soma
        compartment_id =-1 returns lateral dendrite
        compartment_id = 1 returns medial dendrite
        '''
        if(compartment_id < 0):
            return self.lateral_dendrites[0];
        if(compartment_id > 0):
            return self.medial_dendrites[0];
        if(compartment_id == 0):
            return self.soma;       
    
    def __repr__(self):
        return 'Star(medials=%i, laterals=%i)' % ( #, use_passive_channes=%s, gp=%f, E=%f, Ra=%f, soma_Ra=%f, soma_L=%f, soma_diam=%f, soma_nseg=%f, dendrite_Ra=%f, dendrite_length=%f, dendrite_diameter=%f, dendrite_nseg=%i)' % (
            len(self.medial_dendrites), len(self.lateral_dendrites)
        );
        
    def __str__(self):
        return '<Star has %i medial dendrites and %i lateral dendrites>' % (
           len(self.medial_dendrites), len(self.lateral_dendrites)
        );




@morphjongleur.util.auto_string.auto_string
class MorphologyGroups(object):
    pass

@morphjongleur.util.auto_string.auto_string
class CompartmentGroups(object):
    pass

@morphjongleur.util.auto_string.auto_string
class MorphologyInfo(object):
    pass

@morphjongleur.util.auto_string.auto_string
class Morphology_info(object):
    """
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
    """

    import morphjongleur.util.metric_analysis

    @property
    def path_length(self):
        if not vars(self).has_key('_path_length') or self._path_length == None:
            self._path_length   =  morphjongleur.util.metric_analysis.path_length(self.morphology)
        return self._path_length
    @path_length.setter
    def path_length(self, value):raise AttributeError("cannot change calculated information")
    @path_length.deleter
    def path_length(self):       del self._path_length

    @property
    def surface_length_frustum(self):
        if not vars(self).has_key('_surface_length') or self._surface_length == None:
            self._surface_length   = morphjongleur.util.metric_analysis.surface_length_frustum(self.morphology)
        return self._surface_length
    @surface_length_frustum.setter
    def surface_length_frustum(self, value):raise AttributeError("cannot change calculated information")
    @surface_length_frustum.deleter
    def surface_length_frustum(self):       del self._surface_length

    @property
    def cylindric_volume(self):
        if not vars(self).has_key('_cylindric_volume') or self._cylindric_volume == None:
            self._cylindric_volume   = morphjongleur.util.metric_analysis.cylindric_volume(self.morphology)
        return self._cylindric_volume
    @cylindric_volume.setter
    def cylindric_volume(self, value):raise AttributeError("cannot change calculated information")
    @cylindric_volume.deleter
    def cylindric_volume(self):       del self._cylindric_volume

    @property
    def frustum_volume(self):
        if not vars(self).has_key('_frustum_volume') or self._frustum_volume == None:
            self._frustum_volume   = morphjongleur.util.metric_analysis.frustum_volume(self.morphology)
        return self._frustum_volume
    @frustum_volume.setter
    def frustum_volume(self, value):raise AttributeError("cannot change calculated information")
    @frustum_volume.deleter
    def frustum_volume(self):       del self._frustum_volume

    @property
    def cylindric_lateral_area(self):
        if not vars(self).has_key('_cylindric_lateral_area') or self._cylindric_lateral_area == None:
            self._cylindric_lateral_area   = morphjongleur.util.metric_analysis.cylindric_lateral_area(self.morphology)
        return self._cylindric_lateral_area
    @cylindric_lateral_area.setter
    def cylindric_lateral_area(self, value):raise AttributeError("cannot change calculated information")
    @cylindric_lateral_area.deleter
    def cylindric_lateral_area(self):       del self._cylindric_lateral_area

    @property
    def frustum_lateral_area(self):
        if not vars(self).has_key('_frustum_lateral_area') or self._frustum_lateral_area == None:
            self._frustum_lateral_area   = morphjongleur.util.metric_analysis.frustum_lateral_area(self.morphology)
        return self._frustum_lateral_area
    @frustum_lateral_area.setter
    def frustum_lateral_area(self, value):raise AttributeError("cannot change calculated information")
    @frustum_lateral_area.deleter
    def frustum_lateral_area(self):       del self._frustum_lateral_area

    @property
    def cylindric_surface_area(self):
        if not vars(self).has_key('_cylindric_surface_area') or self._cylindric_surface_area == None:
            self._cylindric_surface_area   = morphjongleur.util.metric_analysis.cylindric_surface_area(self.morphology)
        return self._cylindric_surface_area
    @cylindric_surface_area.setter
    def cylindric_surface_area(self, value):raise AttributeError("cannot change calculated information")
    @cylindric_surface_area.deleter
    def cylindric_surface_area(self):       del self._cylindric_surface_area

    @property
    def frustum_surface_area(self):
        if not vars(self).has_key('_frustum_surface_area') or self._frustum_surface_area == None:
            self._frustum_surface_area   = morphjongleur.util.metric_analysis.frustum_surface_area(self.morphology)
        return self._frustum_surface_area
    @frustum_surface_area.setter
    def frustum_surface_area(self, value):raise AttributeError("cannot change calculated information")
    @frustum_surface_area.deleter
    def frustum_surface_area(self):       del self._frustum_surface_area

    @property
    def branches(self):
        if not vars(self).has_key('_branches') or self._branches == None:
            self._branches   = morphjongleur.util.metric_analysis.branching_points(self.morphology)
        return self._branches
    @branches.setter
    def branches(self, value):raise AttributeError("cannot change calculated information")
    @branches.deleter
    def branches(self):       del self._branches


@morphjongleur.util.auto_string.auto_string
class CompartmentInfo(object):
    pass

@morphjongleur.util.auto_string.auto_string
class Compartment_info(object):
    """
 parent_radius          = %f,
 length                 = %f,
 cylindric_volume       = %f,
   frustum_volume       = %f,
 cylindric_lateral_area = %f,
   frustum_lateral_area = %f,
   frustum_length       = %f,
 #children              = %i
    """

    @property
    def parent_radius(self):
        if not vars(self).has_key('_parent_radius') or self._parent_radius == None:
            self._parent_radius   = float('nan')
        return self._parent_radius
    @parent_radius.setter
    def parent_radius(self, value):raise AttributeError("cannot change calculated information")
    @parent_radius.deleter
    def parent_radius(self):       del self._parent_radius

    @property
    def length(self):
        if not vars(self).has_key('_length') or self._length == None:
            self._length   = float('nan')
        return self._length
    @length.setter
    def length(self, value):raise AttributeError("cannot change calculated information")
    @length.deleter
    def length(self):       del self._length

    @property
    def cylindric_volume(self):
        if not vars(self).has_key('_cylindric_volume') or self._cylindric_volume == None:
            self._cylindric_volume   = float('nan')
        return self._cylindric_volume
    @cylindric_volume.setter
    def cylindric_volume(self, value):raise AttributeError("cannot change calculated information")
    @cylindric_volume.deleter
    def cylindric_volume(self):       del self._cylindric_volume

    @property
    def frustum_volume(self):
        if not vars(self).has_key('_frustum_volume') or self._frustum_volume == None:
            self._frustum_volume   = float('nan')
        return self._frustum_volume
    @frustum_volume.setter
    def frustum_volume(self, value):raise AttributeError("cannot change calculated information")
    @frustum_volume.deleter
    def frustum_volume(self):       del self._frustum_volume

    @property
    def cylindric_lateral_area(self):
        if not vars(self).has_key('_cylindric_lateral_area') or self._cylindric_lateral_area == None:
            self._cylindric_lateral_area   = float('nan')
        return self._cylindric_lateral_area
    @cylindric_lateral_area.setter
    def cylindric_lateral_area(self, value):raise AttributeError("cannot change calculated information")
    @cylindric_lateral_area.deleter
    def cylindric_lateral_area(self):       del self._cylindric_lateral_area

    @property
    def frustum_lateral_area(self):
        if not vars(self).has_key('_frustum_lateral_area') or self._frustum_lateral_area == None:
            self._frustum_lateral_area   = float('nan')
        return self._frustum_lateral_area
    @frustum_lateral_area.setter
    def frustum_lateral_area(self, value):raise AttributeError("cannot change calculated information")
    @frustum_lateral_area.deleter
    def frustum_lateral_area(self):       del self._frustum_lateral_area

    @property
    def frustum_length(self):
        if not vars(self).has_key('_frustum_length') or self._frustum_length == None:
            self._frustum_length   = float('nan')
        return self._frustum_length
    @frustum_length.setter
    def frustum_length(self, value):raise AttributeError("cannot change calculated information")
    @frustum_length.deleter
    def frustum_length(self):       del self._frustum_length

    @property
    def children(self):
        if not vars(self).has_key('_children') or self._children == None:
            self._children   = float('nan')
        return self._children
    @children.setter
    def children(self, value):raise AttributeError("cannot change calculated information")
    @children.deleter
    def children(self):       del self._children


if __name__ == '__main__':
#    import morphjongleur.io.swc
    morphologies    = []
#    for swc in ['../../data/test.swc','../../data/H060602DB_10_2_zentai_.swc','../../data/H060602VB_10_2_zentai_.swc','../../data/H060607DB_10_2(zentai).swc','../../data/H060607VB_10_2(zentai).swc']:
#        print swc
#        morphology   = morphjongleur.model.morphology.Morphology.swc_parse(swc)

    from morphjongleur.io.database import Database
    import morphjongleur.orm.morphology
    db = Database(
                db_name='postgresql://hal08.g-node.pri/morphjongleur',
                exec_role='morphjokey_admin',#TODO: select does not find table!
                exec_path='mitsubachi'
            )
    # must be mapped before Object is created
    mapping = morphjongleur.orm.morphology.Mapper( db )
    mapping.orm_map()
    for i in [3,5,4,6]:#[H060602DB_10_2_zentai_','H060607DB_10_2(zentai)','H060602VB_10_2_zentai_','H060607VB_10_2(zentai)']
        morphology   = mapping.load_morphology(i)
        morphologies.append(morphology)
        morphology.plot_endpoints_histogramm(xlim=(0, 1000), ylim=(0, 100), picture_file='/tmp/mitsubachi_endpoints_'+str(morphology.name), picture_formats=['svg'])#
    morphology.plot_all_properties(morphologies=morphologies, picture_file='/tmp/mitsubachi_', picture_formats=['svg'])
