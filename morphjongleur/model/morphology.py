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
        self.neuron_h_Section.diam  = 2 * self.radius
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

    @property
    def parent_distance(self):
        """
        mapped & property
        """
        if not vars(self).has_key('_parent_distance'):
            import numpy
            if self.parent == None:
                self._parent_distance   = 0
            else:
                self._parent_distance   =  numpy.sqrt(
                      ((self.parent.x - self.x) ** 2)
                 +    ((self.parent.y - self.y) ** 2)
                 +    ((self.parent.z - self.z) ** 2)
                 )
        return self._parent_distance

    def path_distance(self, compartment):
        a = self.lca(compartment);
        left    = self
        right   = compartment
        dist = 0;
        while(left  != a):
            dist    +=left.parent_distance
            left    = left.parent
        while(right != a):
            dist    +=right.parent_distance
            right   = right.parent
        return dist

    def plot(self, center={'x':0,'y':0,'z':0}, name='', xlim=None, ylim=None, color='#000000', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        Compartment.plot_distance([self], center=center, name=name, xlim=xlim, ylim=xlim, color=color, picture_file=picture_file, picture_formats=picture_formats)

    @staticmethod
    def plot_distance(compartment_iterable, center, name='', xlim=None, ylim=None, color='#000000', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import scipy.stats.stats
        import matplotlib.pyplot
        radii   = [c.radius for c in compartment_iterable]
        distances   = [c.path_distance(center) for c in compartment_iterable]
        matplotlib.pyplot.scatter(distances, radii, c=color, marker='.', edgecolors=color)

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
        
        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format)
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
                c.length  = ( (c.x - c.parent.x)**2 + (c.y - c.parent.y)**2 + (c.z - c.parent.z)**2 )**0.5   #TODO: use .info from DB?
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
        first_child.neuron_h_Section = neuron.h.Section()
        first_child.neuron_h_Section.L     = first_child.length
        first_child.neuron_h_Section.diam  = 2 * first_child.radius
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
    def terminaltips(self):
        if not vars(self).has_key('_leafs'):
            self._leafs = []
            #self._branching_points = []
        if self._leafs == []:
            self._create_tree()
            for c in self.compartments:
                if len(c.children) == 0:
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
    def branching_points(self):
        if not vars(self).has_key('_branching_points'):
            self._branching_points = []
            #self._leafs = []
        if self._branching_points == []:
            self._create_tree()
            for c in self.compartments:
#                if len(c.children) == 0:
#                    self._leafs.append(c)
                if len(c.children) > 1:
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
        cs  = []
        for compartment in self.compartments:
            cs.append([compartment.x, compartment.y, compartment.z])
        assert len(self.compartments) == len(cs)
        pca_cs  = mdp.pca( numpy.array(cs) )
        assert len(self.compartments) == len(pca_cs)
        compartments    = []
        for i in range(len(pca_cs)):
            compartments.append(    # m.add_compartment(
                Compartment( 
                    self.compartments[i].compartment_id, 
                    self.compartments[i].compartment_parent_id, 
                    self.compartments[i].radius, 
                    x=pca_cs[i][0],
                    y=pca_cs[i][1],
                    z=pca_cs[i][2]
                ) 
            )
        return Morphology(self.name, self.file_origin, self.description, self.datetime_recording, compartments=compartments)

    def plot(self, x=0, y=1, color='#000000', picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        import matplotlib.pyplot    # 00ff00 00800 00ff80 008080
        xs  = [] 
        ys  = []
        ss  = []
        for c in self.compartments:
            xyz = (c.x, c.y, c.z)
            xs.append(xyz[x])
            ys.append(xyz[y])
            ss.append(c.radius)

        matplotlib.pyplot.scatter(xs, ys, s=ss, c=color, marker='.', edgecolors=color)#'. o
        matplotlib.pyplot.axes().set_aspect('equal', 'datalim')
        
        if(picture_file != None):
            for picture_format in picture_formats:
                try:
                    matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format)
                except Exception, e:
                    import traceback
                    print picture_format 
                    print traceback.format_exc()
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()
             
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
        import neuron
        #Morphology.__init__(self);
        # soma (compartment: neuron.h.Section() )
        self.soma = Compartment(1,-1, radius=soma_diameter/2.0);
        self.medial_dendrites = [];
        self.lateral_dendrites = [];
        for k in range(int(medials)): # medial dendrites
            medial_dendrite = Compartment(2*k+2, 1);
            self.medial_dendrites.append(medial_dendrite)
        assert len(self.medial_dendrites) == medials
        for k in range(int(laterals)): # lateral dendrites
            lateral_dendrite = Compartment(2*k+3, 1);
            self.lateral_dendrites.append(lateral_dendrite)
        assert len(self.lateral_dendrites) == laterals

        #http://code.activestate.com/recipes/52235-calling-a-superclasss-implementation-of-a-method/
        
        self.compartments  = []
        self.compartments.extend(self.medial_dendrites)
        self.compartments.extend(self.lateral_dendrites)
        self.soma.neuron_h_Section = neuron.h.Section()

        for c in self.compartments:
            c.neuron_h_Section = neuron.h.Section()
            c.neuron_h_Section.L = dendrite_length;
            c.neuron_h_Section.diam = dendrite_diameter
            c.neuron_h_Section.nseg = dendrite_diameter
            c.neuron_h_Section.Ra = dendrite_Ra

        self.compartments.append(self.soma)
            
        for medial_dendrite in self.medial_dendrites:
            medial_dendrite.neuron_h_Section.connect(self.soma.neuron_h_Section, 1, 0) # connect soma(1) with medial_dendrite(0)
        for lateral_dendrite in self.lateral_dendrites:
            lateral_dendrite.neuron_h_Section.connect(self.soma.neuron_h_Section, 0, 0) # connect soma(0) with lateral_dendrites(0)
            
        self.init_helper()

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
