# -*- coding: utf-8 -*-
'''
:Author: stransky
'''
import morphjongleur.util.auto_string


class Compartment(object):
    '''
    classdocs
    @see http://web.mit.edu/neuron_v7.1/doc/help/neuron/neuron/classes/python.html#Section
    '''
    def __init__(self, compartment_id, compartment_parent_id, radius=1.0, x=float('nan'),y=float('nan'),z=float('nan'), morphology=None):#, compartment_key=None
        '''
        classdocs
        '''
        #self.compartment_key   = compartment_key
        #self.morphology_key    = morphology_key
        self.compartment_id     = int(compartment_id)
        self.compartment_parent_id  = int(compartment_parent_id)
        self.radius  = float(radius)
        self.x  = float(x)
        self.y  = float(y)
        self.z  = float(z)
        self.parent         = None  #parent automatisch mappen lassen
        self._morphology    = morphology
        self._info          = [None]
        self._groups        = []
        self.children       = [] #necessary for Morphology._crate_tree()
        self.synapse        = None
        
    @property
    def info(self):
        if self._info[0] == None:
            self._info[0]   = Compartment_info
        return self._info[0]
    @info.setter
    def info(self, value):  raise AttributeError("cannot change calculated information")
    @info.deleter
    def info(self):         del self._info
    #info = property(get_info, set_info, del_info, "I'm the 'dict' property, containing calculated properties.")
        
    @property
    def parent_distance(self):
        """
        mapped & property
        """
        if not vars(self).has_key('_parent_distance'):
            import numpy
            self._parent_distance   =  numpy.sqrt(
                  ((self.parent.x - self.x) ** 2)
             +    ((self.parent.y - self.y) ** 2)
             +    ((self.parent.z - self.z) ** 2)
             )
        return self._parent_distance
    @parent_distance.setter
    def parent_distance(self, value):  raise AttributeError("cannot change calculated information")
    @parent_distance.deleter
    def parent_distance(self):         del self._parent_distance
    #parent_distance= property(get_info, set_info, del_info, "I'm the 'dict' property, containing calculated properties.")

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
        return dist;

    def lca(self, compartment):
        left    = self
        right   = compartment
        leftp = {};
        rightp = {};
        while True:
            if(left.compartment_parent_id > 0):
                left=left.parent
                leftp[left.compartment_id] = True
                if(rightp.get(left.compartment_id) != None):
                    return left;
            if(right.compartment_parent_id > 0):
                right=right.parent
                rightp[right.compartment_id] = True;
                if(leftp.get(right.compartment_id) != None):
                    return right;

    def neuron_create(self, parent, parent_location=1, self_location=0 ):
        '''
        @see http://web.mit.edu/neuron_v7.1/doc/help/neuron/neuron/geometry.html
        '''
        import neuron
        import numpy
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

    def __init__(self, name, file_origin, description, datetime_recording, compartments=None):
        print compartments
        '''
        Constructor
        '''
        #self.morphology_key     = morphology_key
        self.name               = name
        self.file_origin        = file_origin
        self.description        = description
        self.datetime_recording = datetime_recording
        if compartments == None:#or list(compartments) to prevent static list
            self.compartments   = []
        else:
            self.compartments   = compartments
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
    @info.setter
    def info(self, value):  raise AttributeError("cannot change calculated information")
    @info.deleter
    def info(self):         del self._info
    #info = property(get_info, set_info, del_info, "I'm the 'dict' property, containing calculated properties.")

    @property
    def compartments_map(self):
        if not vars(self).has_key('_compartments_map') or self._compartments_map == {}:
            self._create_tree()
        return self._compartments_map
    @compartments_map.setter
    def compartments_map(self, value):  raise AttributeError("cannot change calculated information")
    @compartments_map.deleter
    def compartments_map(self):         self._compartments_map  = {}#del self._compartments_map
    #compartments_map    = property(get_compartments_map, set_compartments_map, del_compartments_map, "I'm the 'dict' property, containing compartments_map.")

    @property
    def root(self):
        if not vars(self).has_key('_root') or self._root == None:
            self._create_tree()
        return self._root
    @root.setter
    def root(self, value):        raise AttributeError("cannot change calculated information")
    @root.deleter
    def root(self):        self._root        = None #del self._root
    #root = property(get_root, set_root, del_root, "I'm the 'dict' property, containing the root property.")

    @property
    def biggest(self):
        if not vars(self).has_key('_biggest') or self._biggest == None:
            self._create_tree()
        return self._biggest
    @biggest.setter
    def biggest(self, value):        raise AttributeError("cannot change calculated information")
    @biggest.deleter
    def biggest(self):               self._biggest   = None #del self._biggest
    #biggest = property(get_biggest, set_biggest, del_biggest, "I'm the 'dict' property, containing the root property.")

    @property
    def leafs(self):
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
        return self._leafs
    @leafs.setter
    def leafs(self, value):        raise AttributeError("cannot change calculated information")
    @leafs.deleter
    def leafs(self):               self._leafs               = [] #del self._leafs
    #leafs = property(get_leafs, set_leafs, del_leafs, "I'm the 'dict' property, containing the root property.")

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
        return self._branching_points
    @branching_points.setter
    def branching_points(self, value):      raise AttributeError("cannot change calculated information")
    @branching_points.deleter
    def branching_points(self):             self._branching_points               = [] #del self._branching_point
    #branching_points = property(get_branching_points, set_branching_points, del_branching_points, "I'm the 'dict' property, containing the root property.")

    def getCompartments(self):
        '''
        generator over compartments
        '''
        for compartment in self.compartments:
            yield(compartment)

    def getCompartment(self, compartment_id):
        '''
        in order to be conform with Model definition, it is possible to get an compartment by i.
        '''
        return self.compartments_map[ compartment_id ];

    def getSubtree(self, compartment):
        if not vars(self).has_key('_compartments_map') or self._compartments_map == {}:
            self._create_tree()
        if type(compartment) == type(1):
            compartment = self.getCompartment(compartment)

        parts   = []
        todo_stack    = []
        parts.append( compartment )
        todo_stack.append( compartment )
        while len( todo_stack ) > 0:
            c   = todo_stack.pop()
            todo_stack.extend( c.children )
            parts.extend( c.children )
        return parts

    def _create_tree(self):
        if vars(self).has_key('_compartments_map') and self._compartments_map != {}:
            return
        self._compartments_map  = {}
        self._biggest   = self.compartments[0] if self.compartments[0].compartment_parent_id != -1 else self.compartments[1]

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

    def add_compartment(self, compartment):
        """
        doc-string to be added.
        """
        self.compartments.append(compartment)
        #compartment = Compartment(parent, id, length)
        #parent.children[id] = compartment

    def plot_all_properties(self, morphologies=[], picture_file=None, picture_formats=['png', 'pdf', 'svg']):
        self.plot(morphologies, quantity="cylindric_volume",                                                        picture_file=picture_file, picture_formats=picture_formats)
        self.plot(morphologies, quantity="frustum_volume",       yaxis_description=u'volume [µm³]',                 picture_file=picture_file, picture_formats=picture_formats)
        self.plot(morphologies, quantity="cylindric_surface_area",                                                  picture_file=picture_file, picture_formats=picture_formats)
        self.plot(morphologies, quantity="frustum_surface_area", yaxis_description=u'surface area [µm²]',           picture_file=picture_file, picture_formats=picture_formats)
        self.plot(morphologies, quantity="branches",             yaxis_description='#branches',                     picture_file=picture_file, picture_formats=picture_formats)
        self.plot(morphologies, quantity="path_length",          yaxis_description=u'total cell length [µm]',       picture_file=picture_file, picture_formats=picture_formats)
        self.plot(morphologies, quantity="cylindric_mcse",                                                          picture_file=picture_file, picture_formats=picture_formats)
        self.plot(morphologies, quantity="frustum_mcse",         yaxis_description=u"mean cross-section area [µm]", picture_file=picture_file, picture_formats=picture_formats)
        #TODO: self.plot(morphologies, quantity="spatial_strech", yaxis_description='spatial_strech [$\mu$m]',               picture_file=picture_file, picture_formats=picture_formats)

    def plot(self, morphologies=[], quantity='volumes', yaxis_description=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):
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

    def plot_endpoints_histogramm(self, xlim=None, ylim=None, picture_file=None, picture_formats=['png', 'pdf', 'svg']):  
        import matplotlib.pyplot
        import numpy
        #matplotlib.rc('text', usetex=True): error with names

        x   = []
        for c in self.leafs:
            #print "%i/%i\r" % (len(x),len(self.leafs)), 
            x.append(c.path_distance(self.biggest))
        mean   = numpy.mean(x)
        std    = numpy.std(x)

        matplotlib.pyplot.title('Endpoints of %s' % (self.name.replace('_',' ')) )
        print 'Endpoints of %s : mean=%f, std=%f' % (self.name, mean, std)
        
        matplotlib.pyplot.axvline(x=mean)
        matplotlib.pyplot.grid(True, color='lightgrey')
        matplotlib.pyplot.hist(x, 20, normed=0, color='black')#, label='my data'

        matplotlib.pyplot.ylabel('#')#%
        if ylim != None:
            matplotlib.pyplot.ylim(ylim)
        matplotlib.pyplot.xlabel(u'distance [µm]')
        if xlim != None:
            matplotlib.pyplot.xlim(xlim)

        xmin, xmax, ymin, ymax  = matplotlib.pyplot.axis()
        width	= std / (xmax-xmin)
        center	= (mean - xmin) / (xmax-xmin)
        matplotlib.pyplot.axhline(y=.5*(ymax-ymin), xmin=center-.5*width, xmax=center+.5*width, color='red')#TODO: höhe der Line = mittelwert der bins
        matplotlib.pyplot.legend( [ 'mean %f' % ( mean ), 'std %f' % ( std )  ] )

        if(picture_file != None):
            for picture_format in picture_formats:
                matplotlib.pyplot.savefig(picture_file+'.'+picture_format,format=picture_format)
        else:
            matplotlib.pyplot.show()
        matplotlib.pyplot.close()

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

    def getCompartment(self, compartment_id):
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





class Morphology_groups(morphjongleur.util.auto_string.Auto_string):
    pass

class Compartment_groups(morphjongleur.util.auto_string.Auto_string):
    pass

class Morphology_info(morphjongleur.util.auto_string.Auto_string):
    pass
class Morphology_info2(morphjongleur.util.auto_string.Auto_string):
    """
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
    
    
    #TODO: name
    #mean abstand
    #pca
    #konvexen polyeder zur not emal in mathesoftware
    
    


class Compartment_info(morphjongleur.util.auto_string.Auto_string):
    pass
class Compartment_info2(morphjongleur.util.auto_string.Auto_string):
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
