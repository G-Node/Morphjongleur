#http://michelanders.blogspot.de/2012/02/3d-convex-hull-in-python.html
#Z=2
#Y=1
#X=0

debug = False

class Vector:
	def __init__(self,x,y,z):
		self.x=x
		self.y=y
		self.z=z
	
	def __str__(self):
		return str(self.x)+" "+str(self.y)+" "+str(self.z)

	def __sub__(self,other):
		return Vector(other.x-self.x,other.y-self.y,other.z-self.z)
	
	def __add__(self,other):
		return Vector(other.x+self.x,other.y+self.y,other.z+self.z)

	@staticmethod
	def fromArray(na):
		assert len(na) == 3
		return Vector(na[0],na[1],na[2])
	
	def toNumpyArray(self):
		import numpy
		return numpy.array([self.x,self.y,self.z])

class Vertex:
	def __init__(self,v,vnum=None,duplicate=None,onhull=False,mark=False):
		self.v = v
		self.vnum = vnum
		self.duplicate = duplicate # ref to incident cone edge (or None)
		self.onhull = onhull # T iff point on hull.
		self.mark = mark # T iff point already processed.
	
	def __str__(self):
		return str(self.v)
	
	def debug(self):
		return "Vertex: %d %s dup:%s onhull:%s mark:%s\n"%(self.vnum,self.v,self.duplicate,self.onhull,self.mark)
		
	@staticmethod
	def Collinear(a,b,c):
		"""
		Collinear checks to see if the three vertices given are collinear,
		by checking to see if each element of the cross product is zero.
		"""
		return (( c.v.z - a.v.z ) * ( b.v.y - a.v.y ) -
         ( b.v.z - a.v.z ) * ( c.v.y - a.v.y ) == 0
      and ( b.v.z - a.v.z ) * ( c.v.x - a.v.x ) -
         ( b.v.x - a.v.x ) * ( c.v.z - a.v.z ) == 0
      and ( b.v.x - a.v.x ) * ( c.v.y - a.v.y ) -
         ( b.v.y - a.v.y ) * ( c.v.x - a.v.x ) == 0 )

class Edge:
	enum = 0
	
	def __init__(self,adjface=[None,None],endpts=[None,None],newface=None,delete=False):
		self.adjface = []
		self.adjface.extend(adjface)
		self.endpts = []
		self.endpts.extend(endpts)
		self.newface = newface # ref to incident cone face
		self.delete = delete # T iff edge should be deleted
		self.enum = Edge.enum
		Edge.enum+=1
		
	def __str__(self):
		af=[str(f.fnum) if not (f is None) else '.' for f in self.adjface]
		return "edge(%d): %s %s del: %s newf: %d adjf:%s"%(self.enum,self.endpts[0],self.endpts[1],self.delete,self.newface.fnum if not (self.newface is None) else -1," ".join(af))
		
class Face:
	fnum=0
	def __init__(self,edge=[None,None,None],vertex=[None,None,None],visible=False):
		self.edge = []
		self.edge.extend(edge)
		self.vertex = []
		self.vertex.extend(vertex)
		self.visible = visible # T iff face visible from new point
		self.fnum = Face.fnum
		Face.fnum+=1
		
	def __str__(self):
		return """facet normal 0 0 0
outer loop
vertex %s
vertex %s
vertex %s
endloop
endfacet"""%(self.vertex[0],self.vertex[1],self.vertex[2])

	def debug(self):
		return "Face(%d) visble: %s\n\t%s\t%s\t%s"%(self.fnum,self.visible,self.vertex[0].debug(),self.vertex[1].debug(),self.vertex[2].debug())
		
	def InitEdges(self, fold=None):
		v0=self.vertex[0]
		v1=self.vertex[1]
		v2=self.vertex[2]
		newedges=[]
		# Create edges of the initial triangle
		if fold is None:
			e0 = Edge()
			e1 = Edge()
			e2 = Edge()
			newedges=[e0,e1,e2]
		else: # Copy from fold, in reverse order
			e0 = fold.edge[2]
			e1 = fold.edge[1]
			e2 = fold.edge[0]
		e0.endpts[0] = v0
		e0.endpts[1] = v1
		e1.endpts[0] = v1
		e1.endpts[1] = v2
		e2.endpts[0] = v2
		e2.endpts[1] = v0
	
		self.edge[0] = e0
		self.edge[1] = e1
		self.edge[2] = e2
	
		# Link edges to face
		e0.adjface[0] = self
		e1.adjface[0] = self
		e2.adjface[0] = self
		
		return newedges
		
	def MakeCcw(f,e,p): # the customary self is called f here instead of self
		"""
		MakeCcw puts the vertices in the face structure in counterclock wise 
		order.  We want to store the vertices in the same 
		order as in the visible face.  The third vertex is always p.

		Although no specific ordering of the edges of a face are used
		by the code, the following condition is maintained for each face f:
		one of the two endpoints of f.edge[i] matches f.vertex[i]. 
		But note that this does not imply that f.edge[i] is between
		f.vertex[i] and f.vertex[(i+1)%3].  (Thanks to Bob Williamson.)
		"""
		
		# fv	The visible face adjacent to e
		# i     Index of e.endpoint[0] in fv

		if e.adjface[0].visible:      
			fv = e.adjface[0]
		else:
			fv = e.adjface[1]

		# Set vertex[0] & [1] of f to have the same orientation
		# as do the corresponding vertices of fv
		i=0
		while(fv.vertex[i] != e.endpts[0]):
			i+=1
		# Orient f the same as fv
		if fv.vertex[ (i+1) % 3 ] != e.endpts[1] :
			f.vertex[0] = e.endpts[1]  
			f.vertex[1] = e.endpts[0]    
		else:
			f.vertex[0] = e.endpts[0]   
			f.vertex[1] = e.endpts[1]   
			(f.edge[1], f.edge[2] ) = (f.edge[2], f.edge[1] )
			# This swap is tricky. e is edge[0]. edge[1] is based on endpt[0],
			# edge[2] on endpt[1].  So if e is oriented "forwards," we
			# need to move edge[1] to follow [0], because it precedes. */

		f.vertex[2] = p
		
# Define flags
ONHULL = True
REMOVED = True
VISIBLE  = True
PROCESSED = True

class Hull:
	def __init__(self,v):
		self.vertices = []
		self.edges = []
		self.faces = []
		self.ReadVertices(v)
		v=self.DoubleTriangle()
		self.ConstructHull(v)
		self.EdgeOrderOnFaces()

	def ReadVertices(self,v):
		self.vertices = [ Vertex(vc,i) for i,vc in enumerate(v) ]
		
	def EdgeOrderOnFaces(self):
		"""
		EdgeOrderOnFaces: puts e0 between v0 and v1, e1 between v1 and v2,
		e2 between v2 and v0 on each face.  This should be unnecessary, alas.
		"""
		for f in self.faces:
			for i in (0,1,2):
				if ( not (((f.edge[i].endpts[0] == f.vertex[i]) and
					(f.edge[i].endpts[1] == f.vertex[(i+1)%3])) or
					((f.edge[i].endpts[1] == f.vertex[i]) and
					(f.edge[i].endpts[0] == f.vertex[(i+1)%3])))):
					# Change the order of the edges on the face
					for j in (0,1,2):
						# find the edge that should be there
						if (((f.edge[j].endpts[0] == f.vertex[i]) and
							(f.edge[j].endpts[1] == f.vertex[(i+1)%3])) or
							((f.edge[j].endpts[1] == f.vertex[i]) and
							(f.edge[j].endpts[0] == f.vertex[(i+1)%3]))) :
							# Swap it with the one erroneously put into its place
							(f.edge[i],f.edge[j]) = (f.edge[j],f.edge[i])

	def __str__(self):
		str_list	= ['solid points'] 
		str_list.extend( [str(f) for f in self.faces] )
		str_list.append('endsolid points')
		return '\n'.join(str_list)
		
	def write(self, f):
		'''
		http://www.subdude-site.com/WebPages_Local/RefInfo/Computer/Linux/LinuxGuidesByBlaze/apps3Dtools/3D_viewers-converters/3DviewersANDconverters_intro.htm
		'''
		stream  = open(f+'.stl', 'w')
		stream.write(str(self))
		stream.close()

	def debug(self,msg=''):
		s=[msg+'\n']
		for f in self.faces:
			s.append(f.debug())
		s.append('-'*40)
		return "".join(s)
		
	@staticmethod
	def VolumeSign(f,p):
		"""
		VolumeSign returns the sign of the volume of the tetrahedron determined by f
		and p.  VolumeSign is +1 iff p is on the negative side of f,
		where the positive side is determined by the rh-rule.  So the volume 
		is positive if the ccw normal to f points outside the tetrahedron.
		The final fewer-multiplications form is due to Bob Williamson.
		
		This implementation differs from the one in the book in that it does not assume that
		coordinates are integers.
		"""
		a=f.vertex[0].v - p.v
		b=f.vertex[1].v - p.v
		c=f.vertex[2].v - p.v
		
		vol = ( a.x * (b.y*c.z - b.z*c.y)
         + a.y * (b.z*c.x - b.x*c.z)
         + a.z * (b.x*c.y - b.y*c.x) )

		# If the volume should be an integer, make epsilon 0.5
		epsilon = 1e-10
		if vol >  epsilon: return  1
		if vol < -epsilon: return -1
		return 0

	def DoubleTriangle(self):
		"""
		 DoubleTriangle builds the initial double triangle.  It first finds 3 
		 noncollinear points and makes two faces out of them, in opposite order.
		 It then finds a fourth point that is not coplanar with that face.  The  
		 vertices are stored in the face structure in counterclockwise order so 
		 that the volume between the face and the point is negative. Lastly, the
		 3 newfaces to the fourth point are constructed and the data structures
		 are cleaned up. 
		"""

		# Find 3 noncollinear points
		v0 = 0
		nv = len(self.vertices)
		while(Vertex.Collinear(self.vertices[v0%nv],self.vertices[(v0+1)%nv],self.vertices[(v0+2)%nv])):
			v0 = (v0+1)%nv
			if v0 == 0:
				raise Exception("DoubleTriangle:  All points are Collinear!")
				
		v1 = (v0+1)%nv
		v2 = (v1+1)%nv
	
		# Mark the vertices as processed
		self.vertices[v0].mark = PROCESSED
		self.vertices[v1].mark = PROCESSED
		self.vertices[v2].mark = PROCESSED

		# Create the two "twin" faces
		self.faces.append(Face(vertex=[self.vertices[v0],self.vertices[v1],self.vertices[v2]]))
		f0=self.faces[-1]
		self.edges.extend(f0.InitEdges())
		self.faces.append(Face(vertex=[self.vertices[v2],self.vertices[v1],self.vertices[v0]]))
		f1=self.faces[-1]
		self.edges.extend(f1.InitEdges(f0))

		# Link adjacent face fields.
		f0.edge[0].adjface[1] = f1
		f0.edge[1].adjface[1] = f1
		f0.edge[2].adjface[1] = f1
		f1.edge[0].adjface[1] = f0
		f1.edge[1].adjface[1] = f0
		f1.edge[2].adjface[1] = f0
	
		#Find a fourth, noncoplanar point to form tetrahedron
		v3 = (v2+1)%nv
		vol = self.VolumeSign( f0, self.vertices[v3] )
		while vol == 0:
			v3 = (v3+1)%nv
			if v3==0:
				raise Exception("DoubleTriangle:  All points are coplanar!")
			vol = self.VolumeSign( f0, self.vertices[v3] )
	
		if debug: print(self.debug('initial'))
		
		return v3
		
	def ConstructHull(self,v):
		"""
		ConstructHull adds the vertices to the hull one at a time.  The hull
		vertices are those in the list marked as onhull.
		"""

		# vertices is supposed to be a circular list that we traverse once, starting at v
		# however, the call to CleanUp may delete vertices from this list
		ev = v
		while(True):
			if not self.vertices[v].mark:
				self.vertices[v].mark = PROCESSED;
			self.AddOne(self.vertices[v]);
			ev,v=self.CleanUp(ev,v) # cleanup may delete vertices!
			if v == ev : break
			
	def AddOne(self,p):
		"""
		AddOne is passed a vertex.  It first determines all faces visible from 
		that point.  If none are visible then the point is marked as not 
		onhull.  Next is a loop over edges.  If both faces adjacent to an edge
		are visible, then the edge is marked for deletion.  If just one of the
		adjacent faces is visible then a new face is constructed.
		"""

		vis = False;

		# Mark faces visible from p.
		for f in self.faces:
			vol = self.VolumeSign( f, p )
			if vol < 0 : f.visible = VISIBLE
			vis = True;                      

		# If no faces are visible from p, then p is inside the hull
		if not vis:
			p.onhull = not ONHULL  
			return False 

		# Mark edges in interior of visible region for deletion.
		# Erect a newface based on each border edge
		for e in self.edges:
			if e.adjface[0].visible and	e.adjface[1].visible:
				# e interior: mark for deletion
				e.delete = REMOVED;
			elif e.adjface[0].visible or e.adjface[1].visible: 
				# e border: make a new face
				e.newface = self.MakeConeFace( e, p )
		
		if debug : print(self.debug('addone'))
		
		return True

	def MakeConeFace(self,e,p):
		"""
		MakeConeFace makes a new face and two new edges between the 
		edge and the point that are passed to it. It returns a pointer to
		the new face.
		"""
		
		new_edge=[None,None]
		# Make two new edges (if don't already exist)
		for i in (0,1):
			# If the edge exists, copy it into new_edge
			# Otherwise (duplicate is NULL), MakeNullEdge
			d = e.endpts[i].duplicate
			if d is None:
				new_edge[i] = Edge(endpts=[e.endpts[i],p])
				e.endpts[i].duplicate = new_edge[i]
				self.edges.append(new_edge[i])
			else:
				new_edge[i] = d

		# Make the new face
		new_face = Face(edge=[e,new_edge[0],new_edge[1]])
		self.faces.append(new_face)
		new_face.MakeCcw( e, p )

		# Set the adjacent face pointers
		for i in (0,1):
			for j in (0,1):
				# Only one None link should be set to new_face
				if new_edge[i].adjface[j] is None:
					new_edge[i].adjface[j] = new_face
					break
		return new_face

	def CleanUp(self,ev,v):
		"""
		CleanUp goes through each data structure list and clears all
		flags and NULLs out some pointers.  The order of processing
		(edges, faces, vertices) is important.
		"""
		de=self.CleanEdges()
		if debug: print(self.debug('cleanedges '+" ".join(de)))
		self.CleanFaces()
		if debug: print(self.debug('cleanfaces'))
		ev,v=self.CleanVertices(ev,v)
		if debug: print(self.debug('cleanvertices'))
		return ev,v
		
	def CleanEdges(self):
		"""
		CleanEdges runs through the edge list and cleans up the structure.
		If there is a newface then it will put that face in place of the 
		visible face and NULL out newface. It also deletes so marked edges.
		"""
		# Integrate the newface's into the data structure
		for e in self.edges:
			if e.newface:
				if e.adjface[0].visible:
					e.adjface[0] = e.newface 
				else:
					e.adjface[1] = e.newface
			e.newface = None

		# Delete any edges marked for deletion. */
		deleted_edges = [str(e.enum) for e in self.edges if e.delete ]
		self.edges = [e for e in self.edges if not e.delete ]
		return deleted_edges
				
	def CleanFaces(self):
		"""
		CleanFaces runs through the face list and deletes any face marked visible.
		"""

		self.faces = [f for f in self.faces if not f.visible ]
		
	def CleanVertices(self,evi,vi):
		"""
		CleanVertices runs through the vertex list and deletes the 
		vertices that are marked as processed but are not incident to any 
		undeleted edges. 
		"""
		# Mark all vertices incident to some undeleted edge as on the hull
		for e in self.edges:
			e.endpts[0].onhull = ONHULL
			e.endpts[1].onhull = ONHULL
		
		# Delete all vertices that have been processed but are not on the hull
		for i,v in enumerate(self.vertices):
			if v.mark and not v.onhull:
				del self.vertices[i]
				if i<evi : evi -= 1
				vi -= 1
		vi = (vi+1)%len(self.vertices)

		# Reset flags
		for v in self.vertices:
			v.duplicate = None
			v.onhull = not ONHULL
			
		return evi,vi

	def surface_area_and_volume(self):
		'''
		return (surface_area, v)
		
		http://stackoverflow.com/questions/451426/how-do-i-calculate-the-surface-area-of-a-2d-polygon
		http://www.mathopenref.com/coordpolygonarea.html
		http://www.python-kurs.eu/matrix_arithmetik.php
		http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
		
		1/3 * G + h
		http://www.onlinemathe.de/forum/Volumen-und-Grundflaechenberechnung-Pyramide
		ab =numpy.array([8,14,-8]); ac =numpy.array([12,3,-3]); ad = numpy.array([5,5,7]);s = numpy.dot(np.cross(ab,ac), ad); abs(s)/6.
		'''
		import numpy
		assert len(self.faces) > 3#TODO: execpt
		surface_area	= 0
		volume	= 0
		d	= self.faces[0].vertex[0].v.toNumpyArray()
		for face in self.faces:
			assert len(face.vertex) == 3
			a	= face.vertex[0].v.toNumpyArray()
			b	= face.vertex[1].v.toNumpyArray()
			c	= face.vertex[2].v.toNumpyArray()
			ab	= b - a 
			ac	= c - a
			ad	= d - a
			
			#cos	= ab.dot(ac) / (numpy.sqrt(ab.dot(ab)) * numpy.sqrt(ac.dot(ac)))
			#sin	= numpy.sqrt(1-cos*cos)
			#sa		= .5 * numpy.sqrt(ab.dot(ab)) * numpy.sqrt(ac.dot(ac)) * sin
			
			k	= numpy.cross(ab,ac)
			sa	= numpy.sqrt(numpy.dot(k,k))/2.	# norm
			#print sa
			surface_area	+= sa
			
			s	= numpy.dot(k, ad)	#spatprodukt
			v	= abs(s)/6
			volume	+= v

		self._surface_area	= surface_area
		self._volume		= volume
		return (surface_area, volume)

	def surface_area(self):
		if not vars(self).has_key('_surface_area') or self._surface_area == None:
			self.surface_area_and_volume()
		return self._surface_area
	
	def volume(self):
		if not vars(self).has_key('_volume') or self._volume == None:
			self.surface_area_and_volume()
		return self._volume

if __name__ == "__main__":
	from random import random
	# simple cube, integer coordinates
	cube=[Vector(0,0,0),Vector(1,0,0),Vector(0,1,0),Vector(1,1,0),Vector(0,0,1),Vector(1,0,1),Vector(0,1,1),Vector(1,1,1)]
	# irregular tetrahedron
	tetrahedron=[Vector(0,0,0),Vector(1,0,0),Vector(0,1,0),Vector(1,1,1)]
	# hexahedron
	hexahedron=[Vector(0,0,0),Vector(1,0,0),Vector(0,1,0),Vector(1,1,1),Vector(1,1,-1)]
	# cube with an internalpoint, integer coordinates
	cube_internal=[Vector(0,0,0),Vector(2,0,0),Vector(0,2,0),Vector(2,2,0),Vector(0,0,2),Vector(2,0,2),Vector(0,2,2),Vector(2,2,2),Vector(1,1,1)]
	# pyramid, fractional coordinates
	pyramid=[Vector(0,0,0),Vector(0,1,0),Vector(1,0,0),Vector(1,1,0),Vector(0.5,0.5,0.3)]
	# simple cube with many internal points with random fractional coordinates
	cubef=[Vector(0,0,0),Vector(1,0,0),Vector(0,1,0),Vector(1,1,0),Vector(0,0,1),Vector(1,0,1),Vector(0,1,1),Vector(1,1,1)]
	pyramide=[Vector(0,0,0),Vector(8,14,-8),Vector(12,3,-3),Vector(5,5,7)]#http://www.onlinemathe.de/forum/Volumen-und-Grundflaechenberechnung-Pyramide
	
	for i in range(20):
		cubef.append(Vector(random(),random(),random()))
	
	sphere=[]
	for i in range(2000):
		x,y,z = 2*random()-1,2*random()-1,2*random()-1
		if x*x+y*y+z*z < 1.0:
			sphere.append(Vector(x,y,z))

	#h=Hull(sphere)
	#h=Hull(cube)
	#h=Hull(tetrahedron)
	#h=Hull(hexahedron)
	#h=Hull(cube_internal)
	#h=Hull(pyramid)
	h=Hull(pyramide)
	#h=Hull(sphere)
	print(h.debug("#%i faces with %f surface area %f volume" % (len(h.faces), h.surface_area(), h.volume()) ))
	h.write('/tmp/chull')