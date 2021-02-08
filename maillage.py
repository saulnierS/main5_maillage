
import gmsh
import sys
import math
import omega
import numpy as np

gmsh.initialize(sys.argv)
model = gmsh.model
model.add("omega")

##############################################
# recuperation gmsh info dans une classe mesh
##############################################
import numpy as np

class Mesh:
	"""class mesh"""
	def __str__(self):

		chaine="***Points***\n"
		for i in range(0,len(self.points)):
			chaine = chaine + "[id_pts: "+str(self.points[i].id)+"] " 
			chaine = chaine +"(x: "+str(self.points[i].x)+", y: "+str(self.points[i].y)+")\n"
		chaine = chaine + "\n************************************\n"


		chaine = chaine + "\n***Segments***\n"
		for i in range(0, len(self.segments)):
			chaine = chaine +"[id_sgs: "+str(self.segments[i].id)+"]"
			for j in range(0,len(self.segments[i].points)):
				chaine = chaine +"\n \t [id_pts: "+str(self.segments[i].points[j].id)+"] "
				chaine = chaine +" (x: "+str(self.segments[i].points[j].x)+", y: "+str(self.segments[i].points[j].y)+")"
			chaine = chaine +"\n"
		chaine = chaine +"\n************************************\n"

		chaine = chaine +"\n***Triangles***\n"
		for i in range(0, len(self.triangles)):
			chaine = chaine +"[id_tgs: "+str(self.triangles[i].id)+"]"
			for j in range(0,len(self.triangles[i].points)):
				chaine = chaine + "\n \t [id_pts: "+str(self.triangles[i].points[j].id)+"] "
				chaine = chaine +" (x: "+str(self.triangles[i].points[j].x)+", y: "+str(self.triangles[i].points[j].y)+")"
			chaine = chaine +"\n"
		chaine=chaine+"\n************************************\n\n\n"

		return chaine

	def GmshToMesh(self,filename,gmsh):
		"""build the maillage from filename.msh"""

		# name check
		assert filename[-4:]==".msh"

		# initialisation
		self.points = []
		self.segments = []
		self.triangles = []

		#reading
		# gmsh.merge(filename)
		##################
		### save nodes ###
		##################
		# gmsh.model.mesh.getNodes(dim=-1, tag=-1, includeBoundary=False, returnParametricCoord=True)
  		# ---Return `nodeTags', `coord', `parametricCoord'.
		nodes_tag = gmsh.model.mesh.getNodes()[0]
		for i_node in range(0,len(nodes_tag)):
			# gmsh.model.mesh.getNode(nodeTag)
  			# ---Return `coord', `parametricCoord'.
			coord = gmsh.model.mesh.getNode(nodes_tag[i_node])[0]
			self.points.append(Point(coord[0],coord[1],int(nodes_tag[i_node]-1)))
			#NB: keep the tag for id of this code
		
		###################################
		### save segments and triangles ###
		###################################

		# gmsh.model.getPhysicalGroups(dim=-1)
        # ---Return `dim,Tags'.
		physical_group= gmsh.model.getPhysicalGroups()
		for i_group in range(0,len(physical_group)):

			# save dim
			pg_dim = physical_group[i_group][0]
			#save physical tag
			pg_tag = physical_group[i_group][1]

			# gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            # ---Return `tags'.
			entities = gmsh.model.getEntitiesForPhysicalGroup(pg_dim,pg_tag)
			for i_entity in range(0,len(entities)):

				# gmsh.model.mesh.getElements(dim=-1, tag=-1)
                # ---Return `elementTypes', `elementTags', `nodeTags'.
				elements = gmsh.model.mesh.getElements(pg_dim, entities[i_entity])[1][0]
				for e in range(0,len(elements)):

					# msh.model.mesh.getElement(elementTag)
                    # ---Return `elementType', `nodeTags'.
					element=gmsh.model.mesh.getElement(elements[e])[1]
					# save points for each element
					p_tmp = []
					for sommet in range(0,len(element)):
						coord = gmsh.model.mesh.getNode(element[sommet])[0]
						p_tmp.append(Point(coord[0],coord[1],int(element[sommet]-1)))

					if pg_dim==1:
						self.segments.append(Segment(p_tmp,pg_tag,int(elements[e]-1)))

					if pg_dim==2:
						self.triangles.append(Triangle(p_tmp,pg_tag,int(elements[e]-1)))

		################			
		#### display ###
		################
		# for p in range(0,len(self.points)):
		# 	print(self.points[p])
		# for s in range(0, len(self.segments)):
		# 	print(self.segments[s])
		# for t in range(0, len(self.triangles)):
		# 	print(self.triangles[t])
		self.Npts = len(self.points)
		self.Nsgs = len(self.segments)
		self.Ntgs = len(self.triangles)

		
	def getElements(self,dim,physical_tag):
		res=[]
		if dim == 1: 
			for i in range(0,len(self.segments)):
				if self.segments[i].physical_tag == physical_tag:
					res.append(self.segments[i])
		if dim == 2: 
			for i in range(0,len(self.triangles)):
				if self.triangles[i].physical_tag == physical_tag:
					res.append(self.triangles[i])
		return res

	def getPoints(self,dim,physical_tag):
		res=[]
		val=[]
		if dim == 1: 
			for i in range(0,len(self.segments)):
				if self.segments[i].physical_tag == physical_tag:
					for j in range(0,len(self.segments[i].points)):
						if self.segments[i].points[j].id not in val:
							res.append(self.segments[i].points[j])
							val.append(self.segments[i].points[j].id)

		if dim == 2: 
			for i in range(0,len(self.triangles)):
				if self.triangles[i].physical_tag == physical_tag:
					for j in range(0,len(self.triangles[i].points)):
						if self.triangles[i].points[j].id not in val:
							res.append(self.triangles[i].points[j])
							val.append(self.triangles[i].points[j].id)
		return res


		
		
class Point (Mesh):
	"""class point"""
	def __init__(self,x,y,id_point):
		Mesh.__init__(self)
		self.id = id_point
		self.x = x
		self.y = y


	def __str__(self):
		return "[id_pts: "+str(self.id)+"] (x: "+str(self.x)+", y: "+str(self.y)+")"

class Segment (Mesh):
	"""class segment"""
	def __init__(self,points,physical_tag,id_segment):
		Mesh.__init__(self)
		self.physical_tag = physical_tag
		self.id = id_segment
		self.points = points
		self.name = "Segment"


	def __str__(self):
		chaine="[id_sgts: "+str(self.id)+"]"
		for i in range(0,len(self.points)):
			chaine = chaine + "\n \t [id_pts: "+str(self.points[i].id)+"] (x: "+str(self.points[i].x)+", y: "+str(self.points[i].y)+")"
		return chaine

	def area(self):
		return 0

	def jac(self):
		return 0

class Triangle (Mesh):
	"""class triangle"""
	def __init__(self,points,physical_tag,id_trigangle):
		Mesh.__init__(self)
		self.physical_tag = physical_tag
		self.id = id_trigangle
		self.points = points
		self.name = "Triangle"


	def __str__(self):
		chaine="[id_tgs: "+str(self.id)+"]"
		for i in range(0,len(self.points)):
			chaine = chaine + "\n \t [id_pts: "+str(self.points[i].id)+"] (x: "+str(self.points[i].x)+", y: "+str(self.points[i].y)+")"
		return chaine

	def area(self):
		a = self.points[1].x-self.points[0].x
		b = self.points[2].x-self.points[0].x
		c = self.points[1].y-self.points[0].y
		d = self.points[2].y-self.points[0].y
		return np.abs(1/2*(a*d-c*b))

	def jac(self):
		#####################  
		# J= (x1-x0  x2-x0) #
		#    (y1-y0  y2-y0) #
		#####################

		a = self.points[1].x-self.points[0].x
		b = self.points[2].x-self.points[0].x
		c = self.points[1].y-self.points[0].y
		d = self.points[2].y-self.points[0].y
		return np.array([[a,b],[c,d]])
	#pour le calcul de la matrice de rigidite
	
	def B(self):
		#####################
		# B= JT^-1  
		# B= (y2-y0  y0-y1) #
		#    (x0-x2  x1-x0) #
		#####################

		a = self.points[1].x-self.points[0].x
		b = self.points[2].x-self.points[0].x
		c = self.points[1].y-self.points[0].y
		d = self.points[2].y-self.points[0].y
		return np.array([[d,-c],[-b,a]])
	#pour la quadrature
	def gaussPoint(self,order=2):
		poids = [1/6]

		if order==1:
			poids = [1/6]
			param = [(1/3,1/3)]

		if order==2:
			poids = [1/6,1/6,1/6]
			param = [(1/6,1/6),(4/6,1/6),(1/6,4/6)]

		pts=[]
		for i in range(0,len(self.points)):
			pts.append((self.points[i].x,self.points[i].y))
		return (poids, param, pts)


gmsh.finalize()