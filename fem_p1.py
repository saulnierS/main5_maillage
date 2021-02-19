#****************************************
# FEM P1: solveur
#****************************************
# Robin Clément & Saulnier Solène
# MAIN5  02/2021
#****************************************
#-----------
# Packages
#-----------
# system
import sys
# maths
import numpy as np
# maillage
import gmsh
# files  
import common
import maillage


############
# Functions
############

##############
# mass matrix
##############
def loc2glob(element, ind:int):
	"""Convert local coordonate of the mesh in global coordonate
	element: triangle or segment of the mesh
	ind: local ind

    Returns global ind: int."""
	return element.points[ind].id


def mass_ref(i:int,j:int):
	"""reference matrix mass
	i=int
    j=int
    (i,j)=coordonate
	###################
	# 1/12  1/24  1/24
	# 1/24  1/12  1/24
	# 1/24  1/24  1/12
	###################
    Returns global ind: int."""
	if i==j:
		return 1/12
	return 1/24


def mass_elem(element, triplets, alpha =1.):
	"""compute elementary mass matrix
	element : Triangle
    triplets :Triplets
    alpha: scalar optionnal
    Returns 0
    Warning, triplets which is pass in argument is modified during the function"""
	for i in range(0,3):
		I = loc2glob(element,i)

		for j in range(0,3):
			J = loc2glob(element,j)

			#area = 1/2*det(jac) so 2*area=det(jac)
			val = mass_ref(i,j)*element.area()*2*alpha

			triplets.append(I,J,val)

	return 0


def Mass(msh, dim:int, physical_tag:int, triplets):
	""" compute mass matrix
    msh : Mesh
    dim : int
    physical_tag : int
    triplets : Triplets
    Returns 0
    Warning, triplets which is pass in argument is modified during the function"""
	triangles=msh.getElements(dim,physical_tag)
	for i in range(0,len(triangles)):
		mass_elem(triangles[i], triplets)

	return 0

#########################
# stiffness matrix
#########################
 
def B(element):
	""" Compute the matrix B = (jacT)-1
    element : Triangle or segment
    Returns matrix B : array"""
	J = element.jac()
	a = J[0][0]
	b = J[0][1]
	c = J[1][0]
	d = J[1][1]
	return 1/(a*d-c*b)*np.array([[d,-c],[-b,a]])


def gradPhi(element, i:int):
	""" compute phi gradient
	element: Triangle
	i : int chose the right phi wanted
    phi_0 = 1-eta-xsi
    phi_1 = xsi
    phi_2 = eta
    Returns the gradient of phi function : array """
	if i==0:
		return np.array([[-1],[-1]])
	if i==1:
		return np.array([[1],[0]])
	if i==2:
		return np.array([[0],[1]])
	return -1


def stiffness_elem(element, triplets):
	"""compute elementary stiffness matrix
	element : Triangle
    triplets :Triplets
    Returns 0
    Warning, triplets which is pass in argument is modified during the function"""

	for i in range(0,3):
		I = loc2glob(element,i)

		for j in range(0,3):
			J = loc2glob(element,j)

			#area = 1/2*det(jac) so 2*area=det(jac) 
			res_1=np.matmul(np.transpose(gradPhi(element,j)),np.transpose(element.B()))
			res_2 = np.matmul(element.B(), gradPhi(element,i))
			res_final = np.matmul(res_1,res_2)
			# De= det(jac)*(grad phi j)^T*B^T*B*(grad phi i)*aire du triangle
			# val= 2*element.area()*res_final[0][0]*element.area()
			val = element.area()*res_final[0][0]
			# print(val[0][0])
			triplets.append(I,J,val)

	return 0


def Stiffness(msh, dim:int, physical_tag:int, triplets):
	""" compute stiffness matrix
    msh : Mesh
    dim : int
    physical_tag : int
    triplets : Triplets
    Returns 0
    Warning, triplets which is pass in argument is modified during the function"""
	triangles = msh.getElements(dim,physical_tag)
	for i in range(0,len(triangles)):
		stiffness_elem(triangles[i], triplets)

	return 0

################################
# Quadratures
################################
def phiRef(i:int, param):
	""" compute the value of phi ref depends on param
    i : int (the phi wanted)
    param : array (coordonate (xsi,eta) or (x,y))
    phi_0 = 1-eta-xsi
    phi_1 = xsi
    phi_2 = eta
    Returns the value of the phi ref : float."""
	if i==0:
		return 1-param[0]-param[1]
	if i==1:
		return param[0]
	if i==2:
		return param[1] 
	return 0

def interpol_geo(element,m):
	""" compute the point xm used in the quadratures formula
    element : Segment
    m : int
    Returns the value of xm : array."""
	res = [0,0]
	_, pts_param, pts_phys=element.gaussPoint()
	for i in range(0,len(element.points)):
		# x du point de gauss passage de xsi, eta à x, y
		res[0]=res[0] + phiRef(i, pts_param[m])*pts_phys[i][0]
		res[1]=res[1] + phiRef(i, pts_param[m])*pts_phys[i][1]
	return res



def Integrale(msh, dim:int, physical_tag:int, f, B, order=2):
	""" compute the integrale by interpolation
    msh : Mesh
    dim : int
    physical_tag : int
    f : function 
    B : array (second member Ax=b)
    order : int (order of the method)
    Returns 0
    Warning, triplets and B which is pass in argument is modified during the function"""
	elements = msh.getElements(dim,physical_tag)
	# print(B[0])
	for ind_elem in range(0,len(elements)):
		# print(elements[ind_elem])
		w, pts_param, pts_phys = elements[ind_elem].gaussPoint()
		for i in range(0,3):
			I = loc2glob(elements[ind_elem],i)
			for m in range(0,len(w)):
				x_interpol=interpol_geo(elements[ind_elem],m)
				# somme ( somme ( somme (pds_m * f(x(xsi_m,eta_m)* phi(xsi_m,eta_m)))))
				B[I]=B[I]+w[m]*elements[ind_elem].area()*2*f(x_interpol[0],x_interpol[1])*phiRef(i, pts_param[m])

	return 0

###################
# Dirichlet
###################
# def Dirichlet(msh, dim:int , physical_tag:int , g, triplets, B):

# 	segments=msh.getElements(dim,physical_tag)

# 	for i in range(0,len(segments)):
# 		id_seg=segments[i].id
# 		print("segment")
# 		print(segments[i])
# 		for i_s in range(0,len(segments[i].points)):
# 			sommet = segments[i].points[i_s].id
# 			print("sommet")
# 			print(sommet)
# 			for j in range(0,triplets.size_data):
# 				#all val = triplets.data[0]
# 				#all I = triplets.data[1][0]
# 				#all J = triplets.data[1][1]
# 				I=triplets.data[1][0][j]
# 				J=triplets.data[1][1][j]
# 				if sommet==I:
# 					triplets.data[0][j]=0

# 			B[sommet] = g(segments[i].points[i_s].x,segments[i].points[i_s].y)
# 			# B[id_seg] = 0 
# 			triplets.append(sommet,sommet,1)

# 	return 0

def Dirichlet(msh, dim:int , physical_tag:int , g, triplets, B):
	""" add the dirichlet condition
    msh : Mesh
    dim : int
    physical_tag : int
    g : int (dirichlet condition value)
    triplets : Triplets
    B : array (second member Ax=b)
    Returns the new matrix with dirichlet condition: triplets.
    Warning, triplets and B which is pass in argument is modified during the function"""
	pts=msh.getPoints(dim,physical_tag)

	for ind_s in range(0,len(pts)):
		sommet = pts[ind_s].id

		for j in range(0,triplets.size_data):
			#all val = triplets.data[0]
			#all I = triplets.data[1][0]
			#all J = triplets.data[1][1]
			I=triplets.data[1][0][j]
			# J=triplets.data[1][1][j]
			if sommet==I:
				triplets.data[0][j]=0
				B[sommet] = g(pts[ind_s].x,pts[ind_s].y)

		triplets.append(sommet,sommet,1)

	return 0

