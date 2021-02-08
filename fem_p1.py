#################################
# fonction des elements finis P1
#################################

from scipy.sparse import coo_matrix
import maillage
import numpy as np
import common
import gmsh
import sys
from scipy import sparse
from scipy.sparse import linalg


def build_matrix(data):
	return coo_matrix(data).tocsr()

###################
# matrice de masse
###################
def loc2glob(element, ind:int):
	return element.points[ind].id


def mass_ref(i:int,j:int):
	###################
	# 1/12  1/24  1/24
	# 1/24  1/12  1/24
	# 1/24  1/24  1/12
	###################
	if i==j:
		return 1/12
	return 1/24

# element = Segment ou Triangle
# triplets = Triplets
# alpha un scalaire optionnel
def mass_elem(element, triplets, alpha =1.):

	for i in range(0,3):
		I = loc2glob(element,i)

		for j in range(0,3):
			J = loc2glob(element,j)

			#area = 1/2*det(jac) so 2*area=det(jac)
			val = mass_ref(i,j)*element.area()*2

			triplets.append(I,J,val)

	return triplets

# msh = Mesh
# dim = int
# physical_tag = int
# triplets = Triplets
def Mass(msh, dim:int, physical_tag:int, triplets):
	triangles=msh.getElements(dim,physical_tag)
	for i in range(0,len(triangles)):
		mass_elem(triangles[i], triplets)

	return 0

#########################
# matrice de rigidité
#########################

#il y a deja un B dans la classe mesh pour un triangle
def B(element):
	J = element.jac()
	a = J[0][0]
	b = J[0][1]
	c = J[1][0]
	d = J[1][1]
	return np.array([[d,-c],[-b,a]])

def gradPhi(element, i:int):
	# if i==0:
	# 	return np.array([-1,-1])
	# if i==1:
	# 	return np.array([1,0])
	# if i==2:
	# 	return np.array([0,1])
	if i==0:
		return np.array([[-1],[-1]])
	if i==1:
		return np.array([[1],[0]])
	if i==2:
		return np.array([[0],[1]])
	return -1



# element = Segment ou Triangle
# triplets = Triplets
# alpha un scalaire optionnel
def stiffness_elem(element, triplets):

	for i in range(0,3):
		I = loc2glob(element,i)

		for j in range(0,3):
			J = loc2glob(element,j)

			#area = 1/2*det(jac) so 2*area=det(jac)
			#D(i,j)=det(jac)* 
			res_dot_B=np.matmul(np.transpose(element.B()),element.B())
			res_dot_tmp = np.matmul(res_dot_B, gradPhi(element,i))
			res_final_dot = np.matmul(np.transpose(gradPhi(element,j)),res_dot_tmp)
			# res_dot_B=np.dot(np.transpose(element.B()),element.B())
			# res_dot_tmp = np.dot(res_dot_B, gradPhi(element,i))
			# res_final_dot = np.dot(np.transpose(gradPhi(element,j)),res_dot_tmp)
			val = 2*element.area()*res_final_dot*element.area()
			# print(val[0][0])
			triplets.append(I,J,val[0][0])

	return triplets


def Stiffness(msh, dim:int, physical_tag:int, t):

	triangles = msh.getElements(dim,physical_tag)
	for i in range(0,len(triangles)):
		stiffness_elem(triangles[i], t)

	return 0

################################
# Quadratures
################################
def phiRef(i:int, param):
	if i==0:
		return 1-param[0]-param[1]
	if i==1:
		return param[0]
	if i==2:
		return param[1] 
	return 0

def interpol_geo(element,m):
	res = [0,0]
	_, pts_param, pts_phys=element.gaussPoint()
	for i in range(0,len(element.points)):
		# x du point de gauss passage de xsi, eta à x, y
		res[0]=res[0] + phiRef(i, pts_param[m])*pts_phys[i][0]
		res[1]=res[1] + phiRef(i, pts_param[m])*pts_phys[i][1]
	return res



def Integrale(msh, dim:int, physical_tag:int, f, B, order=2):
	elements = msh.getElements(dim,physical_tag)
	# print(B[0])
	for ind_elem in range(0,len(elements)):
		# print(elements[ind_elem])
		w, pts_param, pts_phys = elements[ind_elem].gaussPoint()
		res=0
		for i in range(0,3):
			I = loc2glob(elements[ind_elem],i)
			for m in range(0,len(w)):
				x_interpol=interpol_geo(elements[ind_elem],m)
				# somme ( somme ( somme (pds_m * f(x(xsi_m,eta_m)* phi(xsi_m,eta_m)))))
				res=res+w[m]*elements[ind_elem].area()*2*f(x_interpol[0],x_interpol[1])*phiRef(i, pts_param[m])
			# print(type(elements[ind_elem].points[i].id))
			B[I] = B[I] +res

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

	pts=msh.getPoints(dim,physical_tag)

	for ind_s in range(0,len(pts)):
		sommet = pts[ind_s].id

		for j in range(0,triplets.size_data):
			#all val = triplets.data[0]
			#all I = triplets.data[1][0]
			#all J = triplets.data[1][1]
			I=triplets.data[1][0][j]
			J=triplets.data[1][1][j]
			if sommet==I:
				triplets.data[0][j]=0
			B[sommet] = g(pts[ind_s].x,pts[ind_s].y)

		triplets.append(sommet,sommet,1)

	return 0

