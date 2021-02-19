#****************************************
# solve the equation of the project
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
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.pyplot as plt
# maillage
import gmsh
# files  
import common
import maillage
import fem_p1
import omega
import appartement as appart


gmsh.initialize(sys.argv)
#------------------------
# Data project
# wall 		 dim 1 tag 1
# radiator   dim 1 tag 2
# window     dim 1 tag 3
# omega      dim 2 tag 4
#-------------------------
# dirichlet = radiator's temperature
def T_radiator(x,y):
  return 25
# dirichlet = temperature outside
def T_window(x,y):
  return -10
 
#------------
# Option
#------------
argv = sys.argv[1:]
h = float(argv[0]) if len(argv) > 0 else 0.1

#-----------
# Mesh part
#-----------
msh = maillage.Mesh()
omega=appart.maillage_gmsh(h=h)
msh.GmshToMesh("omega.msh",omega)

#-------------
# Solveur FEM
#-------------
t = common.Triplets()
# mass matrix
fem_p1.Mass(msh=msh, dim=2, physical_tag=4, triplets=t)
# stiffness matrix
fem_p1.Stiffness(msh=msh, dim=2, physical_tag=4, triplets=t)

b = np.zeros((msh.Npts))
# dirichlet application
fem_p1.Dirichlet(msh=msh, dim=1, physical_tag=2, B=b, triplets=t, g=T_radiator)
fem_p1.Dirichlet(msh=msh, dim=1, physical_tag=3, B=b, triplets=t, g=T_window)

#-------------
# Solve
#-------------
A = (sparse.coo_matrix(t.data)).tocsr()
U = linalg.spsolve(A, b)


#---------------
# visualization
#---------------
x= [pt.x for pt in msh.points]
y= [pt.y for pt in msh.points]
connectivity=[]
for tri in msh.triangles:
  connectivity.append([ p.id for p in tri.points]) 

plt.figure("U")
plt.tricontourf(x, y, connectivity, U, 12)
plt.colorbar()
plt.show()

gmsh.finalize()