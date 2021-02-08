######################################
# resolution / solveur elements finis
######################################

#system
import sys
#maths
import math
import numpy as np
from scipy import sparse
from scipy.sparse import linalg
#class and python files  
import common
import maillage
import fem_p1
import gmsh
import matplotlib.pyplot as plt
import omega


gmsh.initialize(sys.argv)
#data
def g(x,y):
  return np.sin(np.pi*x)*np.sin(np.pi*y)

def f(x,y):
  return g(x,y)*(2*np.pi*np.pi +1 )

def diri(x,y):
  return 0

# Maillage
msh = maillage.Mesh()
om=omega.maillage_gmsh(h=0.5)
msh.GmshToMesh("omega.msh",om)

# Triplets
tM = common.Triplets()
tD = common.Triplets()
fem_p1.Mass(msh, 2, 2, tM)
verif_mass=np.ones(msh.Npts)

fem_p1.Stiffness(msh, 2, 2, tD)
b = np.zeros((msh.Npts))
fem_p1.Integrale(msh, 2, 2, f, b, 2)
# fem_p1.Dirichlet(msh, dim=1, physical_tag=3, B=b, triplets=tM, g=diri)


# # Résolution
# A= fem_p1.build_matrix(t.data)
AM = (sparse.coo_matrix(tM.data)).tocsr()
AD = (sparse.coo_matrix(tD.data)).tocsr()
print(verif_mass)
vect1=np.matmul(verif_mass,AM)
print(np.matmul(vect1,np.transpose(verif_mass)))

UM = linalg.spsolve(AM, b)
UD = linalg.spsolve(AD, b)
print("A masse")
print(AM.toarray())
print("A rigidité")
print(AD.toarray())
print("U masse")
print(UM)
print("U rigidité")
print(UD)

# # Visualisation
x= [pt.x for pt in msh.points]
y= [pt.y for pt in msh.points]
connectivity=[]
for tri in msh.triangles:
  connectivity.append([ p.id for p in tri.points]) 

plt.figure("U")
# print(U)
plt.tricontourf(x, y, connectivity, U, 12)
plt.colorbar()

plt.figure("Uref")
# ### U de référence
Uref = np.zeros((msh.Npts))
for pt in msh.points:
  I = int(pt.id)
  Uref[I] = g(pt.x, pt.y)

# print("U",U)
# print("Uref",Uref)
# print("b",b)
plt.tricontourf(x, y, connectivity, Uref, 12)
plt.colorbar()
plt.show()
gmsh.finalize()