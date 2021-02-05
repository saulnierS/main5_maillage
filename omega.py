import gmsh
import sys
import math
import omega
#################################
# maillage de omega avec gmsh
#################################

gmsh.initialize(sys.argv)
		
# Create a model and name it "omega"
model = gmsh.model
model.add("omega")
def maillage_gmsh(h=1):

	# Ask GMSH to display information in the terminal
	gmsh.option.setNumber("General.Terminal", 1)
	
	# Parameters
	R1 = 1  # Radius
	#h = 1# Mesh size
	
	# Create Point 
	p1 = model.geo.addPoint(0,0,0, h)
	p2 = model.geo.addPoint(0,5,0, h)
	p3 = model.geo.addPoint(5,5,0, h)
	p4 = model.geo.addPoint(5,0,0, h)
	
	
	# Create line
	l1 = model.geo.addLine(p1,p2)
	l2 = model.geo.addLine(p2,p3)
	l3 = model.geo.addLine(p3,p4)
	l4 = model.geo.addLine(p4,p1)
	
	# Curveloop and Surface
	curveloop = model.geo.addCurveLoop([l1,l2,l3,l4])
	omega = model.geo.addPlaneSurface([curveloop])
	
	# Physical groups
	# gmsh.model.addPhysicalGroup(dim, list of tags, physical tag)
	gmsh.model.addPhysicalGroup(1, [l2,l3,l4], 1)
	gmsh.model.addPhysicalGroup(1, [l1], 3)
	gmsh.model.addPhysicalGroup(2, [omega], 2)
	
	# This command is mandatory and synchronize CAD with GMSH Model. The less you launch it, the better it is for performance purpose
	gmsh.model.geo.synchronize()
	# Mesh (2D)
	model.mesh.generate(2)
	# Write on disk
	gmsh.write("omega.msh")
	# print(gmsh.model.mesh.getNodes())
	# Launch the GUI (not mandatory at all)
	# gmsh.fltk.run();
	return gmsh

gmsh.finalize()