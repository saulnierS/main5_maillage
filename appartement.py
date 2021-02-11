
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
	window_size = 1
	radiateur_size = 1
	wall_size = 0.5

	# Create Point


	points=[]

	# left
	p1=model.geo.addPoint(0,0,0, h)
	p2=model.geo.addPoint(0,2,0, h)
	p3=model.geo.addPoint(0,4,0, h)
	p4=model.geo.addPoint(0,6,0, h)
	p5=model.geo.addPoint(4,6,0, h)
	p6=model.geo.addPoint(4,6.5,0, h)
	p7=model.geo.addPoint(0,6.5,0, h)
	p8=model.geo.addPoint(0,8,0, h)
	p9=model.geo.addPoint(0,10,0, h)
	p10=model.geo.addPoint(0,12,0, h)

	#top
	p11=model.geo.addPoint(2,12,0, h)
	p12=model.geo.addPoint(4,12,0, h)
	p13=model.geo.addPoint(4.5,12,0, h)
	p14=model.geo.addPoint(4.5,11,0, h)
	p15=model.geo.addPoint(4.5,9,0, h)
	p16=model.geo.addPoint(4.5,0.5,0, h)
	p17=model.geo.addPoint(5,0.5,0, h)
	p18=model.geo.addPoint(5,12,0, h)
	p19=model.geo.addPoint(7,12,0, h)
	p20=model.geo.addPoint(8,12,0, h)
	p21=model.geo.addPoint(10,12,0, h)
	

	#right
	p22=model.geo.addPoint(10,8,0, h)
	p23=model.geo.addPoint(9.5,8,0, h)
	p24=model.geo.addPoint(8.5,8,0, h)
	p25=model.geo.addPoint(5.5,8,0, h)
	p26=model.geo.addPoint(5.5,7.5,0, h)
	p27=model.geo.addPoint(10,7.5,0, h)
	p28=model.geo.addPoint(10,4.5,0, h)
	p29=model.geo.addPoint(5.5,4.5,0, h)
	p30=model.geo.addPoint(5.5,4,0, h)
	p31=model.geo.addPoint(10,4,0, h)
	p32=model.geo.addPoint(10,3,0, h)
	p33=model.geo.addPoint(10,2,0, h)
	p34=model.geo.addPoint(10,0,0, h)

	#bottom
	p35=model.geo.addPoint(4,0,0, h)
	p36=model.geo.addPoint(2,0,0, h)

	
	m=[]#mur
	f=[]#fenetre
	r=[]#radiateur

	#left
	m.append(model.geo.addLine(p1,p2))
	f.append(model.geo.addLine(p2,p3))
	m.append(model.geo.addLine(p3,p4))
	m.append(model.geo.addLine(p4,p5))
	m.append(model.geo.addLine(p5,p6))
	m.append(model.geo.addLine(p6,p7))
	m.append(model.geo.addLine(p7,p8))
	f.append(model.geo.addLine(p8,p9))
	m.append(model.geo.addLine(p9,p10))

	#top
	m.append(model.geo.addLine(p10,p11))
	f.append(model.geo.addLine(p11,p12))
	m.append(model.geo.addLine(p12,p13))
	m.append(model.geo.addLine(p13,p14))
	r.append(model.geo.addLine(p14,p15))
	m.append(model.geo.addLine(p15,p16))
	m.append(model.geo.addLine(p16,p17))
	m.append(model.geo.addLine(p17,p18))
	m.append(model.geo.addLine(p18,p19))
	f.append(model.geo.addLine(p19,p20))
	m.append(model.geo.addLine(p20,p21))
	
	
	#right
	m.append(model.geo.addLine(p21,p22))
	m.append(model.geo.addLine(p22,p23))
	r.append(model.geo.addLine(p23,p24))
	m.append(model.geo.addLine(p24,p25))
	m.append(model.geo.addLine(p25,p26))
	m.append(model.geo.addLine(p26,p27))
	m.append(model.geo.addLine(p27,p28))
	m.append(model.geo.addLine(p28,p29))
	m.append(model.geo.addLine(p29,p30))
	m.append(model.geo.addLine(p30,p31))
	m.append(model.geo.addLine(p31,p32))
	r.append(model.geo.addLine(p32,p33))
	m.append(model.geo.addLine(p33,p34))

	#bottom
	m.append(model.geo.addLine(p34,p35))
	r.append(model.geo.addLine(p35,p36))
	m.append(model.geo.addLine(p36,p1))

	# Create line
	# lines=[]
	# for i in range(0,len(points)-1):
	# 	lines.append(model.geo.addLine(points[i],points[i+1]))

	
	# Curveloop and Surface
	curveloop = model.geo.addCurveLoop(m+r+f)
	omega = model.geo.addPlaneSurface([curveloop])
	
	# Physical groups
	# gmsh.model.addPhysicalGroup(dim, list of tags, physical tag)
	gmsh.model.addPhysicalGroup(1, m, 1)
	gmsh.model.addPhysicalGroup(1, r, 2)
	gmsh.model.addPhysicalGroup(1, f, 3)
	gmsh.model.addPhysicalGroup(2, [omega], 4)
	
	# This command is mandatory and synchronize CAD with GMSH Model. The less you launch it, the better it is for performance purpose
	gmsh.model.geo.synchronize()
	# Mesh (2D)
	model.mesh.generate(2)
	# Write on disk
	gmsh.write("appartement.msh")
	# print(gmsh.model.mesh.getNodes())
	# Launch the GUI (not mandatory at all)
	# gmsh.fltk.run();
	return gmsh

gmsh.finalize()