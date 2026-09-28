import os
import subprocess
import numpy as np
import math as math
import matplotlib.pyplot as plt
import openEMS
from CSXCAD import ContinuousStructure,AppCSXCAD_BIN
#from easyMesh import GenerateMesh, enhance_csx_for_auto_mesh, enhance_FDTD_for_auto_mesh


#constant
f_start=2e9
f_stop=4e9
CuThick=0.0175
gl_dim=[50,50,1.5]
gp_dim=[50,50,CuThick]
gap_width=0.5
feed_width=1
inner_rad=1.0
outer_rad=2.0
observation_plane=0.2
f0=(f_start+f_stop)/2
fc=(f_stop-f_start)/2
speratral_res=[0.1,0.1]
unit=1e-3

CSX = ContinuousStructure()
FDTD = openEMS.openEMS(EndCriteria=1e-4)
FDTD.SetCSX(CSX)
FDTD.SetGaussExcite((f_start+f_stop)/2, (f_stop-f_start)/2)
FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8'])


#mount
ground = CSX.AddMaterial('copper',kappa=2.62e6)
ground.AddBox(start=[-gp_dim[0]/2,-gp_dim[1]/2, -gl_dim[2]],stop=[gp_dim[0]/2, gp_dim[1]/2, -gl_dim[2]-gp_dim[2]])

substrate = CSX.AddMaterial('FR4', epsilon=4.3)
sub=substrate.AddBox(start=[-gl_dim[0]/2,-gl_dim[1]/2,0.0] , stop=[gl_dim[0]/2, gl_dim[1]/2, -gl_dim[2]])
sub.SetPriority(10)

#Feedline
strip1=ground.AddBox(start=[-gp_dim[0]/2,-feed_width/2, 0],stop=[-gap_width/2, feed_width/2, CuThick])
strip2=ground.AddBox(start=[gp_dim[0]/2,-feed_width/2, 0],stop=[gap_width/2, feed_width/2, CuThick])
strip1.SetPriority(20)
strip1.SetPriority(20)

#Circular
phi = np.linspace(0, 2*np.pi, 50)
cyl1= ground.AddPolygon(np.vstack((outer_rad * np.cos(phi), outer_rad * np.sin(phi)+(feed_width/2+(outer_rad-inner_rad)))), norm_dir=2, elevation=CuThick)
cyl1.SetPriority(20)

#etched
etch_material = CSX.AddMaterial('etched_regions', epsilon=1.0)
cut_cyl = etch_material.AddCylinder(start=[0.0, feed_width/2+(outer_rad-inner_rad),0.0],stop=[0.0, feed_width/2+(outer_rad-inner_rad), CuThick],radius=inner_rad)
cut_cyl.SetPriority(30)
cut_box = etch_material.AddBox(start=[-gap_width/2.0, feed_width/2+inner_rad, 0.0], stop=[gap_width/2.0, -feed_width/2, 0.02])
cut_box.SetPriority(30)

port_impedance = 50
port_1=FDTD.AddLumpedPort(1, port_impedance, [-gl_dim[0]/2, feed_width/2, -gl_dim[2]],[-gl_dim[0]/2,  -feed_width/2,CuThick],'z',excite=True)
port_2=FDTD.AddLumpedPort(2, port_impedance, [gl_dim[0]/2, feed_width/2, -gl_dim[2]],[gl_dim[0]/2,  -feed_width/2,CuThick],'z',excite=False) # Passive port
port = [port_1, port_2]

#mesh
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
epsilon_r = 5.7
mu_r = 1
C0=3e8
v = C0 / math.sqrt(epsilon_r * mu_r)
wavelength = v / f_start / unit
max_res = v / (f0+fc) / 1e-3 / 40
box_dim=[gl_dim[0]+wavelength/4,gl_dim[1]+wavelength/4,gp_dim[2]+gl_dim[2]+CuThick+wavelength/4]
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)

x_pos = np.arange(-inner_rad, inner_rad, speratral_res[0])
mesh.AddLine('x', np.unique(np.concatenate((x_pos, -x_pos))))

y_pos = np.arange(0, inner_rad, speratral_res[1])
mesh.AddLine('y', np.unique(np.concatenate((y_pos, -y_pos))+feed_width/2+inner_rad))  

mesh.AddLine('x', [outer_rad,-outer_rad])
mesh.AddLine('y', np.array([outer_rad,-inner_rad-feed_width,-inner_rad])+feed_width/2+inner_rad)
mesh.AddLine('z',[0.0,observation_plane])
#mesh.SmoothMeshLines('x',max_res/8,1.5)
#mesh.SmoothMeshLines('y',max_res/8,1.5)                         
FDTD.AddEdges2Grid(dirs='all', properties = ground)
FDTD.AddEdges2Grid(dirs='all', properties = substrate)
#FDTD.AddEdges2Grid(dirs='all', properties = etch_material)
mesh.SmoothMeshLines('x',max_res,1.5)
mesh.SmoothMeshLines('y',max_res,1.5)
mesh.SmoothMeshLines('z',max_res,1.5)
mesh.AddLine('x', [box_dim[0]/2,-box_dim[0]/2])
mesh.AddLine('y', [box_dim[1]/2,-box_dim[1]/2])
mesh.AddLine('z', [-gl_dim[2]-gp_dim[2],box_dim[2]])
mesh.SmoothMeshLines('x',max_res,1.5)
mesh.SmoothMeshLines('y',max_res,1.5)
mesh.SmoothMeshLines('z',max_res,1.5)



workspace_root = os.getcwd()
sim_dir = os.path.join(workspace_root,'omega_loop')
csx_file = os.path.join(sim_dir, 'antenna_layout.xml')
if not os.path.exists(sim_dir):
    os.makedirs(sim_dir)
CSX.Write2XML(csx_file)
subprocess.Popen([AppCSXCAD_BIN, csx_file])

#os.environ['PATH'] = r"D:\Obsidian\Study\MS Thesis\Master_code_file_rabi_freq\openEMS"
#FDTD.Run(os.path.join(workspace_root,'omega_loop'),cleanup=True)


















