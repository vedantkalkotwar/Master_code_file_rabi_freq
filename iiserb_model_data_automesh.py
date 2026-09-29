import os
import h5py
from parameters import *
import subprocess
import numpy as np
import math as math
import matplotlib.pyplot as plt
import openEMS
from openEMS.physical_constants import *
#from openEMS.utilities import HDF5Dump
from CSXCAD import ContinuousStructure,AppCSXCAD_BIN
from easyMesh import GenerateMesh, enhance_csx_for_auto_mesh, enhance_FDTD_for_auto_mesh




CSX = ContinuousStructure()
FDTD = openEMS.openEMS(EndCriteria=1e-4)
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
v = C0 / math.sqrt(epsilon_r * mu_r)
wavelength = v / f_start / unit
max_res = v / (f0+fc) / 1e-3 / 40
box_dim=[gl_dim[0]+wavelength/4,gl_dim[1]+wavelength/4,gp_dim[2]+gl_dim[2]+wavelength/4]
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
primitives_mesh_setup = {}
properties_mesh_setup = {}
global_mesh_setup = {
    'drawing_unit': unit,
    'start_frequency': f_start,
    'stop_frequency': f_stop,
    'mesh_resolution': 'high',  # 'low', 'medium', 'high', 'very_high'
    'smooth_metal_edge': 'False', # useful for thin metal layers, Options: False, 'one_third_two_thirds', 'extra_lines', 
    'use_circle_detection': True,
    'boundary_distance': ['auto', 'auto', 'auto', 'auto', -gl_dim[2], 'auto'], # value, 'auto' or None
    'handle_closely_placed_edges': False,  # if True, then mesher will try to handle close placed edges by merging them
    'refined_cellsize': speratral_res[0],
    #'min_cellsize': 0.01,
    # 'num_lines': 3,
    # 'max_cellsize': max_cellsize,
    'f0': f0,  # center frequency
    'fc': fc,  # 20 dB corner frequency
}
CSX = enhance_csx_for_auto_mesh(CSX, primitives_mesh_setup)
FDTD = enhance_FDTD_for_auto_mesh(FDTD, primitives_mesh_setup)
FDTD.SetGaussExcite((f_start+f_stop)/2, (f_stop-f_start)/2)
FDTD.SetBoundaryCond(['PML_4', 'PML_4', 'PML_4', 'PML_4', 'PEC', 'PML_4'])


#mount
ground = CSX.AddConductingSheet('copper',conductivity=5.8e7, thickness=CuThick*unit)
gnd_plate=ground.AddBox(start=[-gp_dim[0]/2,-gp_dim[1]/2, -gl_dim[2]],stop=[gp_dim[0]/2, gp_dim[1]/2, -gl_dim[2]])
gnd_plate.SetPriority(30)

substrate = CSX.AddMaterial('FR4', epsilon=4.3)
sub=substrate.AddBox(start=[-gl_dim[0]/2,-gl_dim[1]/2,0.0] , stop=[gl_dim[0]/2, gl_dim[1]/2, -gl_dim[2]])
sub.SetPriority(20)

#Feedline
strip1=ground.AddBox(start=[-gp_dim[0]/2,-feed_width/2, 0],stop=[-gap_width/2, feed_width/2, 0])
strip2=ground.AddBox(start=[gp_dim[0]/2,-feed_width/2, 0],stop=[gap_width/2, feed_width/2, 0])
strip1.SetPriority(40)
strip2.SetPriority(40)

#Circular
phi = np.linspace(0, 2*np.pi, 50)
cyl1= ground.AddPolygon(np.vstack((outer_rad * np.cos(phi), outer_rad * np.sin(phi)+(feed_width/2+(outer_rad-inner_rad)))), norm_dir=2, elevation=0)
cyl1.SetPriority(40)

#Coplainer Waveguide
# Co_up=ground.AddBox(start=[-gp_dim[0]/2,gap_co+(feed_width/2)+inner_rad+outer_rad, 0],stop=[gp_dim[0]/2, gl_dim[1]/2, 0])
# Co_down=ground.AddBox(start=[-gp_dim[0]/2,-gap_co-feed_width/2, 0],stop=[gp_dim[0]/2, -gl_dim[1]/2, 0])
# Co_up.SetPriority(40)
# Co_down.SetPriority(40)

#etched
etch_material = CSX.AddMaterial('etched_regions', epsilon=1.0)
#cut_cyl = etch_material.AddCylinder(start=[0.0, feed_width/2+(outer_rad-inner_rad),0.0],stop=[0.0, feed_width/2+(outer_rad-inner_rad), gap_diamond],norm_dir=2,radius=inner_rad)
cut_cyl = etch_material.AddPolygon(np.vstack((inner_rad * np.cos(phi), inner_rad * np.sin(phi)+(feed_width/2+(outer_rad-inner_rad)))), norm_dir=2, elevation=0)
cut_cyl.SetPriority(50)
cut_box = etch_material.AddBox(start=[-gap_width/2.0, feed_width/2+inner_rad, 0.0], stop=[gap_width/2.0, -feed_width/2, gap_diamond])
cut_box.SetPriority(50)

#Diamond
# NV_center = CSX.AddMaterial('NV_center', epsilon=5.7)
# NV=NV_center.AddBox(start=[-NV_dim[0]/2,-NV_dim[1]/2+feed_width/2+inner_rad,0.0] , stop=[NV_dim[0]/2, NV_dim[1]/2+feed_width/2+inner_rad, NV_dim[2]])
# NV.SetPriority(10)

port_impedance = 50
port_1=FDTD.AddLumpedPort(1, port_impedance, [-gl_dim[0]/2, feed_width/2, -gl_dim[2]],[-gl_dim[0]/2,  -feed_width/2,0],'z',excite=True)
port_2=FDTD.AddLumpedPort(2, port_impedance, [gl_dim[0]/2, feed_width/2, -gl_dim[2]],[gl_dim[0]/2,  -feed_width/2,0],'z',excite=False) # Passive port
port = [port_1, port_2]

#feild-dump
hfield = CSX.AddDump('HField', dump_type=11, frequency=f_meas, file_type=1)
hfield.AddBox(start=[-dump_dim[0], -dump_dim[1]+(feed_width/2+inner_rad), 0.0],stop =[ dump_dim[0],  dump_dim[1]+(feed_width/2+inner_rad), dump_dim[2]])

#mesh
mesh.AddLine('z',[0,observation_plane,-gl_dim[2]])
GenerateMesh(CSX, global_mesh_setup, primitives_mesh_setup, properties_mesh_setup)

workspace_root = os.getcwd()
sim_dir = os.path.join(workspace_root,'omega_loop')
csx_file = os.path.join(sim_dir, 'antenna_layout.xml')
if not os.path.exists(sim_dir):
    os.makedirs(sim_dir)
CSX.Write2XML(csx_file)
#subprocess.Popen([AppCSXCAD_BIN, csx_file])


os.environ['PATH'] = r"D:\Obsidian\Study\MS Thesis\Master_code_file_rabi_freq\openEMS"
FDTD.Run(os.path.join(workspace_root,'omega_loop'),cleanup=True)



freq = np.linspace(f0 - fc, f0 + fc, 501)
port[0].CalcPort(os.path.join(workspace_root,'omega_loop'), freq)

Zin = port[0].uf_tot / port[0].if_tot
s11 = port[0].uf_ref / port[0].uf_inc
P0_in = float(np.interp(f_meas, freq, port[0].P_acc))
with open("calculated_params.py", "w") as f:
    f.write(f"P0_in = {P0_in:.15e}\n")
    f.write(f"s11 = {s11.tolist()!r}\n")
fig, ax = plt.subplots()
ax.plot(freq / 1e6, 20 * np.log10(np.abs(s11)), 'k-', lw=2)
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('|S₁₁| (dB)')
ax.set_title('Reflection coefficient S₁₁')
ax.grid(True)

fig, ax = plt.subplots()
ax.plot(freq / 1e6, np.real(1. / Zin), 'k-', lw=2, label='real')
ax.plot(freq / 1e6, np.imag(1. / Zin), 'r--', lw=2, label='imag')
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Admittance Y_in (S)')
ax.set_title('Feed port admittance')
ax.legend()
ax.grid(True)

from matplotlib.widgets import Slider
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

h5_path = os.path.join(workspace_root, 'omega_loop', 'HField.h5')

with h5py.File(h5_path, 'r') as f:
    x = f['Mesh/x'][:] / unit          # (24,)
    y = f['Mesh/y'][:] / unit          # (24,)
    z = f['Mesh/z'][:] / unit          # (3,)
    real = f['FieldData/FD/f0_real'][:]   # (3, Nz, Ny, Nx)
    imag = f['FieldData/FD/f0_imag'][:]

H_all = real + 1j * imag
mu0 = 4 * np.pi * 1e-7
B_all = mu0 * H_all /np.sqrt(P0_in)

def b1(idx):
    Bx = B_all[0, idx, :, :]
    By = B_all[1, idx, :, :]
    Bz = B_all[2, idx, :, :]
    return np.sqrt(np.abs(Bx)**2 + np.abs(By)**2 + np.abs(Bz)**2) * 1e4

X, Y = np.meshgrid(x, y, indexing='xy')
z_idx = np.argmin(np.abs(z))

fig, ax = plt.subplots(figsize=(5, 6))
plt.subplots_adjust(bottom=0.18)

im = ax.pcolormesh(X, Y, b1(z_idx), shading='nearest')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'$|B_1^+|$ (G/$\sqrt{\mathrm{W}}$)')

ax.set_aspect('equal')
ax.set_xlabel('x (mm)')
ax.set_ylabel('y (mm)')
ax.set_xlim(-2, 2)
ax.set_ylim(Y.min(), Y.max())
title = ax.set_title(rf'$B_1^+$ Field, $z={z[z_idx]:.3f}$ mm')

slider_ax = fig.add_axes([0.15, 0.06, 0.70, 0.035])
z_slider = Slider(ax=slider_ax, label='z (mm)',
                  valmin=z.min(), valmax=z.max(), valinit=z[z_idx])

def update(val):
    new_idx = np.argmin(np.abs(z - z_slider.val))
    data = b1(new_idx)
    im.set_array(data.ravel())
    im.set_clim(data.min(), data.max())
    title.set_text(rf'$B$ Field, $z={z[new_idx]:.3f}$ mm')
    fig.canvas.draw_idle()

z_slider.on_changed(update)
plt.show()