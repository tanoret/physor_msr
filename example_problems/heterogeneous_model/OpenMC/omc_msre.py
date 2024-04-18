'''
MSRE model
2024-03-11
'''

import openmc as omc
import openmc.mgxs as mgxs
import numpy as np
from isotxs import ISOTXS

def PrepareLibrary(domains, domain_type, scatt_order, groups, geometry,
                   by_nuclide=True):
   lib = mgxs.Library(geometry)
   lib.energy_groups = mgxs.EnergyGroups(groups)
   lib.mgxs_types = ["total", "absorption", "nu-fission", "fission", "chi",
                     "consistent nu-scatter matrix", "consistent scatter matrix",
                     "kappa-fission"]
   lib.correction = None
   lib.legendre_order = scatt_order
   lib.domain_type = domain_type
   lib.domains = domains
   lib.by_nuclide = by_nuclide

   lib.build_library()

   return lib

def CreateISOTXS(lib, sp_file, print_it=True):
   sp = omc.StatePoint(sp_file)
   lib.load_from_statepoint(sp)

   n_xn_adjust_absorb = True
   correct_scatt = False
   isotxs_obj = ISOTXS.from_OpenMC_Library(lib, lib.legendre_order,
                                           n_xn_adjust_absorb, correct_scatt)

   if n_xn_adjust_absorb:
      isotxs_obj.write_file("ISOTXS", scatter_format="total")
   else:
      isotxs_obj.write_file("ISOTXS", scatter_format="macro")

   if print_it:
      # Then tell the ISOTXS to write itself to a file
      isotxs_obj.to_csv("xs_data.csv", print_by="reaction")

def fuelBlock(m1,m2,l_cgb,w_cgb):
   '''
   m1 : fuel
   m2 : graphite
   '''
   # dimensional data
   l_cgb_ch = 3.048 # single core channel length
   w_cgb_ch = 1.016 # single core channel width
   r_cgb_ch = 0.508 # single core channel end radius
   l_cgb_sch=l_cgb_ch-r_cgb_ch*2.0 # core graphite block square channel length

   # surfaces
   # --- core graphite block ---
   xl_cgb=omc.XPlane(x0=-l_cgb/2.0)
   xu_cgb=omc.XPlane(x0=l_cgb/2.0 )
   yl_cgb=omc.YPlane(y0=-w_cgb/2.0)
   yu_cgb=omc.YPlane(y0=w_cgb/2.0 )
   # --- core graphite block channel ---
   # left channel
   xu_cgb_lch=omc.XPlane(x0=-l_cgb/2.0+w_cgb_ch/2.0)
   rt_cgb_lch=omc.ZCylinder(x0=-l_cgb/2.0,y0=l_cgb_sch/2.0,r=r_cgb_ch)
   rb_cgb_lch=omc.ZCylinder(x0=-l_cgb/2.0,y0=-l_cgb_sch/2.0,r=r_cgb_ch)
   yu_cgb_vch=omc.YPlane(y0=l_cgb_sch/2.0)
   yl_cgb_vch=omc.YPlane(y0=-l_cgb_sch/2.0)
   # top channel
   yl_cgb_tch=omc.YPlane(y0=w_cgb/2.0-w_cgb_ch/2.0)
   rl_cgb_tch=omc.ZCylinder(x0=-l_cgb_sch/2.0,y0=w_cgb/2.0,r=r_cgb_ch)
   rr_cgb_tch=omc.ZCylinder(x0=l_cgb_sch/2.0,y0=w_cgb/2.0,r=r_cgb_ch)
   xu_cgb_hch=omc.XPlane(x0=l_cgb_sch/2.0)
   xl_cgb_hch=omc.XPlane(x0=-l_cgb_sch/2.0)
   # right channel
   xl_cgb_rch=omc.XPlane(x0=l_cgb/2.0-w_cgb_ch/2.0)
   rt_cgb_rch=omc.ZCylinder(x0=l_cgb/2.0,y0=l_cgb_sch/2.0,r=r_cgb_ch)
   rb_cgb_rch=omc.ZCylinder(x0=l_cgb/2.0,y0=-l_cgb_sch/2.0,r=r_cgb_ch)
   # bottom channel
   yu_cgb_bch=omc.YPlane(y0=-w_cgb/2.0+w_cgb_ch/2.0)
   rl_cgb_bch=omc.ZCylinder(x0=-l_cgb_sch/2.0,y0=-w_cgb/2.0,r=r_cgb_ch)
   rr_cgb_bch=omc.ZCylinder(x0=l_cgb_sch/2.0,y0=-w_cgb/2.0,r=r_cgb_ch)
   
   # region
   # --- core graphite block ---
   # left channel
   rg1 = (+xl_cgb & -xu_cgb_lch & +yl_cgb_vch & -yu_cgb_vch)
   rg2 = (+xl_cgb & -rt_cgb_lch & +yu_cgb_vch)
   rg3 = (+xl_cgb & -rb_cgb_lch & -yl_cgb_vch)
   rg_lch = (rg1 | rg2 | rg3)
   # top channel
   rg1 = (+xl_cgb_hch & -xu_cgb_hch & +yl_cgb_tch & -yu_cgb)
   rg2 = (-yu_cgb & -rl_cgb_tch & -xl_cgb_hch)
   rg3 = (-yu_cgb & -rr_cgb_tch & +xu_cgb_hch)
   rg_tch = (rg1 | rg2 | rg3)
   # right channel
   rg1 = (-xu_cgb & +xl_cgb_rch & +yl_cgb_vch & -yu_cgb_vch)
   rg2 = (-xu_cgb & -rt_cgb_rch & +yu_cgb_vch)
   rg3 = (-xu_cgb & -rb_cgb_rch & -yl_cgb_vch)
   rg_rch = (rg1 | rg2 | rg3)
   # bottom channel
   rg1 = (+xl_cgb_hch & -xu_cgb_hch & +yl_cgb & -yu_cgb_bch)
   rg2 = (+yl_cgb & -rl_cgb_bch & -xl_cgb_hch)
   rg3 = (+yl_cgb & -rr_cgb_bch & +xu_cgb_hch)
   rg_bch = (rg1 | rg2 | rg3)
   # all channels
   rg_ch = ((rg_lch | rg_tch | rg_rch | rg_bch))
   # graphite
   rg_tot = (+xl_cgb & -xu_cgb & +yl_cgb & -yu_cgb)
   rg_gr = (rg_tot & (~rg_ch))

   # cells
   c1=omc.Cell(name='fuel',fill=m1,region=rg_ch)
   c2=omc.Cell(name='graphite',fill=m2,region=rg_gr)
   u0=omc.Universe(cells=[c1,c2])
   return u0

nthread = 78

# define material
# fuel
m1=omc.Material(temperature=900.0)
m1.set_density('g/cm3',2.3243)
m1.add_nuclide('Be9'  ,1.18410E-01,'ao')
m1.add_nuclide('F19'  ,5.94485E-01,'ao')
m1.add_nuclide('Li6'  ,3.07439E-05,'ao')
m1.add_nuclide('Li7'  ,2.63554E-01,'ao')
m1.add_nuclide('U235' ,1.00802E-03,'ao')
m1.add_nuclide('U238' ,2.23617E-03,'ao')
m1.add_nuclide('Zr90' ,1.04320E-02,'ao')
m1.add_nuclide('Zr91' ,2.27496E-03,'ao')
m1.add_nuclide('Zr92' ,3.47733E-03,'ao')
m1.add_nuclide('Zr94' ,3.52397E-03,'ao')
m1.add_nuclide('Zr96' ,5.67727E-04,'ao')
# graphite
m2=omc.Material(temperature=600.0)
m2.set_density('g/cm3',1.86)
m2.add_element('C',1.0,'ao')
m2.add_s_alpha_beta('c_Graphite')
# export material
mAll=omc.Materials([m1,m2])
mAll.export_to_xml()

# define geometry
l_cgb    = 5.080  # lattice length
w_cgb    = 5.080  # lattice width
r_core   = 70.485 # core radius
h_core   = 170.18 # core height
# axial surfaces
zl = omc.ZPlane(z0=0.0, boundary_type='vacuum')
zu = omc.ZPlane(z0=h_core, boundary_type='vacuum')
# fuel lattices
u1 = fuelBlock(m1, m2, l_cgb, w_cgb)
# radial surrounding
r1 = omc.ZCylinder(x0=0.0, y0=0.0, r=l_cgb)
c1 = omc.Cell(name='rblock', fill=m1, region=(-r1))
u2 = omc.Universe(cells=[c1])
# empty control lattice
r1 = 10.0
s1 = omc.ZCylinder(x0=0.0, y0=0.0, r=r1)
c1 = omc.Cell(fill=m1, region=(-s1))
u3 = omc.Universe(cells=[c1])
# core
lat1=omc.RectLattice(name='core_lattice')
lat1.pitch=(l_cgb, w_cgb)
lat1.universes=[[u2,u2,u2,u2,u2,u2,u2,u2,u2,u2,u2,u1,u1, u1, u1,u1,u2,u2,u2,u2,u2,u2,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u2,u2,u2,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u2,u2,u2,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u2,u2,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u2,u2,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u2,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u2,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2,u2,u2],
                [u2,u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2,u2],
                [u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2],
                [u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2],
                [u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2],
                [u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2],
                [u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2],
                [u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u3, u1, u3,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2],
                [u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2],
                [u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u3, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2],
                [u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2],
                [u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2],
                [u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2],
                [u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2],
                [u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2],
                [u2,u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2,u2],
                [u2,u2,u2,u2,u1,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u1,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u1,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u1,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u2,u1,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u1,u2,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u2,u2,u1,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u1,u2,u2,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u2,u2,u2,u1,u1,u1,u1,u1, u1, u1,u1,u1,u1,u1,u2,u2,u2,u2,u2,u2,u2,u2],
                [u2,u2,u2,u2,u2,u2,u2,u2,u2,u2,u2,u1,u1, u1, u1,u1,u2,u2,u2,u2,u2,u2,u2,u2,u2,u2,u2]]
lat1.lower_left=(-27.0/2.0*l_cgb, -27.0/2.0*w_cgb)
lat1.outer=u2
rc = omc.ZCylinder(x0=0.0,y0=0.0,r=r_core,boundary_type='vacuum')
c1 = omc.Cell(fill=lat1,region=(-rc & +zl & -zu))
u0 = omc.Universe(cells=[c1])
# export geometry
geom=omc.Geometry(u0)
geom.export_to_xml()

# get tally cells
talcl=[]
tmp1=u1.get_all_cells()
talcl += [tmp1[x] for x in tmp1]
tmp1=u2.get_all_cells()
talcl += [tmp1[x] for x in tmp1]
tmp1=u3.get_all_cells()
talcl += [tmp1[x] for x in tmp1]

# setting
nbatch = 200
calset = omc.Settings()
calset.source    = omc.Source(space=omc.stats.Point(xyz=(0.0,0.0,h_core/2)))
calset.particles = 100000
calset.batches   = nbatch
calset.inactive  = 40
corner1 = [-r_core,-r_core,0.0]
corner2 = [ r_core, r_core,h_core]
volcal  = omc.VolumeCalculation(talcl,10000000,corner1,corner2)
calset.volume_calculations = [volcal]
calset.export_to_xml()

# plot
plt1          = omc.Plot()
plt1.basis    = 'xy'
plt1.origin   = (0.0,0.0,h_core/2)
plt1.width    = (r_core*2.0,r_core*2.0)
plt1.pixels   = (8000,8000)
plt1.color_by = 'material'
plt1.colors   = {m1:(255,0,0),m2:(0,255,0)}
plts = omc.Plots([plt1])
plts.export_to_xml()

# tally
energy=(1.0000E-05,7.3000E-01,2.9023E+01,9.1188E+03,2.0000E+07)
mgxs_lib = PrepareLibrary(talcl,'cell',1,energy,geom,True)
tals = omc.Tallies()
mgxs_lib.add_to_tallies_file(tals,merge=True)
tals.export_to_xml()

# run code
omc.plot_geometry()
omc.calculate_volumes(threads=nthread)
omc.run(threads=nthread,geometry_debug=False)
vol_c=omc.VolumeCalculation.from_hdf5('volume_1.h5')
for c in talcl:
   c.add_volume_information(vol_c)
CreateISOTXS(mgxs_lib,f'statepoint.{nbatch}.h5',True)
