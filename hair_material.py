import meep as mp
import math
import numpy as np
import matplotlib.pyplot as plt

def get_geometry():
    r1 = 25 # radius of cylinder

    cuticle_thickness = 0.5
    cmc_thickness = 0.1
    num_of_layer = 6
    index_cuticle = 1.55
    index_cmc = 1.42
    index_hair = 1.59



    geometry = []
    for i in range(num_of_layer-1,0,-1):
        base_r = r1 + i*(cmc_thickness+cuticle_thickness);
        geometry.append(mp.Cylinder(center=mp.Vector3(),radius = base_r+cmc_thickness+cuticle_thickness, material=mp.Medium(index = index_cuticle), height=mp.inf));
        geometry.append(mp.Cylinder(center=mp.Vector3(),radius = base_r+cmc_thickness, material=mp.Medium(index = index_cmc), height=mp.inf));
        
    geometry.append(mp.Cylinder(center=mp.Vector3(),radius = r1+cmc_thickness+cuticle_thickness, material=mp.Medium(index = index_cuticle), height=mp.inf));
    geometry.append(mp.Cylinder(center=mp.Vector3(),radius = r1+cmc_thickness, material=mp.Medium(index = index_cmc), height=mp.inf));
    geometry.append(mp.Cylinder(center=mp.Vector3(),radius = r1, material=mp.Medium(index = index_hair), height=mp.inf))
    return geometry



def get_model_plot():
    resolution = 50  
    sxy = 100 # length of simulation domain excluding PMLs
    dft_mon = 4 # length of DFT box
    dpml = 2.0 # PML thickness 
    rot_angle = np.radians(0) # rotation angle of the source

    cell = mp.Vector3(sxy+2*dpml,sxy+2*dpml,0)
    pml_layers = [mp.PML(dpml)]

    # ------------------ Source parameters ----------------------- #

    fcen = 0.5 # center frequency at 1 micron
    df = 0.05 # bandwidth of gaussian source
    k_point = mp.Vector3(fcen*1).rotate(mp.Vector3(z=1), rot_angle)
    sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                                center=mp.Vector3(x=-40),
                                size=mp.Vector3(y=(sxy + 2*dpml)),
                                direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION,
                                eig_kpoint=k_point,
                                eig_band=1,
                                eig_parity=mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
                                eig_match_freq=True)]
    sim = mp.Simulation(cell_size=cell,
                    resolution=resolution,
                    sources=sources,
                    geometry = get_geometry(),
                    boundary_layers=pml_layers)
    sim.run(until = 1)
    eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
    plt.figure()
    sim.plot2D()
    plt.show()