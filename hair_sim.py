import meep as mp
import math
import numpy as np
import hair_material as mat
import multiprocessing
import visualize as vs
import matplotlib.pyplot as plt

resolution = 35  
sxy = 100 # length of simulation domain excluding PMLs
dft_mon = 4 # length of DFT box
dpml = 2.0 # PML thickness 
rot_angle = np.radians(0) # rotation angle of the source

cell = mp.Vector3(sxy+2*dpml,sxy+2*dpml,0)
pml_layers = [mp.PML(dpml)]

 # ------------------ Source parameters ----------------------- #

fcen = 0.5 # center frequency at 1 micron

def do_sim(fcen):
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

    # ------------------ Formulate simulation -------------------- #

    sim = mp.Simulation(cell_size=cell,
                        resolution=resolution,
                        sources=sources,
                        boundary_layers=pml_layers)

    nearfield_box = sim.add_near2far(fcen, 0, 1,
                                    mp.Near2FarRegion(mp.Vector3(y=0.5*dft_mon), size=mp.Vector3(dft_mon)),
                                    mp.Near2FarRegion(mp.Vector3(y=-0.5*dft_mon), size=mp.Vector3(dft_mon), weight=-1),
                                    mp.Near2FarRegion(mp.Vector3(0.5*dft_mon), size=mp.Vector3(y=dft_mon)),
                                    mp.Near2FarRegion(mp.Vector3(-0.5*dft_mon), size=mp.Vector3(y=dft_mon), weight=-1))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(x=-40), 1e-4))
    eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
    ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
    # ------------------ Save flux data ------------------------- #

    # for normalization run, save flux fields data
    nearfield_refl_data = sim.get_near2far_data(nearfield_box)
    sim.reset_meep()

    r = 25 # radius of cylinder
    geometry = mat.get_geometry()
    

    sim = mp.Simulation(cell_size=cell,
                        resolution=resolution,
                        sources=sources,
                        geometry = geometry,
                        boundary_layers=pml_layers)


    nearfield_box = sim.add_near2far(fcen, 0, 1,
                                    mp.Near2FarRegion(mp.Vector3(y=0.5*dft_mon), size=mp.Vector3(dft_mon)),
                                    mp.Near2FarRegion(mp.Vector3(y=-0.5*dft_mon), size=mp.Vector3(dft_mon), weight=-1),
                                    mp.Near2FarRegion(mp.Vector3(0.5*dft_mon), size=mp.Vector3(y=dft_mon)),
                                    mp.Near2FarRegion(mp.Vector3(-0.5*dft_mon), size=mp.Vector3(y=dft_mon), weight=-1))

    # for normal run, load negated fields to subtract incident from total fields
    sim.load_minus_near2far_data(nearfield_box, nearfield_refl_data)
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(x=-40), 1e-4))
    eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
    ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)



    r = 1000/fcen      # half side length of far-field square box OR radius of far-field circle
    res_ff = 1         # resolution of far fields (points/μm)
    far_flux_box = (nearfield_box.flux(mp.Y, mp.Volume(center=mp.Vector3(y=r), size=mp.Vector3(2*r)), res_ff)[0]
                - nearfield_box.flux(mp.Y, mp.Volume(center=mp.Vector3(y=-r), size=mp.Vector3(2*r)), res_ff)[0]
                + nearfield_box.flux(mp.X, mp.Volume(center=mp.Vector3(r), size=mp.Vector3(y=2*r)), res_ff)[0]
                - nearfield_box.flux(mp.X, mp.Volume(center=mp.Vector3(-r), size=mp.Vector3(y=2*r)), res_ff)[0])

    npts = 360         # number of points in [0,2*pi) range of angles
    angles = 2*math.pi/npts*np.arange(npts)

    E = np.zeros((npts,3),dtype=np.complex128)
    H = np.zeros((npts,3),dtype=np.complex128)
    for n in range(npts):
        ff = sim.get_farfield(nearfield_box,
                            mp.Vector3(r*math.cos(angles[n]),
                                        r*math.sin(angles[n])))
        E[n,:] = [np.conj(ff[j]) for j in range(3)]
        H[n,:] = [ff[j+3] for j in range(3)]


    Px = np.real(np.multiply(E[:,1],H[:,2])-np.multiply(E[:,2],H[:,1]))
    Py = np.real(np.multiply(E[:,2],H[:,0])-np.multiply(E[:,0],H[:,2]))
    Pr = np.sqrt(np.square(Px)+np.square(Py))
    far_flux_circle = np.sum(Pr)*2*np.pi*r/len(Pr)
    return Pr


for wave_length in range(380,780,40):
    Pr = do_sim(wave_length*0.001);
    np.savetxt("sim_data/data_"+str(wave_length)+"nm.txt",Pr/max(Pr));
