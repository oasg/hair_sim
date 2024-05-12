import matplotlib.pyplot as plt
import numpy as np
import math
def show_model(eps_data,ez_data):
    plt.figure()
    plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
    plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
    plt.axis('off')
    plt.show()
    
def show_Pr(Pr,name):
    plt.figure(figsize=(10,6))
    npts = 360
    angles = 2*math.pi/npts*np.arange(npts)
    plt.plot(angles*180/np.pi,Pr/max(Pr),'b-')
    plt.grid(True)
    plt.title(name)
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Scattering (a.u.)')
    plt.savefig('visualize/angles'+name+'.png')
    plt.show()