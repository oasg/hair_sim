import matplotlib.pyplot as plt
import numpy as np
import math
import visualize as vs
data = []
for wave_length in range(380,780,40):
    Pr = np.loadtxt("sim_data/data_"+str(wave_length)+"nm.txt");
    length = len(Pr)
    split_index = length // 2
    nPr = Pr[split_index:] + Pr[:split_index]
    data.append(Pr)


for i in range(10):
   vs.show_Pr(data[i],str(380+i*40)+"nm") 


npts = 360
angles = 2*math.pi/npts*np.arange(npts)

plt.figure(figsize=(20, 12))
for i in range(10):
    Pr = data[i]
    plt.plot(angles*180/np.pi,Pr/max(Pr),label = str(i*40+380)+'nm')
plt.grid(True)
plt.title('Hair scattering intensity')
plt.xlabel('Angle (degrees)')
plt.ylabel('Scattering (a.u.)')
plt.legend()
plt.savefig('visualize/angles'+"Hair_scattering_intensity"+'.png')
plt.show()
    
