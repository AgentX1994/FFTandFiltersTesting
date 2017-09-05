#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

samplerate = 0
in_data = []
out_data = []
with open("data.dat", 'r') as f:
    samplerate = int(f.readline())
    in_data = [float(x) for x in f.readline().split()]
    out_data = [float(x) for x in f.readline().split()]


print(samplerate)
x1 = np.array(in_data)
x2 = np.array(out_data)
f = np.array(range(len(x1)))
f = f*samplerate/len(x1)
fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(f, x1)
#ax1.set_xlim(0,len(x2));
#ax1.set_ylim(-1,1);

ax2.plot(f, x2)
#ax2.set_xlim(0,len(x2));
#ax2.set_ylim(-1,1);

fig.tight_layout()
plt.show()
