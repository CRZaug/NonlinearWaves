"""

This code compares two different types of Snodgrass plots: The black
dots representing ridge cuts on figures such as figure 19, and the smooth
curve shown in figures such as figure 21

(Aug 1 Honolulu)

"""

import numpy as np
import matplotlib.pyplot as plt

# Black dot ridge cuts 
x = [40.87491819371728,
46.322603075916234,
51.913857984293195,
56.16758017015707,
61.43825670811518,
67.25161567408378,
72.19588514397907,
77.68958606020942]

y = [ -28.046875000000007,
-9.687500000000007,
0.07812499999999645,
8.281249999999998,
13.9453125,
15.1171875,
12.96875,
10.624999999999998]


# smooth data (digitized by Diane)
snod = np.transpose(np.loadtxt('Aug1Data/Raw/gauge4.out'))

#Plot both to compare
plt.plot(snod[0],snod[1],'.')
plt.plot(x,y,'*',markersize = 10)
plt.show()