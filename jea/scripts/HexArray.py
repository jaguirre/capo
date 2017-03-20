# TEST CODE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon
#
# This code generates a file of antenna positions and a files with baselines, baseline info, and redundancies

import os
import numpy as np

# TEST CODE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon
#
# This code generates a file of antenna positions and a files with baselines, baseline info, and redundancies

import os
import numpy as np

#HARD CODED SETTINGS
Separation = 14.4
hexNum = 11
precisionFactor = 1000000

positions = [];
x = []
y = []
for row in range(hexNum-1,-(hexNum),-1):
	for col in range(2*hexNum-abs(row)-1):
		xPos = ((-(2*hexNum-abs(row))+2)/2.0 + col)*Separation;
		yPos = row*Separation*np.sqrt(3)/2;
                x.append(xPos)
                y.append(yPos)
		positions.append([xPos, yPos, 0])

nAntennas = len(positions)
baselines = []
baselinePairs = []
for ant1 in range(nAntennas):
	for ant2 in range(ant1+1,nAntennas):
		baselines.append((int(np.round(precisionFactor*(positions[ant1][0]-positions[ant2][0]))), int(np.round(precisionFactor*(positions[ant1][1]-positions[ant2][1]))), int(np.round(precisionFactor*(positions[ant1][2]-positions[ant2][2])))))
		baselinePairs.append((ant1, ant2))

baselineDict = {}
for b in baselines:
	if baselineDict.has_key(b):
		baselineDict[b] = baselineDict[b] + 1
	else:
		baselineDict[b] = 1

print "With", len(positions), "antennas there are", len(baselineDict.items()), "unique baselines."

blx = []
bly = []
for i in range(len(baselines)):
    blx.append(baselines[i][0])
    bly.append(baselines[i][1])

pl.plot(x,y,'bo')
pl.axis('equal')
pl.show()
    
if False:
    scriptDirectory = os.path.dirname(os.path.abspath(__file__))
    np.savetxt(scriptDirectory + "/antenna_positions.dat",np.asarray(positions))
    np.savetxt(scriptDirectory + "/all_baselines.dat",np.asarray(baselines)/(1.0*precisionFactor))
    np.savetxt(scriptDirectory + "/all_baseline_pairs.dat",np.asarray(baselinePairs),fmt='%i')
    
    np.savetxt(scriptDirectory + "/unique_baselines.dat",np.asarray([uniqueBaseline[0] for uniqueBaseline in baselineDict.items()])/(1.0*precisionFactor))
    
    np.savetxt(scriptDirectory + "/redundancy.dat",np.asarray([uniqueBaseline[1] for uniqueBaseline in baselineDict.items()]),fmt='%i')
