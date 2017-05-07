# Edited by  Jin Zhang  06/19/2015
# This script can deal with the binary convexhull
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from itertools import compress

CompositionEnthalpy = np.loadtxt('INPUT.txt')
xComposition = CompositionEnthalpy [:, 0].reshape (-1, 1)
yComposition = CompositionEnthalpy [:, 1].reshape (-1, 1)
xyComposition = CompositionEnthalpy [:, 0:2]

Enthalpy = CompositionEnthalpy [:, 2].reshape (-1, 1)
TotalAtom = xComposition + yComposition
EnthalpyPerAtom = Enthalpy / TotalAtom
ratioComposition = yComposition / TotalAtom
TotalNum = ratioComposition.shape [0]
IndexNum = np.arange ( 1,  TotalNum+1, 1)
IndexNum = IndexNum.reshape (-1, 1)

ComRaitoEnthalpy = np.hstack([xyComposition, ratioComposition, EnthalpyPerAtom, IndexNum]) ## put composition, ratio and enthalpy into one matrix
ComRaitoEnthalpyIndex=ComRaitoEnthalpy[ComRaitoEnthalpy[:,2].argsort()] ## sort according to the ratio

xCompositionNew =  ComRaitoEnthalpyIndex [:, 0].reshape (-1, 1)
yCompositionNew =  ComRaitoEnthalpyIndex [:, 1].reshape (-1, 1)

xyRatio = ComRaitoEnthalpyIndex [:, 2].reshape (-1, 1) # find xyRatio
EnthalpyPerAtom = ComRaitoEnthalpyIndex [:, 3].reshape (-1, 1) # find EnthalpyPerAtom
IndexNum = ComRaitoEnthalpyIndex [:, 4].reshape (-1, 1) # find original index number
TotalAtomNew = xCompositionNew + yCompositionNew
EnthalpyPerAtom =  np.hstack ([EnthalpyPerAtom, IndexNum])  # put the xyRation and original index together

##find the minmum enthalpy per atom at two ends
yMinEnthalpy = []; yMaxEnthalpy = []; FormEnthalpy = []
for i in range (0,  TotalNum, 1) :
     if xyRatio [i, 0] == 0 :
        yMinEnthalpy = np.hstack([yMinEnthalpy, EnthalpyPerAtom [i, 0]])
        yMinEnthalpy = min (yMinEnthalpy)
for j in range (0,  TotalNum, 1) :   
      if  xyRatio [j, 0] == 1 :
          yMaxEnthalpy = np.hstack([yMaxEnthalpy, EnthalpyPerAtom [j, 0]])
          yMaxEnthalpy = min (yMaxEnthalpy)
          
##calculating the formation energy at all ratio
for i in range (0,  TotalNum, 1) :
   FormEnthalpy =np.hstack([FormEnthalpy, EnthalpyPerAtom [i, 0] - (xCompositionNew [i, 0] /TotalAtomNew [i, 0])*yMinEnthalpy- (yCompositionNew [i, 0]/TotalAtomNew [i, 0])*yMaxEnthalpy])
FormEnthalpy = FormEnthalpy.reshape(-1,1)


RatioFormEnthalpy = np.hstack([xyRatio, FormEnthalpy]) # put the Ratio formation enthalpy and original index together
RatioFormEnthalpyIndex = np.hstack([xyRatio, FormEnthalpy, IndexNum]) # put the Ratio formation enthalpy and original index together

## using the convexhull function to find the convex hull

hull= ConvexHull (RatioFormEnthalpy)
hull_x = RatioFormEnthalpy[hull.vertices,0]
hull_y = RatioFormEnthalpy[hull.vertices,1]
hull_index = RatioFormEnthalpyIndex[hull.vertices,2]

y_less0 = [n <=  0 for n in hull_y] # remove the enthalpy of formation which more than 0 on the convex hull
hull_Y = list (compress(hull_y, y_less0))
hull_X = list (compress(hull_x, y_less0))
hull_Index = list (compress(hull_index, y_less0))

hull_Y  =np.array(hull_Y) # transform hull_Y from 2d list to 2d numpy array
hull_X  =np.array(hull_X) 
     
for i in range (0, len (hull_Index)):
     hull_Index [i] =int (hull_Index[i]) # transform the final index into integer
     

## plotting the convex hull
plt.plot(RatioFormEnthalpy[:,0], RatioFormEnthalpy[:,1], 'o')     
plt.plot(hull_X, hull_Y, 'go-', lw=2)

#for i, txt in enumerate(hull_Index):
  #  plt.annotate (txt, (hull_X[i], hull_Y[i]))

for i in range (0, len(hull_Index)):
    plt.annotate (hull_Index[i], (hull_X[i], hull_Y[i]), fontsize=12)


print (hull_Index)
print (hull_Y)
print (hull_X)
plt.xlabel ('Composition ratio: B/(A+B)', fontsize=14)
plt.annotate ('A', xy=(-0.008, hull_Y.min()*1.3), xytext=(-0.018, hull_Y.min() * 1.8))
plt.annotate ('B', xy=(1.008, hull_Y.min()*1.3), xytext=(1.004, hull_Y.min() * 1.8))
plt.ylabel ('Enthalpy of formation (eV/atom)', fontsize=14)
plt.title('Convex hull for binary system)')
plt.ylim(hull_Y.min() * 1.3, FormEnthalpy.max() * 1.05)
plt.xlim(-0.008, 1.008)
plt.show()


print (hull_Index)
print (hull_Y)
print (hull_X)






