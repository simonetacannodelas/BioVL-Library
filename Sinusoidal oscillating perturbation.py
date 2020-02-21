#This is the function to add perturbations inside a mechanistic model with concentrations

#import the model

def disturbances(modelInUse):
  modelInUse.solve()
  t = modelInUse.solve()[0]
  C= modelInUse.solve()[1]

  #Initialization of the vectors for the superposition of the perturbations
  PV =[]
  PV2 = []
  PV3 = []
  
  #Creation of the vector to collect the data with the noise/perturbation
  C_noise =  np.zeros((C.shape))
  import numpy as np
  
  #sinusoidal oscillation of each time - Dependance of time
  for i in range(len(t)):
      PV.append((((math.sin(t[i]*83))/47)+(math.cos(t[i]*11))/23))
      PV2.append(0.6*(((math.sin(t[i] * 61)/31) + (math.cos(t[i] * 43)) / 61)))
      PV3.append(0.33*(((math.sin(t[i] * 41))/19) + (math.cos(t[i] * 61)) / 89))
  
  #stablishing the perturbation inside the concentration
  for i in range(len(C[0])):
      C_noise [:, i] = C[:, i] + C[:, i]* PV + C[:, i]*PV2 + C[:, i]*PV3

  return C_noise
