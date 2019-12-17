import numpy as np

S = np.matrix([[-6.0,-3.0,-3.0,0.0,0.0,0.0],[1.0,3.0,2.0,1.0,0.0,0.0],[2.0,2.0,-1.0,0.0,1.0,0.0],[3.0,-1.0,3.0,0.0,0.0,1.0]])
status = "calculating"
while True:
  negRow0 = S[0,1:] >= 0.0
  if np.all(negRow0):
    status = "optimum"
    break
  #TODO: different strategies for choosing pivot column
  pivotColIdx = 1 + np.where(negRow0 == False)[1][0]
  #End TODO
  pivotCol = S[:, pivotColIdx][:,0]
  if np.all(pivotCol[1:] <= 0):
    status = "unbounded"
    break
  pivotRowIdx = 1 + np.argmin(S[1:,0] / np.maximum(pivotCol[1:], 0))
  S[pivotRowIdx] = S[pivotRowIdx] / pivotCol[pivotRowIdx]
  for i in range(S.shape[0]):
    if i == pivotRowIdx:
      continue
    S[i] -= pivotCol[i]*S[pivotRowIdx]
  print(S)
print(status)