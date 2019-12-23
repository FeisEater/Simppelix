import numpy as np

#S = np.matrix([[-6.0,-3.0,-3.0,0.0,0.0,0.0],[1.0,3.0,2.0,1.0,0.0,0.0],[2.0,2.0,-1.0,0.0,1.0,0.0],[3.0,-1.0,3.0,0.0,0.0,1.0]])
S = np.matrix([[0.0,1.0,1.0,1.0,1.0,1.0],[1.0,3.0,2.0,1.0,0.0,0.0],[3.0,5.0,1.0,1.0,1.0,0.0],[4.0,2.0,5.0,1.0,0.0,1.0]])
#Following has infinite loop
#S = np.matrix([[0.0,0.0,2.0,0.0,1.0,0.0,0.0,5.0],[4.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0],[2.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0],[3.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0],[6.0,0.0,3.0,1.0,0.0,0.0,0.0,1.0]])
status = "calculating"
SnoZero = S[1:,1:]
colSize = SnoZero.shape[0]
unit = np.zeros(colSize)
basis = np.full(colSize, -1)
for i in range(colSize):
  unit[i] = 1.0
  for j in range(SnoZero.shape[1]):
    if np.all(np.transpose(SnoZero[:,j]) == unit):
      basis[i] = j
  unit[i] = 0
missingBasis = np.count_nonzero(basis < 0)

def coreSimplex(twoPhase=False):
  objRow = 1 if twoPhase else 0
  for idx,col in enumerate(basis):
    S[objRow,:] = S[objRow,:] - S[objRow,col+1]*S[objRow+idx+1,:]
  while True:
    negRow0 = S[objRow,1:] >= 0.0
    if np.all(negRow0):
      status = "optimum"
      break
    #TODO: different strategies for choosing pivot column
    pivotColIdx = 1 + np.where(negRow0 == False)[1][0]
    #End TODO
    pivotCol = S[:, pivotColIdx][:,0]
    if np.all(pivotCol[objRow+1:] <= 0):
      status = "unbounded"
      break
    pivotRowIdx = objRow + 1 + np.argmin(S[objRow+1:,0] / np.maximum(pivotCol[objRow+1:], 0))
    basis[pivotRowIdx-objRow-1] = pivotColIdx-1
    S[pivotRowIdx] = S[pivotRowIdx] / pivotCol[pivotRowIdx]
    for i in range(S.shape[0]):
      if i == pivotRowIdx:
        continue
      S[i] -= pivotCol[i]*S[pivotRowIdx]
    print(S)
    print(basis)

if missingBasis > 0:
  S = np.insert(S, 1, 0, axis=0)
  unit = np.zeros(S.shape[0])
  unit[1] = 1
  j = 0
  for idx,col in enumerate(basis):
    if col >= 0:
      basis[idx] += missingBasis
      continue
    unit[idx+2] = 1
    S = np.insert(S, j+1, unit, axis=1)
    basis[idx] = j
    unit[idx+2] = 0
    j += 1
  coreSimplex(True)
  if S[1,0] > 0:
    status = "feasible set empty"
    print(status)
    exit()
  while not np.all(basis >= missingBasis):
    print("Driving out artificial variables...")
    pivotRowIdx = 2 + np.where(basis < missingBasis)[0][0]
    pivotColChoose = lambda x: x in basis and x != 0
    pivotColIdx = 1 + missingBasis + np.where(pivotColChoose(S[pivotRowIdx,1+missingBasis:]))[1][0]
    S[pivotRowIdx] = S[pivotRowIdx] / S[pivotRowIdx, pivotColIdx]
    for i in range(S.shape[0]):
      if i == pivotRowIdx:
        continue
      S[i] -= S[i, pivotColIdx]*S[pivotRowIdx]    
  S = np.delete(S, 1, 0)
  S = np.delete(S, range(1, 1+missingBasis), 1)
  basis -= missingBasis
coreSimplex()
print(S)
print(basis)
print(status)
