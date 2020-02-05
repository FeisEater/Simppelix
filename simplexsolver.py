import numpy as np
from enum import Enum, auto

class SolverState(Enum):
    PRESOLVE = auto()
    PHASE1_STEP = auto()
    DRIVEOUT_STEP = auto()
    PHASE2_STEP = auto()
    DONE = auto()

class SimplexSolver:
    
    def __stripReturnChar(self, str):
      if str[-1] == '\n':
        return str[:-1]
      return str

    '''
    Reads MPS file to solve LP problem. Currently does not support RANGES and BOUNDS
    mpsName - file location string of MPS file
    '''
    def readMps(self, mpsName):
      with open(mpsName, encoding='utf16') as f:
        mode = ""
        rowMap = {}
        rowCount = 1
        varCount = 1
        objectiveRow = ""
        for line in f:
          if len(line) == 0:
            continue
          if line[0] != ' ':
            mode = line.strip().upper()
            continue
          if mode == "ROWS":
            boundType = line[0:4].upper()
            rowName = self.__stripReturnChar(line[4:].upper()).strip()
            if boundType != " N  ":
              self.__S = np.insert(self.__S, self.__S.shape[0], 0, axis=0)
              rowMap[rowName] = rowCount
              rowCount += 1
              lessThan = boundType == " L  "
              moreThan = boundType == " G  "
              if lessThan:
                slackVarName = "*slack-{}*".format(self.__slackVars)
                slackVar = 1.0
              elif moreThan:
                slackVarName = "*surplus-{}*".format(self.__slackVars)
                slackVar = -1.0
              if lessThan or moreThan:
                varLen = len(self.__vars) + 1
                self.__vars[slackVarName] = varLen
                self.__S = np.insert(self.__S, self.__S.shape[1], 0, axis=1)
                self.__S[self.__S.shape[0] - 1, self.__S.shape[1] - 1] = slackVar
                self.__slackVars += 1
            else:
              if objectiveRow != "":
                continue
              objectiveRow = rowName
              rowMap[rowName] = 0
          elif mode == "COLUMNS":
            varName = line[4:12].upper().strip()
            if varName not in self.__vars:
              self.__vars[varName] = varCount
              self.__S = np.insert(self.__S, varCount, 0, axis=1)
              varCount += 1
            rowName1 = line[12:22].upper().strip()
            value1 = float(self.__stripReturnChar(line[22:min(36,len(line))].strip()))
            self.__S[rowMap[rowName1], self.__vars[varName]] = value1
            if len(line) >= 47:
              rowName2 = line[36:47].upper().strip()
              value2 = float(self.__stripReturnChar(line[47:].strip()))
              self.__S[rowMap[rowName2], self.__vars[varName]] = value2
          elif mode == "RHS":
            rowName1 = line[12:22].upper().strip()
            value1 = float(self.__stripReturnChar(line[22:min(36,len(line))].strip()))
            self.__S[rowMap[rowName1], 0] = value1
            if len(line) >= 47:
              rowName2 = line[36:47].upper().strip()
              value2 = float(self.__stripReturnChar(line[47:].strip()))
              self.__S[rowMap[rowName2], 0] = value2
          else:
            print("{} not supported".format(mode))
            break
      self.status = "not started"

    '''Split term to coefficient and var'''
    def __splitToCoeffAndVar(self, x):
      split = x.split('*')
      if len(split) == 1:
        x = x.strip()
        if len(x) == 0:
          return ['0', '*null*']
        if x[0] == '-':
          return ['-1', x[1:].strip()]
        return ['1', x]
      elif len(split) == 2:
        if split[0][0] == '-':
          split[0] = '-' + split[0][1:].strip()
        return [split[0].strip(), split[1].strip()]
      else:
        raise ValueError('Illegal amount of multiplication characters')

    '''Add a row to table from string, be it objective row or constraint row'''
    def __setRow(self, str, objective=False):
      #If not objective row, determine equality/inequality sign and bound
      if not objective:
        pair = str.split('=')
        if len(pair) != 2:
          raise ValueError("Malformed constraint expression, should have equality character (=)")
        str = pair[0]
        bound = float(pair[1])
        lessThan = str[-1] == '<'
        moreThan = str[-1] == '>'
        #Remove inequality sign
        if lessThan or moreThan:
          str = str[:-1].strip()
      #Separate linear expression into terms and split coefficients and variables
      coeffAndVar = list( \
        map( \
          lambda x: [float(x[0]), x[1]], \
          map(self.__splitToCoeffAndVar, str.replace("-", "+-").split('+')) \
        ) \
      )
      coeffAndVar = list(filter(lambda x: x[1] != '*null*', coeffAndVar))
      #For one variable constraints, prune infeasible variables
      if not objective:
        varDuplicationRequired = True
        if len(coeffAndVar) == 1:
          posCoeff = coeffAndVar[0][0] >= 0.0
          posBound = bound >= 0.0
          if posCoeff == posBound and posCoeff == moreThan:
            varDuplicationRequired = False
            if "*neg*" + coeffAndVar[0][1] in self.__vars:
              self.__varsToDelete.append(self.__vars["*neg*" + coeffAndVar[0][1]])
          elif posBound == moreThan:
            varDuplicationRequired = False
            if coeffAndVar[0][1] in self.__vars:
              self.__varsToDelete.append(self.__vars[coeffAndVar[0][1]])
        #Simplex table assumes variables are non-negative, so such constraints can be skipped
        if not varDuplicationRequired and bound == 0.0:
          return
      #Split variable to positive and negative variables, add column if not already in table
      for x in coeffAndVar:
        if x[1] not in self.__vars:
          varLen = len(self.__vars) + 2
          self.__vars[x[1]] = varLen - 1
          self.__vars["*neg*" + x[1]] = varLen
          self.__S = np.insert(self.__S, self.__S.shape[1], 0, axis=1)
          self.__S = np.insert(self.__S, self.__S.shape[1], 0, axis=1)
      #Add slack/surplus variable
      if not objective:
        if lessThan:
          slackVarName = "*slack-{}*".format(self.__slackVars)
        elif moreThan:
          slackVarName = "*surplus-{}*".format(self.__slackVars)
        if lessThan or moreThan:
            varLen = len(self.__vars) + 1
            self.__vars[slackVarName] = varLen
            self.__S = np.insert(self.__S, self.__S.shape[1], 0, axis=1)
      row = [0] * (len(self.__vars) + 1)
      row[0] = bound if not objective else 0.0
      for x in coeffAndVar:
        row[self.__vars[x[1]]] = x[0]
        row[self.__vars["*neg*" + x[1]]] = -x[0]
      if not objective:
        if lessThan:
          row[self.__vars[slackVarName]] = 1.0
        if moreThan:
          row[self.__vars[slackVarName]] = -1.0
      #If row is objective, replace first row. Otherwise add new row
      if objective:
        self.__S[0,:] = np.array(row)
      else:
        self.__S = np.insert(self.__S, self.__S.shape[0], row, axis=0)
      if not objective and (lessThan or moreThan):
        self.__slackVars += 1

    '''
    Core execution of Simplex algorithm.
    Set twoPhase to True to execute first phase of two-phase Simplex
    Returns True if simplex loop halts succesfully, False if more iterations needed
      or LP is unbounded, in which case isDone is set to True.
    '''
    def __coreSimplex(self, twoPhase=False):
      divide = np.vectorize(lambda a,b: a/b if b > 10e-6 else np.inf)
      objRow = 1 if twoPhase else 0
      negRow0 = self.__S[objRow,1:] >= 0.0
      if np.all(negRow0):
        self.status = "optimum"
        return True
      #TODO: different strategies for choosing pivot column
      pivotColIdx = 1 + np.where(negRow0 == False)[1][0]
      #End TODO
      pivotCol = self.__S[:, pivotColIdx][:,0]
      #If all pivot column has non-positive entries, problem is unbounded
      if np.all(pivotCol[objRow+1:] <= 0):
        self.status = "unbounded"
        self.__state = SolverState.DONE
        self.isDone = True
        return False
      pivotRowIdx = objRow + 1 + np.argmin(divide(self.__S[objRow+1:,0], pivotCol[objRow+1:]))
      self.__basis[pivotRowIdx-objRow-1] = pivotColIdx-1
      self.__S[pivotRowIdx] = self.__S[pivotRowIdx] / pivotCol[pivotRowIdx]
      for i in range(self.__S.shape[0]):
        if i == pivotRowIdx:
          continue
        self.__S[i] -= pivotCol[i]*self.__S[pivotRowIdx]
      return False
    
    '''Clear SimplexSolver object'''
    def clear(self):
      self.__S = np.matrix([0.0])
      self.__vars = {}
      self.__invVars = []
      self.__slackVars = 0
      self.__varsToDelete = []
      self.__basis = np.array([])
      self.__missingBasis = 0
      self.__state = SolverState.PRESOLVE
      
      self.opt = np.inf
      self.optVars = {}
      self.status = "objective not set"
      self.isDone = False

    '''
    Define objective for LP program using a string.
    Statement should be separated to terms via plus or minus signs,
    each term should have floating number as coefficient on the left side
    and a variable name on the right side, joined by multiplication sign.
    Example: 2.0*x1 + 3.6e10*x2 - x3 + x4 - 7.8*x5
    '''
    def setObjective(self, str):
      self.__setRow(str, objective=True)
      self.status = "not started"

    '''
    Add a constraint to LP program using a string.
    Statement should be separated to terms via plus or minus signs,
    each term should have floating number as coefficient on the left side
    and a variable name on the right side, joined by multiplication sign.
    Constraint should end with ( <= | = | => ) and a single floating number.
    Example: 2.0*x1 + 3.6e10*x2 - x3 + x4 - 7.8*x5 <= 6.6
    '''
    def addConstraint(self, str):
      self.__setRow(str, objective=False)
    
    '''Before starting Simplex algorithm, subtract basis rows from the objective row'''
    def __subtractBasisFromObjectiveRow(self, twoPhase=False):
      objRow = 1 if twoPhase else 0
      for idx,col in enumerate(self.__basis):
        self.__S[objRow,:] = self.__S[objRow,:] - self.__S[objRow,col+1]*self.__S[objRow+idx+1,:]
    
    '''Calculation up to core simplex execution'''
    def __presolveStep(self):
      self.status = "preprocessing"

      #Flip signs of a row if rhs is negative
      negativeRows = 1 + np.where(self.__S[1:,0] < 0)[0]
      self.__S[negativeRows,:] = -1.0 * self.__S[negativeRows,:]

      #Delete unneeded variables and their corresponding columns
      self.__S = np.delete(self.__S, self.__varsToDelete, 1)
      #Invert vars dictionary
      self.__invVars = [""] * (len(self.__vars) + 1)
      for var, idx in self.__vars.items():
        self.__invVars[idx] = var
      self.__invVars = np.array(self.__invVars)
      self.__invVars = np.delete(self.__invVars, self.__varsToDelete, 0)

      #Simplex table without objective row and boundary column
      SnoZero = self.__S[1:,1:]
      colSize = SnoZero.shape[0]
      #Detect unit columns as basis
      unit = np.zeros(colSize)
      self.__basis = np.full(colSize, -1)
      for i in range(colSize):
        unit[i] = 1.0
        for j in range(SnoZero.shape[1]):
          if np.all(np.transpose(SnoZero[:,j]) == unit):
            self.__basis[i] = j
            break
        unit[i] = 0
      self.__missingBasis = np.count_nonzero(self.__basis < 0)
      #Couldn't find whole basis, do two phase Simplex
      if self.__missingBasis > 0:
        self.__S = np.insert(self.__S, 1, 0, axis=0)
        unit = np.zeros(self.__S.shape[0])
        unit[1] = 1
        j = 0
        for idx,col in enumerate(self.__basis):
          if col >= 0:
            self.__basis[idx] += self.__missingBasis
            continue
          unit[idx+2] = 1
          self.__S = np.insert(self.__S, j+1, unit, axis=1)
          self.__basis[idx] = j
          unit[idx+2] = 0
          j += 1
        self.__subtractBasisFromObjectiveRow(True)
        self.__state = SolverState.PHASE1_STEP
        self.status = "calculating phase 1 Simplex"
      else:
        self.__subtractBasisFromObjectiveRow()
        self.__state = SolverState.PHASE2_STEP
        self.status = "calculating"

    '''Phase 1 Simplex in two-phase Simplex execution'''
    def __phase1Step(self):
      if self.__coreSimplex(True):
        #Feasible set empty, end solving
        if self.__S[1,0] > 10e-6:
          self.__state = SolverState.DONE
          self.status = "feasible set empty"
          self.isDone = True
        else:
          self.__state = SolverState.DRIVEOUT_STEP
          self.status = "driving out artificial variables"

    '''
    Iteration of driving out artificial variables in case after
    phase 1 there are still artificial variables in basis
    '''
    def __driveoutStep(self):
      if np.all(self.__basis >= self.__missingBasis):
        self.__S = np.delete(self.__S, 1, 0)
        self.__S = np.delete(self.__S, range(1, 1+self.__missingBasis), 1)
        self.__basis -= self.__missingBasis
        self.__subtractBasisFromObjectiveRow()
        self.__state = SolverState.PHASE2_STEP
        self.status = "calculating phase 2 Simplex"
      else:
        #Artificial variable in basis, need to drive it out
        artificialsInBasis = 2 + np.where(self.__basis < self.__missingBasis)[0]
        for pivotRowIdx in artificialsInBasis:
          nonzeroEntires = np.where(self.__S[pivotRowIdx,(1+self.__missingBasis):] != 0)[1]
          if len(nonzeroEntires) == 0:
            continue
          for x in nonzeroEntires:
            if x not in self.__basis:
              pivotColIdx = 1 + self.__missingBasis + x
              break
          self.__basis[pivotRowIdx-2] = pivotColIdx-1
          self.__S[pivotRowIdx] = self.__S[pivotRowIdx] / self.__S[pivotRowIdx, pivotColIdx]
          for i in range(self.__S.shape[0]):
            if i == pivotRowIdx:
              continue
            self.__S[i] -= self.__S[i, pivotColIdx]*self.__S[pivotRowIdx]
          break

    '''Standard Simplex or phase 2 Simplex for two-phase Simplex execution. Gathers results on succesful execution'''
    def __phase2Step(self):
      if self.__coreSimplex():
        #Gather results
        resultVector = self.__S[:,0]
        self.opt = -resultVector[0,0]
        self.optVars = {}
        for idx, var in enumerate(self.__basis):
          self.optVars[self.__invVars[var+1]] = resultVector[idx+1,0]
        for var in self.__vars:
          if var not in self.optVars:
            self.optVars[var] = 0.0
        self.__state = SolverState.DONE
        self.status = "done"
        self.isDone = True

    '''Dummy method in case iteration called on solver that is done'''
    def __doneStep(self):
      pass

    def __init__(self):
      self.__S = np.matrix([0.0])         #Simplex table
      self.__vars = {}                    #Dictionary of variable names to table column index
      self.__invVars = []                 #Reverse of vars
      self.__slackVars = 0                #Counter of slack/surplus variables
      self.__varsToDelete = []            #Variables with infeasible constraints to delete
      self.__basis = np.array([])         #Basis of current solution as array of variable indices
      self.__missingBasis = 0             #Amount of rows not assigned to basis
      self.__state = SolverState.PRESOLVE #State to enable external loop
      self.__stateMethods = { \
          SolverState.PRESOLVE: self.__presolveStep, \
          SolverState.PHASE1_STEP: self.__phase1Step, \
          SolverState.DRIVEOUT_STEP: self.__driveoutStep, \
          SolverState.PHASE2_STEP: self.__phase2Step, \
          SolverState.DONE: self.__doneStep \
      }
      
      self.opt = np.inf                   #Optimum value after solving LP
      self.optVars = {}                   #Variable substitutions making up optimum value
      self.status = "objective not set"   #Human readable status
      self.isDone = False                 #Solver has ended its execution

    '''One step of solving process, could be used to analyze solving process from outside this class'''
    def solveStep(self):
      self.__stateMethods[self.__state]()

    '''Solve call'''
    def solve(self):
      while self.__state is not SolverState.DONE:
        self.solveStep()
