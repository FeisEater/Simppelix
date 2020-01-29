import simplexsolver
import time

solver = simplexsolver.SimplexSolver()          #Initiate solver object
solver.setObjective("x1 + x2 + x3 + x4 + x5")   #Set objective
solver.addConstraint("3*x1 + 2*x2 + x3 = 1")    #Add constraints
solver.addConstraint("2*x1 - x2 + x4 = 2")
solver.addConstraint("3*x2 - x1 + x5 = 3")
solver.addConstraint("x1 >= 0")                 #Remember to specify non-negative variables
solver.addConstraint("x2 >= 0")
solver.addConstraint("x3 >= 0")
solver.addConstraint("x4 >= 0")
solver.addConstraint("x5 >= 0")
solver.solve()                                  #Call to solve
print("Optimum: {}".format(solver.opt))         #Print optimum value

solver.clear()                                  #Clear solver object to reuse it
solver.setObjective("-x2")
solver.addConstraint("3*x1+2*x2<=6")
solver.addConstraint("-3*x1+2*x2<=0")
solver.addConstraint("x1>=0")
solver.addConstraint("x2>=0")
#Solve iteration can be called by hand
while not solver.isDone:                        #isDone marks when solver shouldn't be iterated anymore
  print(solver.status)                          #status is human readable sentence of solver's current status
  solver.solveStep()                            #One iteration of solver algorithm

print("Optimum: {}".format(solver.opt))
print("with following variable substitutions:")
for var, value in solver.optVars.items():       #optVars is dictionary of variables as keys and their substitutions as values
  print("{}: {}".format(var, value))

solver = simplexsolver.SimplexSolver()
solver.setObjective("x + 4*y + 9*z")
solver.addConstraint("x + y <= 5")
solver.addConstraint("x + z >= 10")
solver.addConstraint("-y + z = 7")
solver.addConstraint("x <= 4")
solver.addConstraint("y >= -1")
solver.addConstraint("y <= 1")
solver.solve()
print("Optimum: {}".format(solver.opt))
for var, value in solver.optVars.items():
  print("{}: {}".format(var, value))

solver.clear()
solver.readMps("../mps/afiro.mps")              #Read MPS file
solver.solve()
print("Optimum: {}".format(solver.opt))

