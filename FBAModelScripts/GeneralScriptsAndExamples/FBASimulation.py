import cobra.io

myModel = cobra.io.read_sbml_model('toyModel.xml')
myModel.solver = 'gurobi'

sol = myModel.optimize()

fvaSol = cobra.flux_analysis.flux_variability_analysis(myModel,loopless=True)
print fvaSol
#print sol.fluxes
for x in myModel.reactions:
    print str(x.id) + ' : ' + str(x.reaction) + ' : '  + str(sol[x.id])

