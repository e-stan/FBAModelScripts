from cobra.core import *
import cobra
import cobra.io


m1 = Metabolite('A','A')
m2 = Metabolite('B','B')
m4 = Metabolite('AB','AB')

#A + B -> AB
r1 = Reaction('r1')
r1.lower_bound = 0
r1.upper_bound = 10
r1.add_metabolites({m1:-1,m2:-1,m4:1})
r1.gene_reaction_rule = 'A'

#AB -> A + B (Remove to see irreversible version of network)
r2 = Reaction('r2')
r2.lower_bound = 2
r2.upper_bound = 10
r2.add_metabolites({m1:1,m2:1,m4:-1})
r2.gene_reaction_rule = 'A'

#->AB
r11 = Reaction('r11')
r11.lower_bound = 0
r11.upper_bound = 1000
r11.add_metabolites({m4:-1})
r11.gene_reaction_rule = 'G'

#->A
r8 = Reaction('r8')
r8.lower_bound = 5
r8.upper_bound = 5
r8.add_metabolites({m1:1})
r8.gene_reaction_rule = 'G'

#->B
r9 = Reaction('r9')
r9.lower_bound =5
r9.upper_bound = 5
r9.add_metabolites({m2:1})
r9.gene_reaction_rule = 'G'

myModel = cobra.Model('toyModel')
myModel.add_reactions([r1,r2,r11,r8,r9])
myModel.objective = 'r11'

cobra.io.write_sbml_model(myModel,'toyModel.xml',use_fbc_package=False)

