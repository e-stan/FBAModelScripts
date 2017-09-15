""""
COBRA model of
E + S <-> ES -> E + P
-> E
-> S
S ->
P->
"""

from cobra.core import *
import cobra.io


m1 = Metabolite('E','E')
m2 = Metabolite('S','S')
m3 = Metabolite('ES','ES')
m4 = Metabolite('P','P')

#E + S -> ES
r1 = Reaction('r1')
r1.lower_bound = 0
r1.upper_bound = 1000
r1.add_metabolites({m1:-1,m2:-1,m3:1})
r1.gene_reaction_rule = 'EC'
#ES -> E + S
r2 = Reaction('r2')
r2.lower_bound = 2
r2.upper_bound = 1000
r2.add_metabolites({m1:1,m2:1,m3:-1})
r2.gene_reaction_rule = 'EC'
#ES -> E + P
r3 = Reaction('r3')
r3.lower_bound = 0
r3.upper_bound = 1000
r3.add_metabolites({m3:-1,m1:1,m4:1})
r3.gene_reaction_rule = 'P'

#->E
r4 = Reaction('r4')
r4.lower_bound = 0
r4.upper_bound = 5
r4.add_metabolites({m1:1})
r4.gene_reaction_rule = 'E'

#->S
r5 = Reaction('r5')
r5.lower_bound = 5
r5.upper_bound = 5
r5.add_metabolites({m2:1})
r5.gene_reaction_rule = 'S'

#S->
r6 = Reaction('r6')
r6.lower_bound = 0
r6.upper_bound = 1000
r6.add_metabolites({m2:-1})
r6.gene_reaction_rule = 'S'

#P->
r7 = Reaction('r7')
r7.lower_bound = 0
r7.upper_bound = 1000
r7.add_metabolites({m4:-1})
r7.gene_reaction_rule = 'PD'





myModel = cobra.Model('toyModel')
myModel.add_reactions([r1,r2,r3,r4,r5,r6,r7])
myModel.objective = 'r7'

cobra.io.write_sbml_model(myModel,'toyModel.xml',use_fbc_package=False)

