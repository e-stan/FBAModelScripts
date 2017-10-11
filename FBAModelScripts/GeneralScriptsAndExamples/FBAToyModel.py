from cobra.core import *
import cobra
import cobra.io

m1 = Metabolite('A','A')
m2 = Metabolite('B','B')
m3 = Metabolite('C','C')
m4 = Metabolite('AB','AB')
m5 = Metabolite('AC','AC')
m6 = Metabolite('ABC','ABC')

"""
->A
->B
A + B -> AB
AB -> A + B
->AB

A + C -> AC
AB + C <=> ABC
AC + B <=> ABC
A + B + C -> ABC

->C
ABC->


"""

r1 = Reaction('r1')
r1.lower_bound = 0
r1.upper_bound = 10
r1.add_metabolites({m1:-1,m2:-1,m4:1})
r1.gene_reaction_rule = 'A'

r2 = Reaction('r2')
r2.lower_bound = 2
r2.upper_bound = 10
r2.add_metabolites({m1:1,m2:1,m4:-1})
r2.gene_reaction_rule = 'A'

"""

r2 = Reaction('r2')
r2.lower_bound = 0
r2.upper_bound = 1000
r2.add_metabolites({m1:-1,m3:-1,m5:1})
r2.gene_reaction_rule = 'B'

r3 = Reaction('r3')
r3.lower_bound = 0
r3.upper_bound = 1000
r3.add_metabolites({m1:-1,m3:-1,m6:1})
r3.gene_reaction_rule = 'C'


r4 = Reaction('r4')
r4.lower_bound = 0
r4.upper_bound = 1000
r4.add_metabolites({m1:1,m3:1,m6:-1})
r4.gene_reaction_rule = 'D'


r5 = Reaction('r5')
r5.lower_bound = 0
r5.upper_bound = 1000
r5.add_metabolites({m5:-1,m2:-1,m6:1})
r5.gene_reaction_rule = 'E'


r6 = Reaction('r6')
r6.lower_bound = 0
r6.upper_bound = 1000
r6.add_metabolites({m5:1,m2:1,m6:-1})
r6.gene_reaction_rule = 'F'


r7 = Reaction('r7')
r7.lower_bound = 0
r7.upper_bound = 1000
r7.add_metabolites({m1:-1,m2:-1,m3:-1,m6:1})
r7.gene_reaction_rule = 'G'
"""
r11 = Reaction('r11')
r11.lower_bound = 0
r11.upper_bound = 1000
r11.add_metabolites({m4:-1})
r11.gene_reaction_rule = 'G'

r8 = Reaction('r8')
r8.lower_bound = 5
r8.upper_bound = 5
r8.add_metabolites({m1:1})
r8.gene_reaction_rule = 'G'

r9 = Reaction('r9')
r9.lower_bound =5
r9.upper_bound = 5
r9.add_metabolites({m2:1})
r9.gene_reaction_rule = 'G'
"""
r10 = Reaction('r10')
r10.lower_bound = 0
r10.upper_bound = 1000
r10.add_metabolites({m3:1})
r10.gene_reaction_rule = 'G'
"""
myModel = cobra.Model('toyModel')
#myModel.add_reactions([r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11])
myModel.add_reactions([r1,r2,r11,r8,r9])
myModel.objective = 'r11'

cobra.io.write_sbml_model(myModel,'toyModel.xml',use_fbc_package=False)

