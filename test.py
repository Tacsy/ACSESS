
import filters
filters.Init()
bred=filters.ActiveFilters['Bredt\'s rule']

from rdkit import Chem
nor = Chem.MolFromSmiles('C\\1=C\\C2CC/1CC2')
nor1 = Chem.MolFromSmiles('C1CC2=CC1CC2')
nor2 = Chem.MolFromSmiles('C\\1=C\\C2CC/1=CC2')
print "nor"
print bred(nor)
print "nor1"
print bred(nor1)
print "nor2"
print bred(nor2)
