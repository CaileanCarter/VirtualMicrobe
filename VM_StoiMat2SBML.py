from libsbml import *


MicrobeContents = {'LPS': 25, 'PG': 25, 'Glucose': 100, 'NAD': 100, 'NADH': 100, 'ADP': 48, 'ATP': 48, 'nitrate': 10, 'amino_acids': 60, 'RNA': 60, 'precursor': 10, 'nucleotides': 50, 
                    'Hydrogen': 40, 'RNAstrand': 0, 'peptides': 0, 'DNA': 0, 'e': 0, 'PMF': 50, 'Genome': 100, 'CellWallIntegrity': 100, 'Oxygen': 100, 'nitrite' : 0}


StoiMat = {'Substrate': ['LPS', 'PG', 'Glucose', 'NAD', 'NADH', 'ADP', 'ATP', 'nitrate', 'amino_acids', 'RNA', 'precursor', 'e', 'Hydrogen', 'Oxygen', 'nitrite', 'PMF', 'DNA', 'RNAstrand', 'peptides', 'CellWallIntegrity', 'nucleotides', 'Genome'], 
'DNARepair': {'ADP': 4, 'ATP': -4, 'precursor': -4, 'Genome': 1}, 
'WallSynthesis': {'LPS': -1, 'PG': -1, 'ADP': 2, 'ATP': -2, 'CellWallIntegrity': 1},
'WallDecay': {'LPS': 1, 'PG': 1, 'CellWallIntegrity': -1},
'AAsynthesis': {'ADP': 1, 'ATP': -1, 'amino_acids': 1, 'precursor': -1},
'RnaSynthesis': {'ADP': 1, 'ATP': -1, 'RNA': 1, 'precursor': -1},
'PrecursorSynthesis': {'Glucose': -1, 'ADP': 1, 'ATP': -1, 'precursor': 1}, 
'WallComponents': {'LPS': 1, 'PG': 1, 'ADP': 1, 'ATP': -1, 'precursor': -2},
'ATPhydrolysis': {'ADP': 1, 'ATP': -1},
'Gluconeogenesis': {'Glucose': 1, 'ADP': 1, 'ATP': -1, 'precursor': -1},
'ADPsynthesis': {'ADP': 1, 'precursor': -1},
'NucleotideSynthesis': {'ADP': 1, 'ATP': -1, 'precursor': -1, 'nucleotides': 1},
'ProteinSynthesis': {'ADP': 20, 'ATP': -20, 'amino_acids': -20, 'RNAstrand': -1, 'peptides': 20},
'RNAsynthesis': {'ADP': 20, 'ATP': -20, 'RNA': -20, 'RNAstrand': 1},
'DNAsynthesis': {'ADP': 1, 'ATP': -1, 'DNA': 1, 'peptides': -1, 'nucleotides': -1},
'Replicate': {'LPS': -10, 'PG': -10, 'ADP': 10, 'ATP': -10, 'DNA': -10},
'ElectronTransportChain': {'NAD': 1, 'NADH': -1, 'PMF': 1},
'Gly_TCA': {'Glucose': -1, 'NAD': -12, 'NADH': 12, 'Hydrogen': -12},
'OxidativePhosphorylation': {'ADP': -1, 'ATP': 1, 'Hydrogen': 2, 'PMF': -1}} #Stoichiometry Matrix



document = SBMLDocument(3, 1)

model = document.createModel()

model.setTimeUnits("second")
model.setExtentUnits("mole")
model.setSubstanceUnits("mole")


per_second = model.createUnitDefinition()
per_second.setId('per_second')
unit = per_second.createUnit()
unit.setKind(UNIT_KIND_SECOND)
unit.setExponent(-1)
unit.setScale(0)
unit.setMultiplier(1)



metabolites = model.createCompartment()
metabolites.setId("metabolites")
metabolites.setConstant(True)
metabolites.setSize(len(StoiMat.keys()))
metabolites.setSpatialDimensions(3)
metabolites.setUnits('micro litre')



for metabolite in StoiMat["Substrate"]:

    exec(metabolite + "= model.createSpecies()")
    exec(metabolite + ".setId(metabolite)")
    exec(metabolite + ".setConstant(False)")
    exec(metabolite + ".setInitialAmount(MicrobeContents[metabolite])")
    exec(metabolite + ".setSubstanceUnits('mole')")
    exec(metabolite + ".setBoundaryCondition(False)")
    exec(metabolite + ".setHasOnlySubstanceUnits(False)")


refNum = 0

for reaction in StoiMat.keys():
    if reaction != "Substrate":
        exec(reaction + " = model.createReaction()")
        exec(reaction + ".setId(reaction)")
        exec(reaction + ".setFast(False)")

        for reaction_substrate_product in StoiMat[reaction].keys():
            
            ref = reaction_substrate_product + "_" + str(refNum)

            if StoiMat[reaction][reaction_substrate_product] < 0:
                exec(ref + "= " + reaction + ".createReactant()")
                exec(ref + ".setSpecies(reaction_substrate_product)")
                exec(ref + ".setConstant(True)")
            
            else:
                exec(ref + "= " + reaction + ".createProduct()")
                exec(ref + ".setSpecies(reaction_substrate_product)")
                exec(ref + ".setConstant(True)")
            
            refNum += 1

print(writeSBMLToString(document))


