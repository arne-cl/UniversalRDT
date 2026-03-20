# UniversalRDT
initial state

Short overview of the workflow and its principles:

1) You have a model/set of reactions
2) You need SMILES notation of all metabolites present in the reactions
3) You calculate INCHi-keys for all metabolites
4) You create SMILES strings for all reactions, each in iits own directory, including the model-IDs and InCHi-keys for the used metabolites
5) You apply RDT (could be replaced by another automatic atom mapping tool) on the reaction
6) You get a .rxn file as result
7) you parse the individual molecules out that file
8) you calculate the InCHi-key for all molecules
9) you identify the model ID for the molecule sequence of the .rxn file
10) you extract the mapping information
11) you calculate InCHi for all individual molecules of the .rxn file, with auxillary information (containing a mapping of which InCHi-ordered atom is equivalent to which input atom
12) You adapt the mapping information, so now all mappings identify atoms with their InCHi-order in the molecule
13) combine the InCHi-ordered mappings from all reactions
