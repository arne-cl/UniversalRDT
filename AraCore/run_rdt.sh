#!/bin/bash
for rxn_folder in reaction_intermediates/*
do
echo $rxn_folder
cd $rxn_folder
# first: run rdt
smiles=$( cat rxn.smiles )
java -jar ${RDT_JAR:-../../../rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar} -Q SMI -q "$smiles" -g -c -b -j AAM -f TEXT
# now split RDT output into mol-files
rm MOL_*
csplit -f MOL_ ECBLAST_smiles_AAM.rxn '/$MOL/' {*}
# handle the stuff before the first $MOL
# it contains in the last line the number of from-side and to-side molecules
from_to=$(grep -m1 -B1 '$MOL' ECBLAST_smiles_AAM.rxn | head -n 1)
from_num=$(echo "$from_to" | head -c3 | sed 's/ //g')
to_num=$(echo "$from_to" | tail -c+4 | sed 's/ //g')
rm MOL_00

rm mapping_lines.txt
counter=1

for fn2 in MOL_??
do
tail -n +2 $fn2 > ${fn2}.mdl
# MOL_nn.mdl now contains one molecule of the rdt .rxn file
# NOTE: We now preserve M CHG lines to maintain charge information for InChIKey matching
# This fixes the InChIKey mismatch bug where charge-stripped keys didn't match references
awk '(NF==16) { print $4 "\t" $14; }; (NF==15) { print $4 "\t" (0+$13); }' ${fn2}.mdl > ${fn2}.rdt_index
#-xT/nochg - RDT sometimes changes protons, leave them out of InChI (for canonical atom ordering)
#-xa we want aux info - it contains the original position in the input - the base for our mapping
obabel -i mdl  ${fn2}.mdl -o inchi -xa -xT/nochg -O ${fn2}.inchi
obabel -i mdl  ${fn2}.mdl -oinchikey -O ${fn2}.inchikey

if [ -s ${fn2}.inchikey ]
then
# Lookup species by InChIKey in the reference table
species_id_without_cmp=$(grep "$(cat ${fn2}.inchikey)" species_id_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')

# Fallback: if exact match fails, try first 14 characters (connectivity layer only)
# This handles cases where InChIKeys differ due to stereochemistry or charge handling
if [ -z "$species_id_without_cmp" ] && [ -s "${fn2}.inchikey" ]; then
    species_id_without_cmp=$(grep "$(head -c14 ${fn2}.inchikey)" species_id_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
fi

if [ -n "$species_id_without_cmp" ]
then
if [ $counter -le $from_num ]
then
species_id=$(grep "$species_id_without_cmp" from_species_with_cmp)
mapping_side="from"
mapping_end="="
else
species_id=$(grep "$species_id_without_cmp" to_species_with_cmp)
mapping_side="to"
mapping_end=","
fi
echo $species_id > ${fn2}.species_id

atom_counter=0
last_element=
#from each line in InChi-Atom order we want the element and the mapping-index
#inchi-order is element-wise, so first all C, then all N, etc.
#so we use a format: <species>:<element>#<index in InChI for this atom per element - i.e. starts with 1 for every element>
if grep -q 'AuxInfo.*/N:' ${fn2}.inchi
then
	inchi_index=$(grep 'AuxInfo' ${fn2}.inchi | sed 's/^.*\/N://; s/\/.*$//; s/,/ /g')
else
	inchi_index=1
fi

#echo $inchi_index
for rdt_line in $inchi_index
do 
element_and_index=$(tail -n +$rdt_line ${fn2}.rdt_index | head -n 1)
element=$(echo "$element_and_index" | cut -f1)
mapping_index=$(echo "$element_and_index" | cut -f2)
#echo $element
#echo $mapping_index

if [ "$last_element" = "$element" ]
then
atom_counter=$(($atom_counter + 1))
else
last_element=$element
atom_counter=1
fi
echo ${mapping_index}"	${mapping_side}	"${species_id}":"${element}"#"${atom_counter}${mapping_end} >> mapping_lines.txt
done

fi

fi #if [ -s ${fn2}.inchikey ]
counter=$(($counter + 1))
done

#now create the actual mapping, remove Hydrogen atoms (which likely will not be matched)
sort -n  mapping_lines.txt | grep -v ':H#' | cut -f3 | tr '\n' ' ' | sed 's/ //g; s/,$//' > mapping.txt
cd -
done
