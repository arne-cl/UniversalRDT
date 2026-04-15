#!/bin/bash
for rxn_id in $(cat rxns_to_correct.txt)
do
rxn_folder=reaction_intermediates/$rxn_id
cd $rxn_folder
# first: run rdt
smiles=$( cat rxn.smiles )
java -jar ${RDT_JAR:-../../rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar} -Q SMI -q "$smiles" -g -c -b -j AAM -f TEXT
# now split RDT output into mol-files
rm MOL_*
csplit -f MOL_ ECBLAST_smiles_AAM.rxn '/$MOL/' {*}
from_to=$(grep -m1 -B1 '$MOL' ECBLAST_smiles_AAM.rxn | head -n 1)
# handle the stuff before the first $MOL
# it contains in the last line the number of from-side and to-side molecules
from_num=$(echo "$from_to" | head -c3 | sed 's/ //g')
to_num=$(echo "$from_to" | tail -c+4 | sed 's/ //g')
rm MOL_00

echo -n "" > mapping_lines.txt
counter=1

for fn2 in MOL_??
do
tail -n +2 $fn2 | grep -v '^M  CHG' > ${fn2}.mdl
# MOL_nn.mdl now contains one molecule of the rdt .rxn file
awk '(NF==16) { print $4 "\t" $14; }; (NF==15) { print $4 "\t" (0+$13); }' ${fn2}.mdl > ${fn2}.rdt_index
#~/Documents/Promotion/INCHI-1-SRC/INCHI_EXE/bin/Linux/inchi-1 ${fn2}.mdl
#tail -n +3 $fn2 | grep -v '^M  CHG' |  sed '3,3 s/^\(...........\)0/\1 /'> ${fn2}.mdl
#-xT/nochg - RDT sometimes changes protons, leave them out of InChI
#-xa we want aux info - it contains the oririgal position in the input - the base for our mapping
obabel -i mdl  ${fn2}.mdl -o inchi -xa -xT/nochg -O ${fn2}.inchi
obabel -i mdl  ${fn2}.mdl -oinchikey -O ${fn2}.inchikey

if [ -s ${fn2}.inchikey ]
then
#now we need the name of the molecule metabolite
#species_id=$(grep $(cat ${fn2}.inchikey) ../../species_id_without_cmp_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
species_id=$(grep $(cat ${fn2}.inchikey) species_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
if [ -z "$species_id" ]
then
species_id=$(grep $(head -c14 ${fn2}.inchikey) species_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
fi
#grep $(cat ${fn2}.inchikey) ../../species_id_without_cmp_inchikey.txt
#cat ${fn2}.inchikey
#echo $species_id_without_cmp
#echo $counter $from_num
if [ $counter -le $from_num ]
then
# if we have the same inchikey from and to, then this will select the from species
species_id=$(echo $species_id | sed 's/ .*//')
mapping_side="from"
mapping_end="="
else
# if we have the same inchikey from and to, then this will select the to species
species_id=$(echo $species_id | sed 's/.* //')
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

fi #if [ -s ${fn2}.inchikey ]
counter=$(($counter + 1))
done

#now create the actual mapping, remove Hydrogen atoms (which likely will not be matched)
sort -n  mapping_lines.txt | grep -v ':H#' | cut -f3 | tr '\n' ' ' | sed 's/ //g; s/,$//' > mapping.txt
cd -
done
