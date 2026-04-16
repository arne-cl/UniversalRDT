#!/bin/bash
# Test to demonstrate InChIKey mismatch bug
#
# Root cause: The script removes M  CHG lines from MDL files,
# causing InChIKeys to differ from reference (which includes charge info)
#
# Reference InChIKeys: generated from SMILES with charge info → ends in -M
# Generated InChIKeys: from MDL without M CHG → ends in -N
# Result: grep fails, species_id_without_cmp becomes empty, script hangs

set -e

echo "=== InChIKey Mismatch Test ==="
echo ""

# Test molecule: orotate (charged)
CHARGED_SMILES="O=C([O-])c1cc(=O)[nH]c(=O)[nH]1"
NEUTRAL_SMILES="O=C(O)c1cc(=O)[nH]c(=O)[nH]1"

echo "1. Reference method (from recreate_data.sh):"
echo "   obabel -:'$CHARGED_SMILES' -oinchikey"
REF_KEY=$(obabel -:"$CHARGED_SMILES" -oinchikey 2>/dev/null)
echo "   Result: $REF_KEY"
echo ""

echo "2. Charged SMILES with -xT/nochg flag:"
echo "   obabel -:'$CHARGED_SMILES' -oinchikey -xT/nochg"
NOCHG_KEY=$(obabel -:"$CHARGED_SMILES" -oinchikey -xT/nochg 2>/dev/null)
echo "   Result: $NOCHG_KEY"
echo ""

echo "3. Neutral SMILES (simulating charge-stripped MDL):"
echo "   obabel -:'$NEUTRAL_SMILES' -oinchikey"
NEUTRAL_KEY=$(obabel -:"$NEUTRAL_SMILES" -oinchikey 2>/dev/null)
echo "   Result: $NEUTRAL_KEY"
echo ""

echo "=== Analysis ==="
echo "Reference file has: $REF_KEY"
echo "Generated has:      $NEUTRAL_KEY"
echo "Match: $([[ "$REF_KEY" == "$NEUTRAL_KEY" ]] && echo "YES" || echo "NO - BUG!")"
echo ""

if [[ "$REF_KEY" != "$NEUTRAL_KEY" ]]; then
    echo "*** BUG CONFIRMED ***"
    echo "The grep for InChIKey will fail because:"
    echo "  - species_id_inchikey.txt contains: $REF_KEY"
    echo "  - Generated from charge-stripped MDL: $NEUTRAL_KEY"
    echo ""
    echo "This causes species_id_without_cmp to be empty,"
    echo "which causes grep to hang when searching for empty pattern."
    exit 1
else
    echo "No bug detected - InChIKeys match!"
    exit 0
fi
