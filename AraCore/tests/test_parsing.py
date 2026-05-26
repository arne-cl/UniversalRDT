import run_rdt


def test_parse_rxn_header(sample_rxn_text):
    from_num, to_num = run_rdt.parse_rxn_header(sample_rxn_text)
    assert from_num == 2
    assert to_num == 1


def test_split_rxn_to_mols(sample_rxn_text):
    mol_blocks = run_rdt.split_rxn_to_mols(sample_rxn_text)
    assert len(mol_blocks) == 3
    for block in mol_blocks:
        assert block.startswith("M0")


def test_parse_mdl_atom_table(sample_mol_01_mdl):
    atoms = run_rdt.parse_mdl_atom_table(sample_mol_01_mdl)
    expected = [
        ("O", 1), ("C", 2), ("C", 3), ("O", 4), ("C", 5),
        ("O", 6), ("P", 7), ("O", 8), ("O", 9), ("O", 10),
    ]
    assert atoms == expected


def test_parse_mdl_atom_table_mol02(sample_rxn_dir):
    mdl = (sample_rxn_dir / "MOL_02.mdl").read_text()
    atoms = run_rdt.parse_mdl_atom_table(mdl)
    expected = [
        ("O", 11), ("C", 12), ("C", 13), ("O", 14), ("C", 15),
        ("O", 16), ("P", 17), ("O", 18), ("O", 19), ("O", 20),
    ]
    assert atoms == expected


def test_parse_mdl_atom_table_mol03(sample_rxn_dir):
    mdl = (sample_rxn_dir / "MOL_03.mdl").read_text()
    atoms = run_rdt.parse_mdl_atom_table(mdl)
    assert len(atoms) == 20
    assert atoms[0] == ("O", 1)
    assert atoms[1] == ("C", 2)


def test_parse_inchi_atom_order(sample_mol_01_inchi):
    order = run_rdt.parse_inchi_atom_order(sample_mol_01_inchi)
    assert order == [2, 5, 3, 1, 4, 8, 9, 10, 6, 7]


def test_parse_inchi_atom_order_mol02(sample_mol_02_inchi):
    order = run_rdt.parse_inchi_atom_order(sample_mol_02_inchi)
    assert order == [3, 5, 2, 4, 1, 8, 9, 10, 6, 7]


def test_parse_inchi_atom_order_mol03(sample_mol_03_inchi):
    order = run_rdt.parse_inchi_atom_order(sample_mol_03_inchi)
    assert order == [12, 8, 13, 4, 2, 15, 3, 1, 16, 9, 10, 18, 5, 6, 20, 11, 7, 14, 17, 19]


def test_parse_inchi_atom_order_no_aux():
    inchi_text = "InChI=1S/H2O/h1H2\n"
    order = run_rdt.parse_inchi_atom_order(inchi_text)
    assert order == [1]
