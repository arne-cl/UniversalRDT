import run_rdt


def test_build_mapping_lines_from_side():
    rdt_index = [
        ("O", 1), ("C", 2), ("C", 3), ("O", 4), ("C", 5),
        ("O", 6), ("P", 7), ("O", 8), ("O", 9), ("O", 10),
    ]
    inchi_order = [2, 5, 3, 1, 4, 8, 9, 10, 6, 7]
    lines = run_rdt.build_mapping_lines(
        rdt_index, inchi_order, "M_GAP[h]", "from"
    )
    assert len(lines) == 10
    assert lines[0] == "2\tfrom\tM_GAP[h]:C#1="
    assert lines[1] == "5\tfrom\tM_GAP[h]:C#2="
    assert lines[4] == "4\tfrom\tM_GAP[h]:O#2="
    assert lines[9] == "7\tfrom\tM_GAP[h]:P#1="


def test_build_mapping_lines_to_side():
    rdt_index = [
        ("O", 1), ("C", 2), ("O", 14), ("C", 13), ("O", 10),
        ("O", 9), ("O", 6), ("C", 5), ("O", 20), ("O", 19),
        ("O", 16), ("C", 15), ("C", 12), ("O", 11), ("C", 3),
        ("O", 4), ("P", 17), ("O", 18), ("P", 7), ("O", 8),
    ]
    inchi_order = [12, 8, 13, 4, 2, 15, 3, 1, 16, 9, 10, 18, 5, 6, 20, 11, 7, 14, 17, 19]
    lines = run_rdt.build_mapping_lines(
        rdt_index, inchi_order, "M_FBP[h]", "to"
    )
    assert len(lines) == 20
    assert lines[0] == "15\tto\tM_FBP[h]:C#1,"
    assert lines[-1] == "7\tto\tM_FBP[h]:P#2,"


def test_assemble_mapping(sample_mapping_txt, sample_mapping_lines):
    result = run_rdt.assemble_mapping(sample_mapping_lines)
    assert result == sample_mapping_txt.strip()


def test_assemble_mapping_removes_hydrogen():
    lines = "1\tfrom\tM_X[h]:C#1=\n2\tfrom\tM_X[h]:H#1=\n3\tto\tM_Y[h]:C#1,\n"
    result = run_rdt.assemble_mapping(lines)
    assert ":H#" not in result
    assert "M_X[h]:C#1=M_Y[h]:C#1" == result
