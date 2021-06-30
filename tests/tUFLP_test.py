"""Tests the tUFLP procedures. (:py:mod:`tUFLP`)"""
from tUFLP import *

@pytest.mark.parametrize("test_inst",
                         [generate_test_instance(15) for _ in range(100)] +
                         [generate_string_instance(15) for _ in range(100)] +
                         [generate_d4_instance(15) for _ in range(100)] +
                         [generate_organic_instance(15) for _ in range(100)])
def test_random_UFL(test_inst):
    """Tests the formulation for typed-UFLP (overlap DD) -- random instance."""
    TOL=1e-3
    S, f, fc, kb = test_inst
    print(f"Running a test with:\nS={S}; f={f}; fc={fc}; kb={kb}")

    m_naive = solve_with_MIP(S, f, fc, kb)
    m = solve_with_BDD_MIP(S, f, fc, kb)
    assert abs(m_naive.objVal - m.objVal) < TOL, f"Naive: {m_naive.objVal} (status {m_naive.status}), while BDD: {m.objVal} (status {m.status})"


@pytest.mark.parametrize("test_inst", [generate_test_instance(10)
                                       for _ in range(200)])
def test_triple(test_inst):
    """Tests the formulation for typed-UFLP."""
    TOL=1e-5
    S, f, fc, kb = test_inst

    m_naive = solve_with_MIP(S, f, fc, kb)
    m_DD_MIP = solve_with_BDD_MIP(S, f, fc, kb)
    m_NF = solve_with_align_BDD(S, f, fc, kb)
    print(f"Naive model: status={m_naive.status}, obj={m_naive.objVal}")
    print(f"BDD-MIP model: status={m_DD_MIP.status}, obj={m_DD_MIP.objVal}")
    print(f"BDD-NF model: status={m_NF.status}, obj={m_NF.objVal}")
    assert abs(m_naive.objVal - m_DD_MIP.objVal)<TOL, f"Naive: {m_naive.objVal} (status {m_naive.status}), while BDD: {m_DD_MIP.objVal} (status {m_DD_MIP.status})"
    assert abs(m_naive.objVal - m_NF.objVal)<TOL, f"Naive: {m_naive.objVal} (status {m_naive.status}), while BDD: {m_NF.objVal} (status {m_NF.status})"


@pytest.mark.parametrize("test_inst", [make_instance(7)
                                       for _ in range(100)])
def test_TypeSorter(test_inst):
    """Tests that TypeSorter finds an optimal order (compares to bruteforce)."""
    f_types, target_order = test_inst
    print(f"The instance is: f_types={f_types}, target={target_order}")
    assert get_score(find_correct_order(f_types, target_order), target_order) == get_score(
        bruteforce_correct_order(f_types, target_order), target_order)


@pytest.mark.parametrize("test_inst", [generate_test_instance(7)
                                       for _ in range(100)])
def test_randomized_cover_DD(test_inst):
    S, f, fc, kb = test_inst
    C, _ = build_cover_DD(S, f)
    Cr, _ = build_randomized_cover_DD(S, f)

    assert C.is_equivalent(Cr)[0], f"Not equivalent:\nS={S}, f={f}"


@pytest.mark.parametrize("test_inst", [generate_test_instance(8)
                                       for _ in range(100)])
def test_randomized_type_DD(test_inst):
    S, f, fc, kb = test_inst
    C, _ = build_type_DD(f, fc, kb)
    Cr, _ = build_randomized_type_DD(f, fc, kb)

    assert C.is_equivalent(Cr)[0], f"Not equivalent:\nS={S}, f={f}"
