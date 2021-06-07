"""
Testing module for `BDD.py`

---
(c) A. Bochkarev, Clemson University, 2021
abochka@g.clemson.edu
"""
import BDD
import pytest
import numpy as np

@pytest.mark.parametrize("D", [BDD.BDD.random(N = 1 + np.random.randint(1, 10),
                                              p = np.random.uniform(),
                                              weighted=np.random.choice([True, False]))
                               for _ in range(5000)])
def test_load_save(D):
    """Tests that load/save functionality works as intended."""
    D.save("tests/BDD_load_save.bdd")
    print(f"D is: {D}")
    D2 = BDD.BDD()
    D2.load("tests/BDD_load_save.bdd")
    assert D.profile() == D2.profile()
