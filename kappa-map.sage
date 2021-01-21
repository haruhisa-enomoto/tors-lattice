def max_element(poset,subset):
    r"""
    Return one maximal element of a subset of a poset.
    It returns a maximal element which is found first.

    INPUT:

    - ``poset`` -- a poset,
      which should be an instance of FinitePoset_with_category

    - ``subset`` -- a subset of ``poset``,
      which can be a set, a tuple, a list, and so on

    OUTPUT:

    an element in ``subset``

    """
    subset = list(subset)
    maximum = subset[0]
    for p in subset:
        if poset.is_greater_than(p,maximum):
            maximum = p
    return maximum

def kappa_map(lattice, j):
    r"""
    Return $\kappa(j)$ for a join-irreducible element $j$
    in a finite semidistributive lattice,
    see [BTZ] or [FJN].

    REFERENCES:

    - [BTZ] E. Barnard, G. Todorov, S. Zhu, Dynamical Combinatorics and Torsion Classes, arXiv:1911.10712.

    - [FJN] R. Freese, J. Jezek, J. B. Nation, Free lattices, volume 42 of Mathematical Surveys and Monographs. American Mathematical Society, Providence, RI, 1995.
    """

    # Ensures that a lattice is facade
    lattice = LatticePoset(lattice, facade=True)
    j_irr = lattice.join_irreducibles()
    m_irr = lattice.meet_irreducibles()

    if  j not in j_irr:
        raise ValueError(f"{j} is not join-irreducible.")

    # j_star is the unique element covered by j.
    j_star = lattice.lower_covers(j)[0]

    # kappa(j) is the unique maximum element x satisfying
    # j_star \leq x and j \not\leq x.
    tmp_subset = [
        x for x in lattice
        if lattice.is_lequal(j_star,x)
        and not lattice.is_lequal(j,x)
        ]
    m = max_element(lattice,tmp_subset)
    assert m in m_irr, f"BUG: {m} is not meet-irreducible."
    return max_element(lattice,tmp_subset)

def extended_kappa_map(lattice,x):

    # First, compute the canonical join representation
    # of x.
