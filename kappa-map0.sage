def kappa(self, j, check = False):
    r"""
    Return `\kappa(j)` for a join-irreducible element `j`
    in a finite lattice ``self`` if it exists.

    INPUT:

    - ``self`` -- an instance of
      :class:`sage.combinat.posets.lattices.FiniteLatticePoset`

    - ``j`` -- an element of ``self``,
      which is expected to be join-irreducible

    - ``check`` -- a Boolean (default: ``False``),
      whether to check that ``j`` is indeed join-irreducible

    OUTPUT:

    an element of ``self``, or None if it does not exist.

    .. SEEALSO::
      :meth:`sage.combinat.posets.hasse_diagram.HasseDiagram.kappa`
    """

    if check and j not in self.join_irreducibles():
        raise ValueError(f"{j} is not join-irreducible.")

    hasse = self._hasse_diagram
    j_vtx = self._element_to_vertex(j)
    m_vtx = hasse.kappa(j_vtx)
    if m_vtx is None:
        return None
    m = self._vertex_to_element(m_vtx)
    return m

def extended_kappa_map(self,x):
    CJR = self.canonical_joinands(self)
    if CJR is None:
        raise ValueError(f"{x} does not have a canonical join representation.")
    kappa_CJR = [kappa_map(j) for j in CJR]
    try:
        return self.meet(kappa_CJR)
    except:
        raise ValueError("There is a canonical joinand for which the kappa map is not defined.")
