r"""
This code adds methods to
:class:`sage.combinat.posets.lattices.FiniteLatticePoset`
which compute the (extended) kappa (dual) map.
"""

from sage.combinat.posets.lattices import FiniteLatticePoset

# Algorithm:
# There is actually a method which computes kappa,
# :meth:`sage.combinat.posets.hasse_diagram.HasseDiagram.kappa`.
# The problem is there is no method for lattice,
# but only for Hasse diagrams.
# Thus, for a given lattice ``L``, we get its Hasse diagram
# ``L._hasse_diagram``, and just apply the existing method.


def kappa(self, j, check = True):
    r"""
    Return `\kappa(j)` for a join-irreducible element `j`
    in a finite lattice ``self`` if it exists.

    INPUT:

    - ``j`` -- an element of ``self``,
      which is expected to be join-irreducible

    - ``check`` -- a Boolean (default: ``True``),
      whether to check that ``j`` is indeed join-irreducible

    OUTPUT:

    an element of ``self``, or ``None`` if it does not exist.

    .. SEEALSO::
      :meth:`sage.combinat.posets.hasse_diagram.HasseDiagram.kappa`
    """

    if check and j not in self.join_irreducibles():
        raise ValueError("{0} is not join-irreducible.".format(j))

    hasse = self._hasse_diagram
    j_vtx = self._element_to_vertex(j)
    m_vtx = hasse.kappa(j_vtx)
    if m_vtx is None:
        return None
    m = self._vertex_to_element(m_vtx)
    return m

def extended_kappa(self, x):
    r"""
    Return `\overline{\kappa(x)}` for an element `x`
    in a finite lattice ``self`` if it exists.
    This first computes the canonical joinands of `x`,
    and then computes the meet of kappa of them.
    This returns ``None`` if ``x`` admits no canonical join reprensetation
    or kappa of some canonical joinand does not exist.

    INPUT:

    - ``x`` -- an element of ``self``

    OUTPUT:

    an element of ``self``, or ``None`` if it does not exist.

    .. SEEALSO::
      :func:`kappa`

    REFERENCE:

    .. [BCZ] E. Barnard, G. Todorov, S. Zhu,
      Dynamical combinatorics and torsion classes,
      J. Pure Appl. Algebra 225 (2021), no. 9, 106642.

    """
    CJR = self.canonical_joinands(x)
    if CJR is None:
        return None
    kappa_CJR = [self.kappa(j, check = False) for j in CJR]
    try:
        return self.meet(kappa_CJR)
    except:
        return None

def kappa_dual(self, m, check = True):
    r"""
    Return `\kappa^d(m)` for a meet-irreducible element `m`
    in a finite lattice ``self`` if it exists.

    INPUT:

    - ``m`` -- an element of ``self``,
      which is expected to be meet-irreducible

    - ``check`` -- a Boolean (default: ``True``),
      whether to check that ``m`` is indeed join-irreducible

    OUTPUT:

    an element of ``self``, or ``None`` if it does not exist.

    .. SEEALSO::
      :meth:`sage.combinat.posets.hasse_diagram.HasseDiagram.kappa_dual`
    """

    if check and m not in self.meet_irreducibles():
        raise ValueError("{0} is not join-irreducible.".format(m))

    hasse = self._hasse_diagram
    m_vtx = self._element_to_vertex(m)
    j_vtx = hasse.kappa_dual(m_vtx)
    if j_vtx is None:
        return None
    j = self._vertex_to_element(j_vtx)
    return j

def extended_kappa_dual(self, x):
    r"""
    Return `\overline{\kappa^d(x)}` for an element `x`
    in a finite lattice ``self`` if it exists.
    This first computes the canonical meetands of `x`,
    and then computes the join of kappa_dual of them.
    This returns ``None`` if `x` admits no canonical meet reprensetation
    or kappa_dual of some canonical meetand does not exist.

    INPUT:

    - ``x`` -- an element of ``self``

    OUTPUT:

    an element of ``self``, or ``None`` if it does not exist.

    .. SEEALSO::
      :func:`kappa_dual`

    REFERENCE:

    .. [BCZ] Emily Barnard, Gordana Todorov, Shijie Zhu,
      Dynamical combinatorics and torsion Classes,
      J. Pure Appl. Algebra 225 (2021), no. 9, 106642.
    """
    CMR = self.canonical_meetands(x)
    if CMR is None:
        return None
    kappa_d_CMR = [self.kappa_dual(m, check = False) for m in CMR]
    try:
        return self.join(kappa_d_CMR)
    except:
        return None

def add_kappa():
    r"""
    This function adds these functions as methods of a Sage class ``FiniteLatticePoset``.
    """
    FiniteLatticePoset.kappa = kappa
    FiniteLatticePoset.extended_kappa = extended_kappa
    FiniteLatticePoset.kappa_dual = kappa_dual
    FiniteLatticePoset.extended_kappa_dual = extended_kappa_dual
