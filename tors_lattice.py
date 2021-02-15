r"""
This module implements lattices of torsion classes of an abelian length category
in case they are finite.

This module defines :class: `FiniteTorsLattice`, which is a subclass of
:class:`sage.combinat.posets.lattices.FiniteLatticePoset`.
"""

from sage.combinat.posets.lattices import FiniteLatticePoset, LatticePoset

def kappa(lattice, j, *, check = True):
    r"""
    Return `\kappa(j)` for a join-irreducible element `j`
    in a finite lattice ``self`` if it exists.

    INPUT:

    - ``lattice`` -- a finite lattice

    - ``j`` -- an element of ``self``,
      which is expected to be join-irreducible

    - ``check`` -- a Boolean (default: ``True``),
      whether to check that ``j`` is indeed join-irreducible

    OUTPUT:

    an element of ``lattice``, or ``None`` if it does not exist.

    .. SEEALSO::
      :meth:`sage.combinat.posets.hasse_diagram.HasseDiagram.kappa`
    """

    if check:
        if j not in lattice:
            raise ValueError("{0} is not an element of {1}".format(j, lattice))
        elif j not in lattice.join_irreducibles():
            raise ValueError("{0} is not join-irreducible.".format(j))

    hasse = lattice._hasse_diagram
    j_vtx = lattice._element_to_vertex(j)
    m_vtx = hasse.kappa(j_vtx)
    if m_vtx is None:
        return None
    m = lattice._vertex_to_element(m_vtx)
    return m


def TorsLattice(data = None, *args, **kwargs):
    if isinstance(data, FiniteTorsLattice) and not args and not kwargs:
        return data
    L = LatticePoset(data, *args, **kwargs)
    if not L.is_semidistributive():
        raise ValueError("This lattice is not semidistributive,
                         thus cannot be a lattice of torsion classes.")
    return FiniteTorsLattice(L)

class FiniteTorsLattice(FiniteLatticePoset):

    def _kappa(self, j):
        return kappa(self, j, check = False)

    def is_brick_in(self, B, itv):
        U, T = itv
        return self.is_lequal(B,T) and self.is_gequal(self._kappa(B),U)

    def bricks(self, arg = None, *, tors = True):
        all_bricks = self.join_irreducibles()
        if arg is None:
            return all_bricks
        if isinstance(arg, (list,tuple)):
            itv = arg
        elif tors:
            itv = (self.bottom(), arg)
        else:
            itv = (arg, self.top())
        return [B for B in all_bricks if self.is_brick_in(B, itv)]

    def is_wide_itv(self, itv, *, check = True):
        U,T = itv
        if check and not self.is_lequal(U,T):
            raise ValueError("The second argument is not an interval.")
        if U == T:
            return True
        covers = [x for x in self.upper_covers(U) if self.is_lequal(x,T)]
        return T == self.join(covers)

    def itv_lequal(self, itv1, itv2):
        return set(self.bricks(itv1)).issubset(self.bricks(itv2))

    def minus(self, T):
        if T == self.bottom():
            return T
        lower = self.lower_covers(T)
        return self.meet(lower)

    def plus(self, U):
        if U == self.top():
            return U
        upper = self.upper_covers(U)
        return self.join(upper)

    def wide_simple(self, itv, *, check = True):
        if check and not self.is_wide_itv(itv):
            raise ValueError("This interval is not a wide interval.")

        U, T = itv
        covers = [x for x in self.upper_covers(U) if self.is_lequal(x,T)]
        return [self.bricks((U,x))[0] for x in covers]

    def wide_lequal(self, itv1, itv2):
        return all([self.is_brick_in(S,itv2) for S in self.wide_simple(itv1)])

    def wide_lattice(self):
        wide_itvs = [(self.minus(T), T) for T in self]
        return LatticePoset((wide_itvs, self.wide_lequal))
