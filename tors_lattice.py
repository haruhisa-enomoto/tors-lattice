import types
from kappa import kappa
from sage.combinat.posets.lattices import FiniteLatticePoset, LatticePoset

def TorsLattice(data = None, *args, **options):
    if isinstance(data, FiniteTorsLattice) and not args and not options:
        return data
    L = LatticePoset(data, *args, **options)
    if not L.is_semidistributive():
        raise ValueError("This lattice is not semidistributive.")
    return FiniteTorsLattice(L)


class FiniteTorsLattice(FiniteLatticePoset):

    def _kappa(self, j):
        return kappa(self, j, check = False)

    def is_brick_in(self, B, itv):
        U, T = itv
        return self.is_lequal(B,T) and self.is_gequal(self._kappa(B),U)

    def bricks(self, arg = None, *, torf = False):
        all_bricks = self.join_irreducibles()
        if arg is None:
            return all_bricks
        if isinstance(arg, (list,tuple)):
            itv = arg
        elif not torf:
            itv = (self.bottom(), arg)
        else:
            itv = (arg, self.top())
        return [B for B in all_bricks
                if self.is_brick_in(B, itv)]

    def is_wide_itv(self, itv, *, check = True):
        U,T = itv
        if not self.is_lequal(U,T):
            return False
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

    def wide(self):
        return [(self.minus(T), T) for T in self]

    def wide_lattice(self):
        wide_itvs = [(self.minus(T), T) for T in self]
        return LatticePoset((wide_itvs, self.itv_lequal))

    def wide_simple(self, itv):
        U, T = itv
        covers = [x for x in self.upper_covers(U) if self.is_lequal(x,T)]
        return [self.bricks((U,x))[0] for x in covers]

    def wide_lequal(self, itv1, itv2):
        return all([self.is_brick_in(S,itv2) for S in self.wide_simple(itv1)])

    def wide_lattice2(self):
        wide_itvs = [(self.minus(T), T) for T in self]
        return LatticePoset((wide_itvs, self.wide_lequal))


def temp(lattice, check = False):
    torsA = LatticePoset(lattice, facade = True)
    if check and not torsA.is_semidistributive():
        raise ValueError("This lattice is not semidistributive!")
