r"""
This module implements lattices of torsion classes of an abelian length category
in case they are finite.

This module defines :class: `FiniteTorsLattice`, which is a subclass of
:class:`sage.combinat.posets.lattices.FiniteLatticePoset`,
together with various methods which are useful in the study of
the representation theory of algebras.
"""
import sage.all
from sage.misc.cachefunc import cached_method
from sage.combinat.posets.lattices import FiniteLatticePoset, LatticePoset

def _kappa(lattice, j):
    r"""
    Return `\kappa(j)` for a join-irreducible element `j`
    in a finite lattice ``lattice`` if it exists.

    INPUT:

    - ``lattice`` -- a finite lattice,
      which should be an instance of
      :class:`sage.combinat.posets.lattices.FiniteLatticePoset`

    - ``j`` -- an element of ``lattice``,
      which is expected to be join-irreducible

    OUTPUT:

    an element of ``lattice``, or ``None`` if it does not exist.

    .. SEEALSO::
      :meth:`sage.combinat.posets.hasse_diagram.HasseDiagram.kappa`
    """
    hasse = lattice._hasse_diagram
    j_vtx = lattice._element_to_vertex(j)
    m_vtx = hasse.kappa(j_vtx)
    if m_vtx is None:
        return None
    m = lattice._vertex_to_element(m_vtx)
    return m

def _extended_kappa(lattice, x):
    r"""
    Return `\overline{\kappa(x)}` for an element `x`
    in a finite lattice ``lattice`` if it exists.

    This first computes the canonical joinands of `x`,
    and then computes the meet of kappa of them.
    This returns ``None`` if ``x`` admits no canonical join reprensetation
    or kappa of some canonical joinand does not exist.

    INPUT:

    - ``lattice`` -- a finite lattice,
      which should be an instance of
      :class:`sage.combinat.posets.lattices.FiniteLatticePoset`

    - ``x`` -- an element of ``lattice``

    OUTPUT:

    an element of ``lattice``, or ``None`` if it does not exist.

    REFERENCES:

    .. [BCZ] E. Barnard, G. Todorov, S. Zhu,
       Dynamical combinatorics and torsion classes,
       J. Pure Appl. Algebra 225 (2021), no. 9, 106642.

    """
    CJR = lattice.canonical_joinands(x)
    if CJR is None:
        return None
    kappa_CJR = [_kappa(lattice, j) for j in CJR]
    try:
        return lattice.meet(kappa_CJR)
    except:
        return None

def myshow(poset, label = True, vertex_size = 100, **kwargs):
    """
    A variant of ``show`` method with vertex labels,
    which looks nicer for the Hasse diagram of a poset.

    NOTE that the directions of Hasse arrows are opposite to
    the representation-threorist's convention, that is,
    every arrow starts at a smaller one and ends at the larger one.

    INPUT:

    - ``poset`` -- an object which has ``show`` method,
      which we expect to be an instance of
      :class:`sage.combinat.posets.posets.FinitePosets`

    - ``label`` -- a Boolean (default: ``True``), whether to label vertices

    """
    if label:
        poset.show(vertex_color = "white", vertex_shape = "_",
                   vertex_size = vertex_size, aspect_ratio = "automatic", **kwargs)
    else:
        poset.show(label_elements= False, vertex_size = vertex_size,
                   aspect_ratio = "automatic", **kwargs)

def TorsLattice(data = None, *args, **kwargs):
    """
    Construct a lattice of torsion classes from various forms of input data

    This raises an error if the constructed lattice is not semidistributive,
    since the lattice of torsion classes is semidistributive.

    INPUT:

    - ``data``, ``*args``, ``**kwargs`` -- data and options that will
      be passed down to :func:`LatticePoset` to construct a poset that is
      also a lattice.

    OUTPUT:

    An instance of :class:`FiniteTorsLattice`

    """
    if isinstance(data, FiniteTorsLattice) and not args and not kwargs:
        return data
    L = LatticePoset(data, *args, **kwargs)
    if not L.is_semidistributive():
        raise ValueError("This lattice is not semidistributive.")
    return FiniteTorsLattice(L)

class FiniteTorsLattice(FiniteLatticePoset):
    """
    A subclass of :class:`FiniteLatticePoset`,
    which we regard as the class of lattices of all torsion classes
    of an abelian length category.
    The argument passed to FiniteTorsLattice is assumed to be
    a finite semidistributive lattice.
    """

    def _repr_(self):
        return "Finite lattice of torsion classes of some abelian length category containing %s elements" % self._hasse_diagram.order()

    @cached_method
    def zero(self):
        """
        Return the smallest torsion class `0`
        """
        return self.bottom()

    @cached_method
    def whole(self):
        """
        Return the largest torsion class, i.e. the whole abelian category
        """
        return self.top()

    @cached_method
    def all_itvs(self):
        """
        Return the set of all intervals in the torsion poset.
        """
        return {(U,T) for U in self for T in self if self.is_lequal(U,T)}

    @cached_method
    def simples(self):
        """
        Return the list of simple torsion classes.

        Here a simple torsion class is a Serre subcategory which contains
        exactly one simple module, or equivalently,
        torsion classes from which there are arrows to `0`.
        We can use this to represent the list of simple modules.
        """
        return self.upper_covers(self.zero())

    @cached_method
    def all_bricks(self):
        """
        Return the set of all bricks represented by join-irreducibles.
        We always use join-irreducible torsion classes
        to represent bricks by a bijection in [DIRRT].

        REFERENCES:

        .. [DIRRT] L. Demonet, O. Iyama, N. Reading, I. Reiten, H. Thomas,
           Lattice theory of torsion classes, arXiv:1711.01785.

        """
        return set(self.join_irreducibles())

    @cached_method
    def kappa(self, T):
        r"""
        Return the (extended) kappa map of ``T``

        This is computed as follows:
        Let $B_1,\dots, B_k$ be brick labels of all arrows starting from $T$.
        Then $T$ is the join of $T(B_1), \dots, T(B_k)$,
        where $T(B)$ is the smallest torsion classes containing $B$.
        The kappa map $\kappa(T)$ is defined to be the intersection of
        $^\perp B_1, \dots, {}^\perp B_k$.
        This operation maps the canonical join representation
        to the canonical meet representation.
        See [BCZ] for the detail.

        INPUT:

        - ``T`` -- an element (torsion class) of ``self``

        OUTPUT:

        an element of ``self``

        REFERENCES:

        .. [BCZ] E. Barnard, G. Todorov, S. Zhu,
           Dynamical combinatorics and torsion classes,
           J. Pure Appl. Algebra 225 (2021), no. 9, 106642.
        """
        return _extended_kappa(self, T)

    @cached_method
    def bricks_in_tors(self, T):
        r"""
        Return the frozenset of bricks contained in a torsion class ``T``

        INPUT:

        - ``T`` -- an element (torsion class) of ``self``

        OUTPUT:

        the frozenset of bricks (represented by join-irreducibles in ``self``)
        contained in ``T``

        """
        return frozenset({j for j in self.all_bricks() if self.is_lequal(j,T)})

    @cached_method
    def bricks_in_torf(self, T):
        r"""
        Return the frozenset of bricks contained in a torsion-free class $T^\perp$

        INPUT:

        - ``T`` -- an element (torsion class) of ``self``

        OUTPUT:

        the frozenset of bricks (represented by join-irreducibles in ``self``)
        contained in the torsion-free class corresponding to ``T``,
        i.e. $T^\perp$
        """
        return frozenset(j for j in self.all_bricks()
                         if self.is_gequal(self.kappa(j),T) )

    @cached_method
    def bricks(self, itv, *, check = True):
        r"""
        Return the set of bricks in the heart of an interval of torsion classes

        For two torsion classes $U,T$ with $U \subseteq T$,
        its heart is $T \cap U^\perp$ (see [ES]).

        INPUT:

        - ``itv`` -- a pair (tuple or list) of torsion classes

        - ``check`` -- a Boolean (default: ``True``),
          whether to check ``itv`` is actually an interval

        OUTPUT:

        the frozenset of bricks (represented by join-irreducibles)
        contained in the heart of the given interval

        REFERENCES:

        .. [ES] H. Enomoto, A. Sakai,
           ICE-closed subcategories and wide $\tau$-tilting modules,
           to appear in Math. Z.
        """
        U, T = itv
        if check and not self.is_lequal(U,T):
            raise ValueError("This is not an interval.")
        return self.bricks_in_tors(T) & self.bricks_in_torf(U)

    @cached_method
    def label(self, itv, *, check = True):
        r"""
        Return the brick label of an Hasse arrow in the lattice of torsion classes
        For a Hasse arrow $T \to U$, its label is a unique brick
        contained in $T \cap U^\perp$ [DIRRT].

        INPUT:

        - itv -- a pair (tuple or list) of torsion classes (U,T),
          which we expect that U is covered by T

        - ``check`` -- a Boolean (default: ``True``),
          whether to check ``itv`` actually gives a covering relation

        REFERENCES:

        .. [DIRRT] L. Demonet, O. Iyama, N. Reading, I. Reiten, H. Thomas, Lattice theory of torsion classes, arXiv:1711.01785.
        """
        bricks_in = self.bricks(itv)
        if check and len(bricks_in) > 1:
            raise ValueError("The heart contains more than one brick, \
                             so not a covering relation.")
        return list(bricks_in)[0]

    @cached_method
    def plus(self, U):
        """
        Return the join of all Hasse arrows ending at ``U``

        For a torsion class $U$, its plus $U^{+}$ satisfies that
        $[U,U^{+}]$ is a wide interval which is the largest wide interval
        of the form $[U,T]$.

        INPUT:

        - U -- an element (torsion class) of ``self``
        """
        if U == self.whole():
            return U
        upper = self.upper_covers(U)
        return self.join(upper)

    @cached_method
    def minus(self, T):
        """
        Return the meet of all Hasse arrows starting at ``T``

        For a torsion class $T$, its minus $T^{-}$ satisfies that
        $[T^{-},T]$ is a wide interval which is the largest wide interval
        of the form $[U,T]$.

        INPUT:

        - T -- an element (torsion class) of ``self``
        """
        if T == self.zero():
            return T
        lower = self.lower_covers(T)
        return self.meet(lower)

    def is_wide_itv(self, itv, *, check = True):
        r"""
        Return ``True`` if ``itv`` is a wide interval, and ``False`` otherwise

        An interval $[U,T]$ is a wide interval if its heart
        $T \cap U^\perp$ is a wide subcategory.
        This function uses the characterization of wide intervals
        given in [AP].

        INPUT:

        - itv -- a pair (list or tuple) of torsion classes,
          which is expected to be an interval

        - ``check`` -- a Boolean (default: ``True``),
          whether to check ``itv`` is actually an interval

        REFERENCES:

        .. [AP] S. Asai, C. Pfeifer,
           Wide subcategories and lattices of torsion classes,
           arXiv:1905.01148.
        """
        U,T = itv
        if check and not self.is_lequal(U,T):
            raise ValueError("This is not an interval.")
        if U == T:
            return True
        covers = [x for x in self.upper_covers(U) if self.is_lequal(x,T)]
        return T == self.join(covers)

    def is_ice_itv(self, itv, *, check = True):
        r"""
        Return ``True`` if ``itv`` is an ICE interval, and ``False`` otherwise

        An interval $[U,T]$ is a wide interval if its heart
        $T \cap U^\perp$ is an ICE-closed subcategory, that is,
        closed under taking images, cokernels, and extensions.
        This function uses the characterization of ICE intervals
        given in [ES].

        INPUT:

        - itv -- a pair (list or tuple) of torsion classes,
          which is expected to be an interval

        - ``check`` -- a Boolean (default: ``True``),
          whether to check ``itv`` is actually an interval

        REFERENCES:

        .. [ES] H. Enomoto, A. Sakai,
           ICE-closed subcategories and wide $\tau$-tilting modules,
           to appear in Math. Z.
        """
        U,T = itv
        if check and not self.is_lequal(U,T):
            raise ValueError("This is not an interval.")
        return self.is_lequal(T,self.plus(U))

    def is_ike_itv(self, itv, *, check = True):
        r"""
        Return ``True`` if ``itv`` is an IKE interval, and ``False`` otherwise

        An interval $[U,T]$ is a wide interval if its heart
        $T \cap U^\perp$ is an IKE-closed subcategory, that is,
        closed under taking images, kernels, and extensions.
        This function is just a dual of :func:`is_ice_itv`.

        INPUT:

        - itv -- a pair (list or tuple) of torsion classes,
          which is expected to be an interval

        - ``check`` -- a Boolean (default: ``True``),
          whether to check ``itv`` is actually an interval
        """
        U,T = itv
        if check and not self.is_lequal(U,T):
            raise ValueError("This is not an interval.")
        return self.is_lequal(self.minus(T),U)

    def itv_lequal(self, itv1, itv2):
        r"""
        Return whether the heart of ``itv1`` is contained in ``itv2``

        The heart of an interval $[U,T]$ is a subcategory $T \cap U^\perp$.
        By [DIRRT], the heart is recovered from bricks contained in it,
        hence this function compare the sets of bricks in two hearts.

        INPUT:

        - ``itv1``, ``itv2`` -- pairs (list or tuple) of torsion classes,
          which are assumed to be intervals

        OUTPUT:

        ``True`` if the heart of ``itv1`` is contained in that of ``itv2``,
        and ``False otherwise.
        """
        return set(self.bricks(itv1)).issubset(self.bricks(itv2))

    def wide_simples(self, itv):
        """
        Return the list of simple objects in a wide subcategory corresponding to ``itv``

        INPUT:

        - ``itv`` -- a pair of torsion class, which we assume is a wide interval

        OUTPUT:

        the list of simple objects (bricks represented by join-irreducibles) of a wide subcategory which is the heart of ``itv``
        """
        if not self.is_wide_itv(itv):
            raise ValueError("This interval is not a wide interval.")

        U, T = itv
        covers = [x for x in self.upper_covers(U) if self.is_lequal(x,T)]
        return [self.label((U,x), check = False) for x in covers]

    def wide_lequal(self, U, T):
        r"""
        Compare two wide subcategories corresponding to two torsion classes

        If there are only finitely many torsion classes, then there is a bijection
        between the set of torsion classes and the set of wide subcategories
        by [MS] for finite-dimensional algebras and [E] for an abelian length category.
        Write $W_L(T)$ for the wide subcategory corresponding to $T$,
        which is a filtration closure of the brick labels of all Hasse arrows starting at $T$.
        Then this function returns whether $W_L(U) \subseteq W_L(T)$.

        INPUT:

        - U, T -- elements of ``self`` (torsion classes)

        OUTPUT:

        ``True`` if $W_L(U)$ is contained in $W_L(T)$, and ``False`` otherwise

        REFERENCES:

        .. [MS] F. Marks, J. Stovicek,
           Torsion classes, wide subcategories and localisations,
           Bull. London Math. Soc. 49 (2017), Issue 3, 405–416.

        .. [E] H. Enomoto,
           Monobrick, a uniform approach to torsion-free classes and wide subcategories,
           arXiv:2005.01626.
        """
        # U, T represents W_L(U)
        kappa_U, kappa_T = self.kappa(U), self.kappa(T)
        return self.is_lequal(U, T) and self.is_gequal(kappa_U, kappa_T)

    def wide_lattice(self):
        """
        Return the lattice of wide subcategories
        """
        return LatticePoset((self, self.wide_lequal))

    def ice_lattice(self):
        """
        Return the lattice of ICE-closed subcategories, that is,
        subcategories closed under images, cokernels, and extensions.

        REFERENCES:

        .. [ES] H. Enomoto, A. Sakai,
           ICE-closed subcategories and wide $\tau$-tilting modules,
           to appear in Math. Z.
        """
        ice_bricks = {self.bricks(itv, check = False) for itv in self.all_itvs()
                      if self.is_ice_itv(itv, check = False)}
        return LatticePoset((ice_bricks, attrcall("issubset")))

    def ike_lattice(self):
        """
        Return the lattice of IKE-closed subcategories, that is,
        subcategories closed under images, kernels, and extensions.

        This is just a dual of :func:`ice_lattice`.
        """
        ike_bricks = {self.bricks(itv, check = False) for itv in self.all_itvs()
                      if self.is_ike_itv(itv, check = False)}
        return LatticePoset((ike_bricks, attrcall("issubset")))

    def heart_poset(self):
        """
        Return the poset of torsion hearts ordered by inclusion

        A torsion heart is a subcategory which arises as a heart of some interval
        of torsion classes. For example, every wide subcategory, ICE-closed subcategory are
        torsion hearts.

        This poset is not a lattice in general.
        """
        brick_set = {self.bricks(itv, check = False) for itv in self.all_itvs()}
        return Poset((brick_set, attrcall("issubset")))

    def num_of_proj(self, arg):
        r"""
        Return the number of indecomposable Ext-projective objects in a subcategory

        If ``arg`` is an element of ``self`` (i.e. a torsion class),
        then the considered category is ``arg`` itself.
        If ``arg`` is an interval $[U, T]$ of torsion classes, then this considers
        the heart of this interval, i.e. $T \cap U^\perp$.
        So far, we can only deal with the case where the heart is ICE-closed or IKE-closed
        (closed under images, (co)kernels, and extensions).

        INPUT:

        - ``arg`` -- either an element of ``self``, or an interval in ``self`` as a tuple or list
        """
        if isinstance(arg, (list,tuple)):
            U, T = arg
        else:
            U, T = self.zero(), arg
        if not self.is_ice_itv([U,T]) and not self.is_ike_itv([U,T]):
            raise NotImplementedError("The heart of this interval is not ICE-closed or IKE-closed.")
        if U == T:
            return 0
        if self.is_ice_itv([U,T]):
            # Here by Asai-Pfeifer's reduction,
            # we pretend that U = 0 and compute the number of projs in a torsion class T.
            # The idea is that the number of projs is equal to the number of supports,
            # hence we compute the number of supports by constructing Serre subcategories.
            simples = set(self.upper_covers(U))
            n = len(simples)
            for i in range(1,n+1):
                # First compute Serre subcats consisting of i simples
                for set_of_simples in combinations(simples, i):
                    serre = self.join(set_of_simples)
                    # If T is contained in this Serre subcats,
                    # then the number of supports of T is i,
                    # hence there are i projectives.
                    if self.is_lequal(T, serre):
                        return i
        else:
            # This is IKE case, and is completely dual to ICE case.
            simples = set(self.lower_covers(T))
            n = len(simples)
            for i in range(1,n+1):
                for set_of_simples in combinations(simples, i):
                    serre = self.meet(set_of_simples)
                    if self.is_lequal(serre, U):
                        return i

    def indec_tau_rigid(self):
        r"""
        Return the set of indecomposable $\tau$-rigid modules,
        represented by join-irreducible torsion classes.

        For a $\tau$-tilting finite algebra, there is a bijection by [DIJ] between
        indecomposable $\tau$-rigid module $M$ and join-irreducible torsion classes.
        The correspondence is $T(M) = \mathsf{Fac} M$ for $\tau$-rigid $M$, and
        the unique split projective object in $T$ for a join-irreducible torsion class $T$.

        Since this is the same as :func:`all_bricks`,
        this function is only needed for the readability reason.
        """
        return self.all_bricks()

    @cached_method
    def has_tau_rigid_summand(self, M, *, check = True):
        r"""
        Return the set of $\tau$-tilting pairs which has ``M`` as a $\tau$-tilting summand

        We consider ``self`` as the set of support $\tau$-tilting pairs.
        Then this returns the set of support $\tau$-tilting pairs
        which contain $(M,0)$ as a direct summand.
        We use join-irreducible torsion classes to represent indecomposable $\tau$-rigid modules.
        See :meth:`indec_tau_rigid`.

        INPUT:

        - ``M`` -- an element of ``self``, which is expected to be join-irreducible

        - ``check`` -- a Boolean (default: ``True``),
          whether to check ``M`` is actually join-irreducible
        """
        if check and M not in self.indec_tau_rigid():
            raise ValueError("The first argument is not join-irreducible.")
        M_plus = self.plus(M)
        # then [M, M_plus] is the set of tau-tilting pairs
        # containing (M,0) as a summand.
        return {T for T in self if self.is_lequal(M,T) and self.is_lequal(T,M_plus) }

    @cached_method
    def has_support_summand(self, S, *, check = True):
        r"""
        Return the set of $\tau$-tilting pairs which have the projective cover of ``S`` as a support summand

        We consider ``self`` as the set of support $\tau$-tilting pairs.
        Then this returns the set of support $\tau$-tilting pairs
        which contain $(0,P)$ as a direct summand,
        where $P$ is the projective cover of a simple module ``S``.
        We use simple torsion classes to represent a simple module,
        hence ``S`` is expected to be a simple torsion class, that is,
        a Serre subcategory with one simple module.

        INPUT:

        - ``S`` -- an element of ``self``, which is expected to be a simple torsion class

        - ``check`` -- a Boolean (default: ``True``),
          whether to check ``S`` is actually a simple torsion class
        """
        if check and S not in self.simples():
            raise ValueError("This is not a simple torsion class (doesn't cover 0).")
        # We compute the Serre subcategory ``non_S_Serre`` consisting of modules
        # such that ``S`` don't appear as composition factors.
        non_S_simple = set(self.simples()) - {S}
        non_S_Serre = self.join(non_S_simple)
        # Then a $\tau$-tilting pair contains $(0,P)$ as a summand
        # if there's no non-zero map from P to any element in the corresponding torsion class $T$,
        # that is, $T$ is contained in ``non_S_Serre``.
        return {T for T in self if self.is_lequal(T,non_S_Serre)}

    def tau_rigid_summand(self, T):
        r"""
        Return the set of indecomposable $\tau$-tilting pairs which are direct summands of ``T``

        We represent indecomposable $\tau$-tilting pairs as follows.
        - For a pair $(M,0)$ with $M$ being indecomposable $\tau$-rigid,
          we use ``(M,0)`` for ``M`` in ``self.indec_tau_rigid()``,
          that is, ``M`` is the join-irreducible torsion class corresponding to $M$.
        - For a pair $(0,P)$ with $P$ being indecomposable projective,
          we use ``(S,1)``, where ``S`` is a simple module $\mathrm{top} P$
          represented by the simple torsion class.
        """
        return {(M,0) for M in self.indec_tau_rigid()
                if T in self.has_tau_rigid_summand(M, check = False)} | \
               {(S,1) for S in self.simples()
                if T in self.has_support_summand(S, check = False)}

    def support_tau_tilting_complex(self):
        r"""
        Return the support $\tau$-tilting simplicial complex of the algebra

        OUTPUT:

        an instance of :clas:`sage.homology.simplicial_complex.SimplicialComplex`
        """
        return SimplicialComplex([self.tau_rigid_summand(T) for T in self],
                                 maximality_check = False)

    def num_of_proj(self, arg):
        r"""
        Return the number of indecomposable Ext-projective objects in a given subcategory

        If ``arg`` is an element of ``self`` (i.e. a torsion class),
        then the considered category is ``arg`` itself.
        If ``arg`` is an interval $[U, T]$ of torsion classes, then this considers
        the heart of this interval, i.e. $T \cap U^\perp$.
        This is based on [ES, Corollary 4.28]

        INPUT:

        - ``arg`` -- either an element of ``self``, or an interval in ``self`` as a tuple or list

        REFERENCES:

        .. [ES] H. Enomoto, A. Sakai,
           ICE-closed subcategories and wide $\tau$-tilting modules,
           to appear in Math. Z.
        """
        if isinstance(arg, (list,tuple)):
            U, T = arg
        else:
            U, T = self.zero(), arg
        projs_in_T = {M for M in self.indec_tau_rigid()
                      if T in self.has_tau_rigid_summand(M, check = False)}
        T_minus_U = {M for M in projs_in_T if not self.is_lequal(M,U)}
        return len(T_minus_U)
