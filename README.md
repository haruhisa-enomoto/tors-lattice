# Lattices of torsion classes in SageMath
[`tors_lattice.py`](tors_lattice.py) in this repository is for [SageMath](https://www.sagemath.org/).
This module can deal with the lattice of torsion classes of &tau;-tilting finite algebras and construct various objects from it.

## Overview

This module defines a class `FiniteTorsLattice`, which is a subclass of a SageMath class for finite lattices: [`FiniteLatticePoset`](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/posets/lattices.html#sage.combinat.posets.lattices.FiniteLatticePoset).
`FiniteTorsLattice` is a class for lattices of torsion classes over &tau;-tilting fintie artin algebras. Once you input the lattice of torsion classes (e.g. using my [StringApplet-to-SageMath-converter](https://github.com/haruhisa-enomoto/StringApplet-to-SageMath-converter)), then this program can compute various objects naturally arising in the representation theory of algebras, such as the lattice of wide subcategories and support &tau;-tilting simplicial complex, and the number of indecomposable Ext-projectives of each torsion class, and so on.

### Author
[Haruhisa Enomoto](http://haruhisa-enomoto.github.io/), a postdoc. at Osaka Prefecture University in Japan.

### Requirements
[SageMath](https://www.sagemath.org/) version 9.x or later (since this code is based on Python 3)

## Definitions and Usage

- [Manual](https://nbviewer.jupyter.org/github/haruhisa-enomoto/tors-lattice/blob/main/Manual.ipynb).

### References

- [E] H. Enomoto,
  *Computing various objects of an algebra from the poset of torsion classes*,
  in preparation.
