--------------------------------
Computing with adèles and idèles
--------------------------------

This is a `SageMath <https://www.sagemath.org/>`__ package for computing with
adèles and idèles. It is based on and part of the master's thesis [Her2021].

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

In the root of this repository you can find [Her2021] as a PDF-file.


Contents of the package
-----------------------

The package can be seen to consist out of four parts.

Part 1 corresponds to Chapters 3--6 of [Her2021] and provides the functionality
to compute with adèles and idèles over number fields. It consists out of these
files:

- ``profinite_integer.py`` -- profinite integers over number fields
- ``profinite_number.py`` -- profinite numbers over number fields
- ``completion.py`` -- infinite completions of number fields
- ``adele.py`` -- adèles over number fields
- ``multiplicative_padic.py`` -- multiplicative `p`-adics
- ``idele.py`` -- idèles over number fields
- ``ray_class_group.py`` - ray class groups of number fields

Part 2 corresponds to Chapter 7 of [Her2021] and implements profinite graphs,
which visualize graphs of functions from and to the ring of rational profinite
integers. In particular, the profinite Fibonacci function is implemented. Part 2
consists of out two files:

- ``profinite_function.py`` -- profinite functions, including Fibonacci
- ``profinite_graph.py`` -- graphs of profinite functions

Part 3 corresponds to Chapter 8 of [Her2021] and implements the adèlic matrix
factorization algorithms discussed there. This resides in the file:

- ``matrix.py`` -- adèlic matrix factorization algorithms

Part 4 corresponds to Chapter 9 of [Her2021] and implements the computation of
Hilbert class fields of imaginary quadratic number fields using Shimura's
reciprocity law. It consists of the files:

- ``modular.py`` -- modular functions and their actions
- ``shimura.py`` -- Shimura's connecting homomorphism
- ``hilbert.py`` -- example hilbert class field computations


Getting acquainted with the package
-----------------------------------

Instead of browsing through the source code files, we recommend browsing the
documentation, which is nicer formatted. It contains many examples to illustrate
the functionality.


Documentation
-------------

The documentation resides in the folder ``docs`` and is also hosted online at
the following webpage: `<https://mathehertogh.github.io/adeles>`__.


Installing the package
----------------------

First of all you should make sure you have a recent version of `SageMath
<https://www.sagemath.org/download.html>`__ installed, specifically *SageMath
version 9.2 or newer*.

Now run the command ::

	$ sage -pip install adeles

To use the package, from anywhere on your computer, open ``sage`` ::

		$ sage

and within the ``sage`` prompt, load the package::

		sage: from adeles.all import *

Now you will have all functionality available, for example::

		sage: Adeles(QQ)
		Adèle Ring of Rational Field


Updating the package
--------------------

To update to the latest stable version of this package, run ::

	$ sage -pip install --upgrade adeles

It might be the case that the GitHub repository
`<https://github.com/mathehertogh/adeles>`__ contains an ever newer version.
To install that version, clone the repository ::

	$ git clone https://github.com/mathehertogh/adeles.git

change to the root directory of the package ::

	$ cd adeles

and build the package using ::

	$ make


Background information
----------------------

For more detailed information on this implementation of adèles and idèles, we
refer to [Her2021]. There we elaborate on properties of our representations of
adèles and idèles, design choices we made and implementation details.

For questions you can contact the author via email (see below).


Copyright
---------
::

	# **************************************************************************
	#       Copyright (C) 2021 Mathé Hertogh <m.c.hertogh@gmail.com>
	#
	# This program is free software: you can redistribute it and/or modify
	# it under the terms of the GNU General Public License as published by
	# the Free Software Foundation, either version 2 of the License, or
	# (at your option) any later version.
	#                  https://www.gnu.org/licenses/
	# **************************************************************************