--------------------------------
Computing with adèles and idèles
--------------------------------

This is a `SageMath <https://www.sagemath.org/>`_ package for computing with
adèles and idèles. It is based on and part of the master's thesis [Her2021].

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.


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
integers. In particular, the profinite Fibonacci function is implemnted. Part 2
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


How to view the contents and discover the functionality
-------------------------------------------------------

First of all, we recommomend the user to read/get acquainted to the files in the
same order as we listed them above. Secondly, we recommend reading the
documentation instead of the source code. The documentation can be generated
from the source code as follows.

Make sure you have `(a recent version of) SageMath
<https://www.sagemath.org/download.html>`_ installed.

Download the contents of this repository and save it somewhere on your drive,
say in the directory ``/path/to/adeles``.

In order to generate the HTML documentation (which we prefer) of
``profinite_integer.py``, perform the command ::

	$ sage --docbuild "file=/path/to/adeles/profinite_integer.py" html

The output of this command will look like this::

	$ ...
	$ Build finished. The built documents can be found in /home/mathe/.sage/docbuild/profinite_integer/output/html

Now open the file ``/home/mathe/.sage/docbuild/profinite_integer/output/html/index.html``
in a browser to view the documentation of ``profinite_integer.py``.

How to use this package
-----------------------

We assume you have a working SageMath installation.

To use the




Copyright
---------
::

	# ****************************************************************************
	#       Copyright (C) 2021 Mathé Hertogh <m.c.hertogh@gmail.com>
	#
	# This program is free software: you can redistribute it and/or modify
	# it under the terms of the GNU General Public License as published by
	# the Free Software Foundation, either version 2 of the License, or
	# (at your option) any later version.
	#                  https://www.gnu.org/licenses/
	# ****************************************************************************