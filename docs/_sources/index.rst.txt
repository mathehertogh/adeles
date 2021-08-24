--------------------------------
Computing with adèles and idèles
--------------------------------

A `SageMath <https://sagemath.org>`_ package for computing with adèles and
idèles.

To use this package, you need to import it::

    sage: from adeles.all import *

After this, all functionality of this package is available::

    sage: Adeles(QQ)
    Adèle Ring of Rational Field

This package is based on and part of the following master's thesis:

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.


Contents of the package
-----------------------

The package can be seen to consist out of four parts.

Part 1 corresponds to Chapters 3--6 of [Her2021] and provides the functionality
to compute with adèles and idèles over number fields. It consists out of these
files:

    - :doc:`profinite_integer`
    - :doc:`profinite_number`
    - :doc:`completion`
    - :doc:`adele`
    - :doc:`multiplicative_padic`
    - :doc:`idele`
    - :doc:`ray_class_group`

Part 2 corresponds to Chapter 7 of [Her2021] and implements profinite graphs,
which visualize graphs of functions from and to the ring of rational profinite
integers. In particular, the profinite Fibonacci function is implemented. Part 2
consists of out two files:

    - :doc:`profinite_function`
    - :doc:`profinite_graph`

Part 3 corresponds to Chapter 8 of [Her2021] and implements the adèlic matrix
factorization algorithms discussed there. This resides in the file:


    - :doc:`matrix`

Part 4 corresponds to Chapter 9 of [Her2021] and implements the computation of
Hilbert class fields of imaginary quadratic number fields using Shimura's
reciprocity law. It consists of the files:

    - :doc:`modular`
    - :doc:`shimura`
    - :doc:`hilbert`


Table of Contents
----------------------------------

.. toctree::

    profinite_integer
    profinite_number
    completion
    adele
    multiplicative_padic
    idele
    ray_class_group
    profinite_function
    profinite_graph
    matrix
    modular
    shimura
    hilbert


Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
