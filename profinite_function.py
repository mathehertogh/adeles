r"""
Profinite Functions

This file describes an interface for implementations of profinite functions
`\hat{\ZZ} \to \hat{\ZZ}` in the form of the abstract class
:class:`ProfiniteFunction`.

It also implements a good example of such a profinite function, namely the
Profinite Fibonacci function, in :class:`ProfiniteFibonacci`.

REFERENCES:

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

These profinite functions, aspecially the profinite Fibonacci function, is
introduced in Chapter 7 of [Her2021].

AUTHORS:

- Mathé Hertogh (2021-7): initial version based on [Her2021]

"""

# ****************************************************************************
#       Copyright (C) 2021 Mathé Hertogh <m.c.hertogh@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.matrix.constructor import matrix
from sage.arith.misc import gcd
from sage.rings.integer_ring import ZZ
from profinite_integer import Zhat
from sage.misc.abstract_method import abstract_method


class ProfiniteFunction:
    r"""
    Abstract class representing a function from and to the profinite integers

    An instance `P` of this class represents a function `F: \hat{\ZZ} \to
    \hat{\ZZ}` and should implement the abstract method :meth:`__call__` to
    evaluate this function in a point. For a profinite `\QQ`-integer `x` the
    call `P(x)` should return a profinite `\QQ`-inter and this evaluation must
    satisfy the following conditions:

    1. for every profinite `\QQ`-integer `x` and for every `\alpha \in
       \hat{\ZZ}` that `x` represents, `F(\alpha)` is represented by `P(x)`;
    2. for every `k \in \ZZ_{>0}` there exists `m \in \ZZ_{\geq k}` such
       that for all `n \in \ZZ_{\geq m}` and for all profintie `\QQ`-integers
       `x` of modulus divisible by `n!`, the modulus of `P(x)` is divisible by
       `k!`.

    Condition 1 ensures that the evaluation is "correct" and Condition 2 states
    that the profinite function can be evulated with "arbitrarily high
    precision".
        
    .. SEEALSO::
        
        See the class :class:`profinite_graph.ProfiniteGraph` for visualizations
        of the graphs of profinite functions.

    .. automethod:: __call__
    """

    @abstract_method
    def __call__(self, x, des_mod=0):
        r"""
        Evaluate this function in ``x``

        This is an abstract method that should be implemented by subclasses.

        The evulation should always satisfy Condition 1 stated in the docstring
        of the class :class:`ProfiniteFunction`.
        For ``des_mod == 0`` it should also satisfy Condition 2.

        INPUT:

        - ``x`` -- profinite `\QQ`-integer; point to evaluate this profinite
          function in
        - ``des_mod`` -- integer (default: `0`); the desired modulus of
          the output. Set to 0 for "as high as possible".

        OUTPUT:
        
        The profinite integer ``self(x)``. If ``des_mod`` is set to a non-zero
        value, then this function returns the profinite integer
        ``self(x) + Zhat(0, m)`` for some multiple `m` of ``des_mod``.

        .. NOTE::

            The ``des_mod`` option can be used by implementers of profinite
            functions to speed up the evaluation: the user can specify a 
            precision up to which he/she is interested in the result. If this
            desired output precision is much lower than the maximal precision
            (returend for ``des_mod == 0``), the evaluation of the profinite
            function may be much faster. For a good example of such a situation,
            see :class:`ProfiniteFibonacci`.

        EXAMPLES:

        We create a profinite function implementing the constant zero map
        `\hat{\ZZ} \to \hat{\ZZ}, x \mapsto 0`. We only return zero modulo the
        desired modulus, although we know the result is exactly zero::

            sage: class ConstantZero(ProfiniteFunction):
            ....:     def __call__(self, x, des_mod=0):
            ....:         return Zhat(0, des_mod)
            sage: zero = ConstantZero()
            sage: zero(Zhat(5,20), 0)
            0 mod 0
            sage: zero(Zhat(5,20), 30)
            0 mod 30
        
        Next we create a profinite function implementing the squaring map
        `\hat{\ZZ} \to \hat{\ZZ}, x \mapsto x^2`. This time we simply ignore
        the desired input modulus given by the user. ::

            sage: class Square(ProfiniteFunction):
            ....:     def __call__(self, x, des_mod=0):
            ....:         return x * x
            sage: square = Square()
            sage: square(Zhat(2, 20), 0)
            4 mod 40
            sage: square(Zhat(2, 20), 15)
            4 mod 40

        In the last example, ``4 mod 5`` would also have been an acceptable
        answer, since ``gcd(40, 15)=5``. But instead, we choose to return an
        answer with strictly more precision: the multiple 40 of 5.

        For a good "real world" example, take a look at the
        :class:`profinite Fibonacci function <ProfiniteFibonacci>`.
        """
        raise NotImplementedError("no general function evaluation implemented")


class ProfiniteFibonacci(ProfiniteFunction):
    r"""
    The profinite Fibonacci function

    The Fibonacci function `F` is defined on the integers by `F(0)=0`, `F(1)=1`
    and `F(n) = F(n-1) + F(n-2)` for all integers  `n \in \ZZ`.
    There exists a unique continuous extension `F: \hat{\ZZ} \to \hat{\ZZ}` of
    the Fibonacci function, called the *profinite Fibonacci function*. This is
    an implementation of it.

    .. automethod:: __call__
    """
    def __call__(self, x, des_mod=0):
        r"""
        Return the ``x``-th profinite Fibonacci number, for ``x`` a profinite
        integer

        We refer to :meth:`ProfiniteFunction.__call__
        <sage.rings.adeles.profinite_functions.ProfiniteFunction.__call__>`
        for a description of the input and output.

        The modulus of the output will always be a divisor of ``des_mod``, so we
        never compute up to a higher precision than the user asks for.

        EXAMPLES::

            sage: fibonacci = ProfiniteFibonacci()
            sage: fibonacci(Zhat(3, 10))
            2 mod 11
            sage: fibonacci(Zhat(4, 18))
            3 mod 76
            sage: fibonacci(Zhat(5, 100))
            5 mod 12586269025
            sage: fibonacci(Zhat(5, 100), 100)
            5 mod 25
            sage: fibonacci(Zhat(6, 10^7), 10^10)
            8 mod 78125

        All computations above finish "instantly". But if we remove the
        ``des_mod=10^10`` in the last example, the computation takes multiple
        seconds. See the warning below.

        ALGORITHM:
        
        Write `F` for the Fibonacci function `\ZZ \to \ZZ`. Let `m` be the
        modulus of ``x``. Then the modulus of the output is computed as
        `\gcd(F(m), F(m+1)-1)`. The value of the output is computed as `F(v)`
        where `v` denotes the value of ``x``.

        These Fibonacci numbers in `\ZZ` are computed using the formula

        .. MATH::

            \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}^n
            =
            \begin{pmatrix} F(n+1) & F(n) \\ F(n) & F(n-1) \end{pmatrix}

        The matrix exponentiation is done using a square-and-multiply method.
        For `F(m)` and `F(m+1)` we do this exponentiation in the ring of
        matrices over `\ZZ/\texttt{des_mod}\ZZ`; for `F(v)` over `\ZZ/N\ZZ`,
        where `N` denotes the output modulus.

        .. WARNING::

            Not specifying ``des_mod``, or setting it to 0, requires the exact
            computation of ``F(x.modulus())``, which may be huge (`F(m) \approx
            1.6^m`). Hence the computation with bounded desired modulus is
            significantly faster than for unbounded.
        
        .. SEEALSO::
            
            See the class :class:`profinite_graph.ProfiniteGraph` for a
            visualization of the graph of this function.
        """
        if des_mod is None:
            des_mod = x.modulus()
        R = ZZ if des_mod == 0 else IntegerModRing(des_mod)
        # We compute the (input modulus)-th and the (input modulus + 1)-th
        # Fibonacci numbers (F_im0 and F_im1 respectively) in R.
        A = matrix(R, [[1, 1], [1, 0]])
        F_input_mod = A**x.modulus()
        F_im0 = F_input_mod.coefficient((1, 0))
        F_im1 = F_input_mod.coefficient((0, 0))
        if des_mod != 0:
            # If R=ZZ/des_mod*ZZ, lift our results back to ZZ.
            F_im0 = F_im0.lift()
            F_im1 = F_im1.lift()
        # Based on these we calculate the modulus of our output.
        output_modulus = gcd(des_mod, F_im0 - 0)
        output_modulus = gcd(output_modulus, F_im1 - 1)
        # Now we calculate the (x.value())-th Fibonacci number modulo the
        # output_modulus.
        A = matrix(IntegerModRing(output_modulus), [[1, 1], [1, 0]])
        output_value = (A**x.value()).coefficient((1, 0)).lift()
        return Zhat(output_value, output_modulus)
