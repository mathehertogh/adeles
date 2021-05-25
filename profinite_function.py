r"""
Profinite Function

AUTHORS:

- Mathé Hertogh (2021-01-15): initial version

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
    """
    Abstract class representing a function from and to the profinite integers


    This class has one abstract method, the __call__() method.

    .. automethod:: __call__
    """

    @abstract_method
    def __call__(self, x, des_mod=0):
        r"""
        Evaluate this function in ``x``

        This is an abstract method.

        Denote the function that we implement by `F: \hat{\ZZ} \to \hat{\ZZ}`.

        INPUT:

        - ``x`` -- profinite integer; point to evaluate `F` in
        - ``des_mod`` -- integer (default: `0`); the desired modulus of
          the output. Set to 0 for "as high as possible".

        OUTPUT:

        `F(x)` in the following sense: we return a profinite integer ``y``
        satisfying:

        1. any profinite integer that ``x`` represents is mapped by `F` to a
           profinite integer that is represented by ``y``.
        2. Let `m` be the maximal divisor of ``des_mod`` for which ``y`` can
           have modulus `m` under condition 1, if it exists. Then the modulus
           of ``y`` is a multiple of `m`.

           If ``des_mod`` is zero, such `m` might not exist. In this case, the
           modulus of ``y`` should be zero if this is possible under condition
           1. Else, no guarantees are made about the modulus of ``y``. An
           implementer of such a function should clearly specify its own
           precision guarantees in its documentation.

        Condition 1 above says that the output should be "correct". If we give
        ``a mod m`` as input and we receive ``b mod n`` as output, than the
        following must hold:

        .. MATH::

            F(a + m \cdot \hat{\ZZ}) \subset b + n \cdot \hat{\ZZ}.

        Condition 2 says that the output should have a "high" precision. One
        could always output ``0 mod 1`` to satisfy condition 1, since this
        represents the whole ring `\hat{\ZZ}`. But a user wants more precision
        of course. When a user specifies a non-zero ``des_mod``, the precision
        should be as high as possible modulo ``des_mod``. This is the main
        content of condition 2. We allow a modulus that is *a multiple of* `m`,
        not just `m`: returning a multiple of `m` gives strictly more
        information. Note that this always includes zero. So an exact answer may
        always be returned as well.

        For non-zero ``des_mod``, a maximal `m` as mentioned always exists,
        since ``des_mod`` has finitely many divisors and 1 is a candidate
        divisor that always satisfies condition 1.

        Note that if ``des_mod`` is zero, there may still be a maximal `m`: the
        input is a bounded subset of `\hat{\ZZ}`, so in many cases the output
        will also be some bounded subset. In such a case, a maximal `m` could
        exist.

        EXAMPLES:

        This constant-zero map only returns the result modulo the desired
        modulus, although we know the result exactly::

            sage: class ConstantZero(ProfiniteFunction):
            ....:     def __call__(self, x, des_mod=0):
            ....:         return Zhat(0, des_mod)
            sage: zero = ConstantZero()
            sage: zero(Zhat(5,20), 0)
            0 mod 0
            sage: zero(Zhat(5,20), 30)
            0 mod 30

        The constant zero function above is also an example of a
        ProfiniteFunction that returns an exact answer (output modulus 0) when
        ``des_mod`` is zero.

        We may also always return an answer with the maximal modulus::

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
        answer with strictly more precision: a multiple of 5.

        For a good "real world" example, take a look at the
        :class:`profinite Fibonacci function <sage.rings.adeles.profinite_functions.ProfiniteFibonacci>`.

        For a (theoretical) example of a function for which we cannot give
        precision guarantees, consider the prime decomposition
        `\hat{\ZZ} = \prod_p \ZZ_p` and the function
        `G: \hat{\ZZ} \to \hat{\ZZ}` which is given on the `p`-th coordinate as
        follows. For `p \neq 2`, `G` is the `p`-adic identity function on
        `\ZZ_p`, while on the `2`-adics `G` is the constant 1 map.
        Suppose the user asks for `G(2 \mod 6)` with desired modulus 0. Then for
        any positive integer `n`, we know the answer modulo `2^n` (namely 1).
        Hence the output moduli satisfying condition 1 are unbounded. But we do
        *not* know the answer exactly: we don't know the answer modulo 5.

        .. NOTE::

            The purpose of ``des_mod`` is to allow implementations of
            ProfiniteFunctions to speed up their computation, knowing that the
            user is only interested in the result modulo ``des_mod``.

            Hence if function evaluation is very slow, the user is adviced to
            set ``des_mod`` to a non-zero value.
        """
        pass



class ProfiniteFibonacci(ProfiniteFunction):
    r"""
    The Fibonacci function from and to the profinite integers

    The Fibonacci function `F` is defined on the integers by `F(0)=0`, `F(1)=1`
    and `F(n) = F(n-1) + F(n-2)` for all integers  `n \in \ZZ`.
    There exists a unique continuous extension `F: \hat{\ZZ} \to \hat{\ZZ}` of
    the Fibonacci function to the profinite integers.
    This is an implementation of that function.

    .. automethod:: __call__
    """
    def __call__(self, x, des_mod=0):
        r"""
        Return the ``x``-th profinite Fibonacci number

        We refer to :meth:`ProfiniteFunction.__call__
        <sage.rings.adeles.profinite_functions.ProfiniteFunction.__call__>`
        for a description of the input and output.

        The modulus of the output will always be a divisor of ``des_mod``.
        So we never compute up to a higher precision than the user asks for.

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

        Let `m` and `n` be the modulus of the input (``x``) and the output 
        respectively. Then clearly we have the upper bound
        
        .. MATH::

            n \leq gcd(F(m)-F(0), F(m+1)-F(1)).

        Using the Fibonacci-recurrence relation and the fact that the Fibonacci
        function is continous (as a function `F: \hat{\ZZ} \to \hat{\ZZ}`) one
        can show that this is an equality. This enables us to directly compute
        the output modulus.

        Hence for ``x = v mod m``, it suffices to compute F(m), F(m+1) and F(v).
        This is done using the formula

        .. MATH::

            \left(\begin{matrix} 1&1\\1&0 \end{matrix}\right)^n
            =
            \left(\begin{matrix} F(n+1)&F(n) \\F(n)&F(n-1) \end{matrix}\right)

        The matrix exponentiation is done using a square-and-multiply method in
        the ring of matrices over `\ZZ/\texttt{des_mod}\ZZ`.

        .. WARNING::

            Not specifying ``des_mod``, or setting it to 0, requires the exact
            computation of F(x.modulus), which may be huge
            (`F(n) \approx 1.6^n`). Hence the computation with bounded desired
            modulus is significantly faster than for unbounded.
        
        .. SEEALSO::
            
            The class
            :class:`ProfiniteGraph <sage.rings.adeles.profinite_graph.ProfiniteGraph>`
            for a visualization of the graph of this function.
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
