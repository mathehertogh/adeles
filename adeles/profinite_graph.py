r"""
Profinite Graph

Let `F: \hat{\ZZ} \to \hat{\ZZ}` be a function and let `P` be a
:class:`profinite_function.ProfiniteFunction` implementing `F`. This file
implements interactive profinite graphs in the class :class:`ProfiniteGraph`:
a user can view arbitrarily high precision approximations of the graph
`\{(x, F(x)) \mid x \in \hat{\ZZ}\}` of `F` by zooming in and out.

We use the class :class:`profinite_integer.ProfiniteIntegers` as our
implementation of `\hat{\ZZ}` and by a profinite integer we shall mean an
instance of ``Zhat``, i.e. ``ProfiniteIntegers(QQ)``. Also by ``3 mod 6`` we
shall mean the profinite integer ``Zhat(3, 6)``.

In this file we use the character '!' to mean factorial (e.g. `4!` is `24`).

Factorial digits
--------------------------------

Every `\alpha \in \hat{\ZZ}` is uniquely determined by its *factorial digits*,
which are the unique integers `d_1, d_2, d_3, d_4, ... \in \ZZ` satisfying
`0 \leq d_k \leq k` and `\alpha \equiv \sum_{i=1}^k d_i \cdot i! \mod
(k+1)!\hat{\ZZ}` for every `k \in \ZZ_{\geq 1}`. Knowing the first `k` factorial
digits of `\alpha` corresponds to knowing `\alpha` modulo `(k+1)!`. These
factorial digits are implemented in ``Zhat``::

    sage: a = Zhat([1, 2, 1, 1]); a
    35 mod 120
    sage: print(a.str(style='factorial'))
    1*1! + 2*2! + 1*3! + 1*4! + O(5!)
    sage: a.factorial_digits()
    [1, 2, 1, 1]

::

    sage: b = Zhat([1, 2, 1, 1, 0, 0]); b
    35 mod 5040
    sage: print(b.str(style='factorial'))
    1*1! + 2*2! + 1*3! + 1*4! + O(7!)
    sage: b.factorial_digits()
    [1, 2, 1, 1, 0, 0]

The visualization function
----------------------------------------------------

We define the *visualization function* to be the map `\phi: \hat{\ZZ} \to
[0, 1]` (where `[0, 1] \subset \RR` denotes the unit interval) defined by
`\phi(\alpha) = \sum_{i=1}^\infty d_i / (i+1)!` for `d_i` the factorial digits
of `\alpha`.

It is illustrative to verify the equalities `\phi(1+2\hat{\ZZ}) = [1/2, 1]` and
`\phi(2+3\hat{\ZZ}) = [1/6, 1/3] \cup [5/6, 1]`. Can you determine
`\phi(8+24\hat{\ZZ})`?

The visualization function enables us to visualize the ring `\hat{\ZZ}` as the
unit interval and `\hat{\ZZ} \times \hat{\ZZ}` as the unit square. It is
implemented in ``Zhat``::

    sage: Zhat(1, 2).visual()
    (1/2, 1)
    sage: Zhat(2, 3).visual()
    (1/6, 1)
    sage: a.visual()
    (53/60, 107/120)
    sage: b.visual()
    (53/60, 4453/5040)

Graphing a profinite function
----------------------------------------------------------

Let `k` be a positive integer which we call our *precision*. Then we can draw
the graph of a function `F: \hat{\ZZ} \to \hat{\ZZ}` at precision `k` as
follows. For each profinite integer `x` of modulus `k!`, we compute the image of
the represented subset of `x` under `F` and we take all profinite integers `y`
of modulus `k!` whose represented subset intersects this image. Each of
these profinite integers `x` and `y` are mapped by the visualization function to
a closed interval. Hence we can draw the points `(x, y)` as squares in the unit
interval. This is precisely what the class :class:`ProfiniteGraph` does. See
the documentation of that class for examples.

REFERENCES:

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

[Len2005] Hendrik Lenstra, Profinite Fibonacci numbers, Niew Archief voor
Wiskunde, 5/6(4):297-300, december 2005.
http://www.nieuwarchief.nl/serie5/pdf/naw5-2005-06-4-297.pdf

This implementation is based on and part of [Her2021]. For details, see Chapter
7 of [Her2021]. The idea of these kinds of graphs is taken from [Len2005].

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

try:
    from tkinter import Tk, Canvas, TclError
    from tkinter.font import Font
except ModuleNotFoundError:
    tkinter_available = False
else:
    tkinter_available = True
from sage.functions.other import factorial
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ

from .profinite_integer import Zhat



class ProfiniteGraph:
    r"""
    Interactive graph of a function from and to profinite integers
    
    Every :class:`profinite_function.ProfiniteFunction` can be graphed, as well
    as callables that behave like such profinite functions (see below for
    examples).

    The graph allows the user to zoom in (and back out) to arbitrary high
    precisions. Left-click on an area in the graph to zoom into that area.
    Right-click anywhere on the graph to zoom out.

    The profinite integers are drawn based on the
    :meth:`profinite_integer.ProfiniteInteger.visual` method of profinite
    integers. So `\hat{\ZZ} \times \hat{\ZZ}` is identified with the unit square
    `[0,1] \times [0,1]` and a pair of profinite integers (``x mod m``, ``y mod
    n``) is identified with a rectangle `[a,b] \times [c,d]` in the unit square.

    The idea of this kind of graph is taken from the paper "Profinite Fibonacci
    numbers" by Hendrik Lenstra, cf. [Len2005]. Some default settings such as
    the colors to draw in are taken such that the resulting graph looks like the
    picture in that paper.

    .. WARNING::

        This class uses the standard Python interface package `tkinter
        <https://docs.python.org/3/library/tkinter.html>`_ for its graphical
        interface. Make sure that tkinter is installed (which should already be
        the case on Windows and most Unix systems) and that the command ``import
        tkinter`` works.

    EXAMPLES:

    We plot the graph of the profinite Fibonacci function::

        sage: g = ProfiniteGraph(ProfiniteFibonacci())
        sage: g.set_window_sizes(720, 720)
        sage: g.set_title("Graph of the profinite Fibonacci function")
        sage: g.plot() # optional - tkinter

    This produces the following image:

    .. image:: Fibonacci_graph_0mod1x0mod1.png

    What do we see here? Lets first look only at the big orange blocks. Above
    ``3 mod 6`` on the `x`-axis, there are two orange blocks, which on the
    `y`-axis represent ``2 mod 6`` and ``4 mod 6``. This means that the image of
    `3 + 6 \cdot \hat{\ZZ}` under the Fibonacci function lies in
    `(2 + 6 \cdot \hat{\ZZ}) \cup (4 + 6 \cdot \hat{\ZZ})`. It also means that
    both ``2 mod 6`` and ``4 mod 6`` are "hit" by ``3 mod 6``:
    `F(3 + 6 \cdot \hat{\ZZ}) \cap (2 + 6 \cdot \hat{\ZZ}) \neq \emptyset` and
    `F(3 + 6 \cdot \hat{\ZZ}) \cap (4 + 6 \cdot \hat{\ZZ}) \neq \emptyset`.

    The orange blocks represent the approximation of the graph of the profinite
    Fibonacci function of *precision* 3: it computes the images of the profinite
    integers of modulus `3!`. In pink, the approximation of precision 4 is
    drawn. Hence the pink rectangles all represents subsets of `\hat{\ZZ} \times
    \hat{\ZZ}` of the form `(x+4!\hat{\ZZ}) \times (y+4!\hat{\ZZ})`.

    By left clicking we zoom in to the area `(3+6\hat{\ZZ}) \times
    (2+6\hat{\ZZ})`, i.e. ``(3 mod 6) x (2 mod 6)``:

    .. image:: Fibonacci_graph_3mod6x2mod6.png

    Now we see that within `3+6\hat{\ZZ}`, the open subsets `3+24\hat{\ZZ}` and
    `21+24\hat{\ZZ}` are actually mapped to `2+24\hat{\ZZ}`. Apparently, the
    rest of `3+6\hat{\ZZ}` (i.e. `9+24\hat{\ZZ}` and `15+24\hat{\ZZ}`) is mapped
    into `4+6\hat{\ZZ}`. One could see this by right-clicking anywhere on the
    graph to zoom out and than zooming in to the area ``(3 mod 6) x (4 mod 6)``.

    Now within the pink squares, we see brown squares. This is the precision 5
    approximation. From this picture one can already see that ``3 mod 120`` is
    mapped to ``2 mod 120`` and that ``27 mod 120`` is mapped to ``98 mod 120``.
    But a user can also zoom in again to see this more clearly of course.

    Any profinite function (cf. :class:`profinite_function.ProfiniteFunction`)
    can be graphed. We can even graph any callable that behaves like a
    profinite function. Here is an example of how to produce a graph of the
    square function `\hat{\ZZ} \to \hat{\ZZ}, x \mapsto x^2`::

        sage: graph = ProfiniteGraph(lambda x, des_mod: x*x)
        sage: graph.set_title("Square function")
        sage: graph.plot()  # optional - tkinter
    """


    class View:
        r"""
        This class stores the state of the "view" of a ProfiniteGraph.

        A ProfiniteGraph allows the user to change the view by zooming in and
        out to open subsets of `\hat{\ZZ} \times \hat{\ZZ}` of the form
        `(x + p! \cdot \hat{\ZZ}) \times (y + p! \cdot \hat{\ZZ})`, i.e.
        ``(x mod p!) x (y mod p!)``.
        Here x, y and p are integers. We call p above the "precision".
        """
        def __init__(self):
            r"""
            Construct a fully zoomed out view

            Initial view is ``(0 mod 1!) x (0 mod 1!)``, i.e. the whole
            `\hat{\ZZ} \times \hat{\ZZ}`.

            TESTS::

                sage: view = ProfiniteGraph.View()
            """
            self.prec = ZZ(1)
            self.x = ZZ(0)
            self.y = ZZ(0)

        def __repr__(self):
            """
            Return a string representation of this view

            EXAMPLES::

                sage: view = ProfiniteGraph.View()
                sage: view.prec += 1
                sage: view.x = ZZ(1)
                sage: view
                (1 mod 2!) x (0 mod 2!)
            """
            x, y, p = self.x, self.y, self.prec
            return "({} mod {}!) x ({} mod {}!)".format(x, p, y, p)



    def __init__(self, function):
        r"""
        Create the graph of the ProfiniteFunction ``function``

        INPUT:

        - ``function`` -- a profinite function (cf.
          :class:`profinite_function.ProfiniteFunction`), or a callable that
          behaves like a profinite function.


        EXAMPLES::

            sage: fibonacci = ProfiniteFibonacci()
            sage: graph = ProfiniteGraph(fibonacci)

        An example where we pass in a callable behaving like a profinite
        function, which is not an instance of
        :class:`profinite_function.ProfiniteFunction`. ::

            sage: cube_graph = ProfiniteGraph(lambda x, des_mod: x*x*x)

        This creates the graph of the cube map `\hat{\ZZ} \to \hat{\ZZ}, x
        \mapsto x^3`.

        .. NOTE::

            This *creates* the graph, but does not *display* is yet. Call
            :meth:`plot` to actually view the graph.

        TESTS::

            sage: ProfiniteGraph(97)
            Traceback (most recent call last):
            ...
            ValueError: function is not a profinite function (and does not behave like one)
        """
        if not tkinter_available:
            raise ModuleNotFoundError("tkinter is not available")
        try:
            if function(Zhat(5, 24), 6) not in Zhat:
               raise ValueError("function is not a profinite function (and does not behave like one)")
        except:
            raise ValueError("function is not a profinite function (and does not behave like one)")
        self.function = function
        self.view = ProfiniteGraph.View()
        self.canvas = None # the canvas on which we draw everything
        self.highlight_rect = None # yellowish rectangle that follows the mouse
        self._init_settings()

    def plot(self):
        """
        Open the interactive graph in a new window

        EXAMPLES::

            sage: cube_graph = ProfiniteGraph(lambda x, des_mod: x*x*x)
            sage: cube_graph.plot()  # optional - tkinter

        .. NOTE::

            Control flow of the program is handed over to the drawing library
            used (``Tkinter``) until the user closes the graph-window.
        """
        self._create_window()
        self._draw()
        self.window.mainloop()

    def set_title(self, title):
        """
        Set the title of the graph-window to ``title``

        INPUT:
        
        - ``title`` -- string; the new title

        EXAMPLES::

            sage: graph = ProfiniteGraph(lambda x, des_mod: x*x)
            sage: graph.set_title("Graph of square-function")

        TESTS::

            sage: graph = ProfiniteGraph(lambda x, des_mod: x*x)
            sage: graph.set_title(79)
            Traceback (most recent call last):
            ...
            TypeError: title should be a string
        """
        if not isinstance(title, str):
            raise TypeError("title should be a string")
        self.title = title

    def set_colors(self, approx=None, highlight=None, identity_line=None, axis=None):
        """
        Set the colors of the graph

        INPUT:

        - ``approx`` -- non-empty list of strings (optional, default is keep the
          current setting); colors of the different approximations
        - ``highlight`` -- string (optional, default is keep the current
          setting); color of the transparent box that follows the mouse
        - ``identity_line`` -- string (optional, default is keep the current
          setting); color of the identity line, see also
          :meth:`set_identity_line`
        - ``axis`` -- string (optional, default is keep the current setting);
          color of the axis and the corresponding coordinate labels (e.g. "1/6",
          "2/3")
        
        Strings to denote a color should be formatted in one of the following
        ways:

        - A hexadecimal string specifying the red, green and blue portions of
          the color, in that order: "#rrggbb". Examples: "#00ff00" is green,
          "#00ffff" is cyan and "#ff0077" is pink.
        - A locally defined standard color name. The colors "white", "black",
          "red", "green", "blue", "cyan", "yellow", and "magenta" are always
          available. Usually more colors work, like "brown", but this depends
          on your local installation.

        EXAMPLES::

            sage: graph = ProfiniteGraph(lambda x, des_mod: Zhat(0, 0))
            sage: graph.set_colors(approx=["blue", "green", "yellow", "red"])
            sage: graph.set_colors(highlight="#efefef", identity_line="#0033aa")
            sage: graph.set_colors(axis="black")
        """
        if approx is not None:
            if not isinstance(approx, list):
                raise TypeError("approx should be a list")
            if len(approx) == 0:
                raise ValueError("approx should be non-empty")
            for color in approx:
                if not self._is_valid_color(color):
                    raise ValueError("{} does not define a color".format(color))
            self.approx_colors = approx
        if highlight is not None:
            if not self._is_valid_color(highlight):
                raise ValueError("{} does not define a color".format(highlight))
            self.highlight_color = highlight
        if identity_line is not None:
            if not self._is_valid_color(identity_line):
                raise ValueError("{} does not define a color".format(identity_line))
            self.identity_line_color = identity_line
        if axis is not None:
            if not self._is_valid_color(axis):
                raise ValueError("{} does not define a color".format(axis))
            self.axis_color = axis

    def set_identity_line(self, draw):
        r"""
        Set whether or not to draw the identity line

        By the identity line we mean de graph of the identity function on
        `\hat{\ZZ}`, i.e. "`x = y`". We draw this as an actual line (with zero
        surface area). This can be useful to find fixed points of the function
        you are graphing. This is for example done in the paper [Len2005] for
        the graph of the Fibonacci function.

        INPUT:

        - `draw` -- boolean; whether or not to draw the identity line

        EXAMPLES::

            sage: graph = ProfiniteGraph(lambda x, des_mod: x*x)
            sage: graph.set_identity_line(draw=False)

        .. SEEALSO::

            To set the color of the identity line, see :meth:`set_colors`.
        """
        if draw is not None:
            if not isinstance(draw, bool):
                raise TypeError("parameter `draw` should be a boolean")
            self.identity_line = draw

    def set_minimal_draw_area(self, min_area):
        r"""
        Set the minimal draw area of this graph to ``min_area``

        We only draw approximations up to a precision such that we can still
        see the rectangles drawn for pairs `(x, y)` in `\hat{\ZZ} \times
        \hat{\ZZ}`. The minimal area in pixels a rectangle representing such a
        point `(x, y)` must have to be drawn is called the "minimal draw area".
        This method sets this minimal draw area.

        INPUT:

        - ``min_area`` -- integer; the new minimal draw area, in pixels

        EXAMPLES::

            sage: graph = ProfiniteGraph(lambda x, des_mod: x)
            sage: graph.set_minimal_draw_area(25)
        """
        if min_area not in ZZ:
            raise TypeError("{} is not an integer".format(min_area))
        self.min_draw_area = min_area

    def set_minimal_axis_distance(self, min_axis_distance):
        """
        Set the minimal axis distance

        We usualy draw axes at two different precisions. But if you zoom in very
        far, the axes with the highest precision start taking up (almost) the
        whole graph. Hence if two such consecutive axes are less than a certain
        number of pixels apart, we don't draw them any more; we only draw the
        lowest precision axis. This minimal number of pixels that these axes
        should be apart we call the "minimal axis distance".

        INPUT:

        - ``min_axis_distance`` -- integer; the new minimal axis distance

        EXAMPLES::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph.set_minimal_axis_distance(10)

        .. SEEALSO::

            To set the color of the axis, see :meth:`set_colors`.
        """
        if min_axis_distance not in ZZ:
            raise TypeError("{} is not an integer".format(min_axis_distance))
        self.min_axis_distance = min_axis_distance

    def set_window_sizes(self, width=None, height=None, border=None, footer=None):
        """
        Set the sizes in pixels of the graph window

        We define sizes in pixels of parts of the window, as depicted in
        ASCII-art below. ::

            -----------------------------------------------------------
            |            ^                                            |
            |            |border                1 mod 2               |
            |         0  v              1/2                 1         |
            |        1-------------------|-------------------1        |
            | border  |     ^                               | border  |
            |<------->|     |            width              |<------->|
            |         |<----------------------------------->|         |
            |         |     |                               | 1 mod 2 |
            |         |     |                               |         |
            |         |     |                               |         |
            |      1/2|     |                               |1/2      |
            |         |     |                               |         |
            |         |     |height                         |         |
            | 0 mod 2 |     |                               | 0 mod 2 |
            |         |     |                               |         |
            |         |     v                               |         |
            |        0-------------------|-------------------0        |
            |         0  ^              1/2                 1         |
            |            |border                1 mod 2               |
            |   _   _   _v  _   _   _   _   _   _   _   _   _   _   _ |
            |            ^                                            |
            |            |footer    currenty viewing: ...             |
            |            v                                            |
            -----------------------------------------------------------

        Note that width and height are *not* the width and height of the full
        window.
        
        For the best drawing results, make sure k! divides width and height,
        for a "high" value of k (at least 4!, but strive for 5!).

        EXAMPLES::

            sage: graph = ProfiniteGraph(lambda x, des_mod: Zhat(2)*x)
            sage: graph.set_window_sizes(720, 720, 30, 15)
        """
        if width is not None:
            if width not in ZZ:
                raise TypeError("{} is not an integer".format(width))
            self.width = width
        if height is not None:
            if height not in ZZ:
                raise TypeError("{} is not an integer".format(height))
            self.height = height
        if border is not None:
            if border not in ZZ:
                raise TypeError("{} is not an integer".format(border))
            self.border = border
        if footer is not None:
            if footer not in ZZ:
                raise TypeError("{} is not an integer".format(footer))
            self.footer = footer



    def _init_settings(self):
        """
        Initialize the settings of the graph

        All of these settings have "public" setters to be used by the user.
        See those functions for an explanation of the variables defined below.

        For example: the meaning of ``self.min_draw_area`` is explained in the
        documentation of the method set_minimal_draw_area().

        TESTS::

            sage: graph = ProfiniteGraph(lambda x, des_mod: x-Zhat(2))
            sage: graph._init_settings()
        """
        self.title = "Interactive profinite graph - left-click to zoom in, right-click to zoom out"
        self.approx_colors = ["orange", "pink", "sienna3"] # sienna3 is light brown
        self.highlight_color = "yellow"
        self.identity_line = False
        self.identity_line_color = "green"
        self.min_draw_area = ZZ(1)
        self.min_axis_distance = ZZ(8)
        self.axis_color = "cornflower blue"
        self.width = ZZ(600)
        self.height = ZZ(600)
        self.border = ZZ(42)
        self.footer = ZZ(20)

    def _create_window(self):
        """
        Creates the ``Tk``-window with canvas and sets up event bindings

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._create_window()
        """
        self.window = Tk()
        self.window.title(self.title)
        window_height = self.height + 2*self.border + self.footer
        window_width = self.width + 2*self.border
        self.window.geometry("{}x{}+{}+{}".format(window_width, window_height, 0, 0))
        self.canvas = Canvas(self.window, width=window_width, height=window_height)
        self.canvas.configure(background="white")
        self.canvas.pack()
        self.canvas.bind("<Button-1>", self._zoom_in) # left mouse click
        self.canvas.bind("<Button-3>", self._zoom_out) # right mouse click
        self.canvas.bind("<Motion>", self._update_highlight) # mouse movement



    def _zoom_in(self, event):
        """
        Zoom in to the region that was clicked durent event ``event``

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph.view
            (0 mod 1!) x (0 mod 1!)
            sage: import tkinter
            sage: ev = tkinter.Event()
            sage: ev.x = 2*graph.border
            sage: ev.y = 2*graph.border
            sage: graph._zoom_in(ev)
            sage: graph.view
            (0 mod 2!) x (1 mod 2!)
        """
        x, y = self._find_region(event)
        if x is None:
            return
        self.view.x, self.view.y = x, y
        self.view.prec += 1
        self._draw()

    def _zoom_out(self, event):
        """
        Zoom out one precision level

        Do nothing if we are already maximally zoomed out, that is the zoom is 
        (0 mod 1!)x(0 mod 1!).

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph.view
            (0 mod 1!) x (0 mod 1!)
            sage: graph._zoom_out(None)
            sage: graph.view
            (0 mod 1!) x (0 mod 1!)
            sage: import tkinter
            sage: ev = tkinter.Event()
            sage: ev.x, ev.y = 2*graph.border, 2*graph.border
            sage: graph._zoom_in(ev)
            sage: graph.view
            (0 mod 2!) x (1 mod 2!)
            sage: graph._zoom_out(ev)
            sage: graph.view
            (0 mod 1!) x (0 mod 1!)
        """
        if self.view.prec <= 1:
            return
        self.view.prec -= 1
        self.view.x = self.view.x % factorial(self.view.prec)
        self.view.y = self.view.y % factorial(self.view.prec)
        self._draw()

    def _update_highlight(self, event):
        """
        Reposition the highlight-rect to follow the mouse
        
        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph.canvas.coords(graph.highlight_rect)
            [0.0, 0.0, 1.0, 1.0]
            sage: import tkinter
            sage: ev = tkinter.Event()
            sage: ev.x, ev.y = graph.border, graph.border
            sage: graph._update_highlight(ev)
            sage: graph.canvas.coords(graph.highlight_rect)
            [42.0, 42.0, 342.0, 342.0]
            sage: ev.x, ev.y = graph.width, graph.height
            sage: graph._update_highlight(ev)
            sage: graph.canvas.coords(graph.highlight_rect)
            [342.0, 342.0, 642.0, 642.0]
            sage: ev.x, ev.y = 0, 0
            sage: graph._update_highlight(ev)
            sage: graph.canvas.coords(graph.highlight_rect)
            [-1.0, -1.0, -1.0, -1.0]
        """
        v, w = self._find_region(event)
        if v is None:
            # Mouse is outside of the actual graph (so on the border or footer),
            # so we make the highlight-rectangle invisible.
            self.canvas.coords(self.highlight_rect, -1, -1, -1, -1)
            return
        a = Zhat(v, factorial(self.view.prec+1))
        b = Zhat(w, factorial(self.view.prec+1))
        x1, x2 = self._pixel_coordinates(a, is_x=True)
        y1, y2 = self._pixel_coordinates(b, is_x=False)
        x1, y1 = self._in_canvas_coordinates(x1, y1)
        x2, y2 = self._in_canvas_coordinates(x2, y2)
        self.canvas.coords(self.highlight_rect, x1, y1, x2, y2)



    def _is_valid_color(self, color):
        """
        Check whether or not ``color`` defines a valid color for Tkinter

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._is_valid_color("blabla")
            False
            sage: graph._is_valid_color("yellow")
            True
            sage: graph._is_valid_color("#08e1ab")
            True
        """
        try:
            test = self.canvas.create_text(0, 0, fill=color)
            self.canvas.delete(test)
        except TclError:
            return False
        return True

    def _depth(self):
        """
        Return the "depth": the number of approximations we draw

        All approximations are drawn for which the individual blocks have an
        area of at least ``min_draw_area``.

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._depth()
            3
            sage: graph.set_window_sizes(100, 100)
            sage: graph._depth()
            2
            sage: graph.set_minimal_draw_area(50)
            sage: graph._depth()
            1
        """
        depth = ZZ(1)
        while self._area_on_screen(self._prec() + depth) >= self.min_draw_area:
            depth += ZZ(1)
        return depth

    def _prec(self):
        """
        Return the precision of the lowest precision approximation we draw

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._prec()
            3
            sage: graph.view.prec += 1
            sage: graph._prec()
            3
            sage: graph.view.prec += 1
            sage: graph._prec()
            4
        """
        return max(ZZ(3), self.view.prec+1)

    def _color(self, p):
        """
        Return the color of the approximation of precision ``p``

        TESTS::
        
            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._color(1) in graph.approx_colors
            True
        """
        return self.approx_colors[p % len(self.approx_colors)]



    def _pixel_coordinates(self, a, is_x):
        """
        Return the pixel-range in which ``a`` lies

        INPUT:

        - a -- profinite integer
        - is_x -- boolean; ``True`` for `x`-coordinates, ``False`` for
                  `y`-coordinates
        
        OUTPUT:

        A pair of integers ``(s, t)``, such that ``a`` is drawn in the current
        view at pixels between ``s`` and ``t`` (inclusive).

        If ``a`` is (partly) outside of the view, return ``(None, None)``.

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._pixel_coordinates(Zhat(2, 6), True)
            (100, 200)
            sage: graph.view.prec += 1
            sage: graph._pixel_coordinates(Zhat(2, 6), True)
            (200, 400)
            sage: graph._pixel_coordinates(Zhat(5, 6), True)
            (None, None)
            sage: graph._pixel_coordinates(Zhat(4, 72), True)
            (400, 450)
        """
        l, r = a.visual()
        value = self.view.x if is_x else self.view.y
        view_l, view_r = Zhat(value, factorial(self.view.prec)).visual()
        if not (view_l <= l <= r <= view_r):
            return None, None
        view_size = view_r - view_l
        length = self.width if is_x else self.height
        left = int(length * (l - view_l) / (view_r - view_l))
        right = int(length * (r - view_l) / (view_r - view_l))
        return left, right

    def _in_canvas_coordinates(self, x, y):
        """
        Convert the graph-coordinates ``(x,y)`` to canvas coordinates

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._in_canvas_coordinates(100, 100)
            (142, 542)
            sage: graph._in_canvas_coordinates(graph.width, graph.height)
            (642, 42)
        """
        return self.border+x, self.border+self.height-y

    def _area_on_screen(self, prec):
        """
        Return the area in pixels of an element of precision ``prec``
        
        INPUT:

        - ``prec`` -- positive integer; the precision
        
        OUTPUT:

        The area in pixels of the block that will be drawn for an open subset of
        the form ``(a mod prec!) x (b mod prec!)`` in the current view.

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._area_on_screen(3)
            10000
            sage: graph._area_on_screen(4)
            625
            sage: graph._area_on_screen(5)
            25
        """
        x = Zhat(self.view.x, factorial(prec))
        y = Zhat(self.view.y, factorial(prec))
        left, right = self._pixel_coordinates(x, is_x=True)
        bottom, top = self._pixel_coordinates(y, is_x=False)
        return (right - left) * (top - bottom)

    def _find_region(self, event):
        """
        Compute the open subset in which ``(event.x, event.y)`` lies
        
        INPUT:
        
        - ``event`` -- ``Tkinter.Event``
        
        OUTPUT:

        A pair of integers ``(v, w)`` such that the event was inside the region
        (``v`` mod ``p``!) x (``w`` mod ``p``!), where ``p=view.prec+1``.

        If the event was outside the actual graph (so in the border or the
        footer), then ``(None, None)`` is returned.

        ALGORITHM:

        For each x-coordinate in the current view, we calculate the pixel-range
        where it is/should be drawn and check whether ``x`` lies within this
        range.
        We do the same for the y-coordinate.

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: import tkinter
            sage: ev = tkinter.Event()
            sage: ev.x, ev.y = graph.border, graph.border
            sage: graph._find_region(ev)
            (0, 1)
            sage: ev.x = 0
            sage: graph._find_region(ev)
            (None, None)
        """
        v, w = None, None
        for i in range(self.view.prec+1):
            v_cand = self.view.x + i*factorial(self.view.prec)
            a = Zhat(v_cand, factorial(self.view.prec+1))
            s, t = self._pixel_coordinates(a, is_x=True)
            x = event.x - self.border
            if s <= x <= t:
                v = v_cand
                break
        for i in range(self.view.prec+1):
            w_cand = self.view.y + i*factorial(self.view.prec)
            b = Zhat(w_cand, factorial(self.view.prec+1))
            s, t = self._pixel_coordinates(b, is_x=False)
            y = self.height + self.border - event.y
            if s <= y <= t:
                w = w_cand
                break
        if v is None or w is None:
            return None, None
        return v, w

    def _draw_pair(self, x, y, color):
        """
        Draw the point (``x``, ``y``) in color ``color``

        If (the rectangle representing) the point (``x``, ``y``) lies (partly)
        outside of the view, nothing is drawn.

        INPUT:

        - ``x``, ``y`` -- profinite integers
        - ``color`` -- string; indicating the color to draw with

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: len(graph.canvas.find_all())
            274
            sage: graph._draw_pair(Zhat(1,24), Zhat(2,6), "red")
            sage: len(graph.canvas.find_all())
            275
            sage: graph.view.x, graph.view.prec = 1, 2
            sage: graph._draw_pair(Zhat(0,6), Zhat(2,6), "black")
            sage: len(graph.canvas.find_all())
            275
        """
        # x, y: profinite integers
        x1, x2 = self._pixel_coordinates(x, is_x=True)
        y1, y2 = self._pixel_coordinates(y, is_x=False)
        if x1 is None or y1 is None:
            return
        x1, y1 = self._in_canvas_coordinates(x1, y1)
        x2, y2 = self._in_canvas_coordinates(x2, y2)
        self.canvas.create_rectangle(x1, y1, x2, y2, fill=color, width=0)

    def _compute_image(self, x, prec):
        r"""
        Compute the image of ``x mod prec!`` under ``self.function``

        INPUT:

        - ``x`` -- integer; value of the element we want to know the image of
        - ``prec`` -- integer; precision of the input and output

        OUTPUT:

        A minimal set `S` of integers ``y`` such that the image of each profinite
        integer in ``(x mod prec!)`` under ``self.function`` is contained in
        ``(y mod prec!)`` for some ``y`` in `S`.

        ALGORITHM:
        
        Evaluates ``self.function`` in ``z mod p!`` for increasing values of
        ``p``, until the precision of the output is at least ``prec``.
        For ``z`` we use all preimages of ``x`` under the canonical projection
        `\ZZ/p\ZZ \to \ZZ/``prec``\ZZ`.

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._compute_image(1, 3)
            {1, 5}
            sage: graph._compute_image(4, 3)
            {1, 3, 5}
            sage: graph._compute_image(0, 4)
            {0}
            sage: graph._compute_image(7, 5)
            {13}
        """
        output_modulus = factorial(prec)
        input_prec = prec
        input_modulus = factorial(input_prec)
        result = set()
        while x < input_modulus:
            image = self.function(Zhat(x, input_modulus), output_modulus)
            if image.modulus() < output_modulus and image.modulus() != ZZ(0):
                input_prec += 1
                input_modulus = factorial(input_prec)
            else:
                # If ``function`` returns higher precision than requested, we
                # should make sure the returned value is at most output_modulus.
                value = image.value() % output_modulus
                result.add(value)
                x += output_modulus
        return result

    def _draw_approximations(self):
        """
        Draw all ``_depth()`` approximations within the current view

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: len(graph.canvas.find_all())
            274
            sage: graph._draw_approximations()
            sage: len(graph.canvas.find_all())
            433
        """
        view_modulus = factorial(self.view.prec)
        for p in range(self._prec(), self._prec()+self._depth()):
            for x in range(self.view.x, factorial(p), view_modulus):
                for y in self._compute_image(x, p):
                    if y % factorial(self.view.prec) == self.view.y:
                        pfx = Zhat(x, factorial(p))
                        pfy = Zhat(y, factorial(p))
                        self._draw_pair(pfx, pfy, self._color(p))

    def _draw_identity_line(self):
        """
        Draw the identity line (given by the equation `x=y`)

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: len(graph.canvas.find_all())
            274
            sage: graph._draw_identity_line()
            sage: len(graph.canvas.find_all())
            275
            sage: graph.set_identity_line(draw=False)
            sage: graph._draw_identity_line()
            sage: len(graph.canvas.find_all())
            275
            sage: graph.set_identity_line(draw=True)
            sage: graph._draw_identity_line()
            sage: len(graph.canvas.find_all())
            276
            sage: graph.view.x, graph.view.prec = 1, 2
            sage: graph._draw_identity_line()
            sage: len(graph.canvas.find_all())
            276
        """
        if self.identity_line and self.view.x == self.view.y:
            # If the identity line is within our view, then it is the diagonal.
            self.canvas.create_line(self.border, self.border+self.height,
                                    self.border+self.width, self.border,
                                    fill=self.identity_line_color)

    def _label_font(self):
        """
        Return the appropriately sized font based on the current view

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._label_font().actual()["size"]
            10
            sage: graph.view.prec += 5
            sage: graph._label_font().actual()["size"]
            8
            sage: graph.view.prec += 5
            sage: graph._label_font().actual()["size"]
            3
            sage: graph.view.prec += 5
            sage: graph._label_font().actual()["size"]
            1
            sage: graph.view.prec += 5
            sage: graph._label_font().actual()["size"]
            1
        """
        n_blocks = prod(range(self.view.prec+1, self._prec()+1))
        max_width = self.width / n_blocks
        p = self._prec()
        subset_label = " {} mod {}! ".format(factorial(p), p)
        coordinate_label = " {}/{}! ".format(factorial(p), p)
        size = ZZ(10)
        while size > 0:
            font = Font(family="Calibri", size=size)
            if (font.measure(subset_label) <= max_width
                    and font.measure(coordinate_label) <= max_width):
                return font
            size -= 1
        return Font(family="Calibri", size=1)

    def _draw_axis(self):
        """
        Draw the axis of the graph (all horizontal and diagonal lines)

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: len(graph.canvas.find_all())
            274
            sage: graph._draw_axis()
            sage: len(graph.canvas.find_all())
            338
            sage: graph.view.prec += 2
            sage: graph._draw_axis()
            sage: len(graph.canvas.find_all())
            390
            sage: graph.view.prec += 1
            sage: graph._draw_axis()
            sage: len(graph.canvas.find_all())
            464
            sage: graph.view.prec += 30
            sage: graph._draw_axis()
            sage: len(graph.canvas.find_all())
            536
        """
        line_width = ZZ(1)
        n_blocks = prod(range(self.view.prec+1, self._prec()+1+1))
        size_between_lines = self.width / n_blocks
        if size_between_lines >= self.min_axis_distance:
            precs = [self._prec()+1, self._prec()]
        else:
            precs = [self._prec()]
        for p in precs:
            n_blocks = prod(range(self.view.prec+1, p+1))
            for i in range(n_blocks+1):
                x = self.border + int(self.width * i / n_blocks)
                self.canvas.create_line(x, self.border, x, self.border+self.height,
                                        fill=self.axis_color, width=line_width)
            for i in range(n_blocks+1):
                y = self.border + int(self.height * i / n_blocks)
                self.canvas.create_line(self.border, y, self.border+self.width, y,
                                        fill=self.axis_color, width=line_width)
            line_width += 1

    def _axis_label_text(self, v):
        """
        Return a string representing the fraction v/p!, with p=self._prec()
        
        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._axis_label_text(3)
            '1/2'
            sage: graph.view.prec += 3
            sage: graph._axis_label_text(3)
            '1/40'
            sage: graph.view.prec += 3
            sage: graph._axis_label_text(3)
            '3/8!'
        """
        if self._prec() <= 6:
            return str(v / factorial(self._prec()))
        return "{}/{}!".format(v, self._prec())

    def _draw_axis_labels(self):
        """
        Draw the coordinate axis labels (e.g. 0, 1/6, 1/3, 1/2, 2/3, ...)
        
        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: len(graph.canvas.find_all())
            274
            sage: graph._draw_axis_labels()
            sage: len(graph.canvas.find_all())
            298
            sage: graph.view.prec += 1
            sage: graph._draw_axis_labels()
            sage: len(graph.canvas.find_all())
            310
            sage: graph.view.prec += 1
            sage: graph._draw_axis_labels()
            sage: len(graph.canvas.find_all())
            326
        """
        n_blocks = prod(range(self.view.prec+1, self._prec()+1))
        y_up = ZZ(2) * self.border / ZZ(3)
        y_down = self.border + self.height + self.border / ZZ(3)
        font = self._label_font()
        for i in range(n_blocks+1):
            x = self.border + int(self.width * i / n_blocks)
            text = self._axis_label_text(self.view.x*n_blocks + i)
            if i != 0:
                self.canvas.create_text(x, y_up, text=text, fill=self.axis_color, font=font)
            if i != n_blocks:
                self.canvas.create_text(x, y_down, text=text, fill=self.axis_color, font=font)
        x_left = ZZ(2) * self.border / ZZ(3)
        x_right = self.border + self.width + self.border / ZZ(3)
        for i in range(n_blocks+1):
            y = self.border + int(self.height * i / n_blocks)
            text = self._axis_label_text((self.view.y+1)*n_blocks - i)
            if i != n_blocks:
                self.canvas.create_text(x_left, y, text=text, fill=self.axis_color, font=font, angle=90)
            if i != 0:
                self.canvas.create_text(x_right, y, text=text, fill=self.axis_color, font=font, angle=-90)

    def _subset_label_text(self, v):
        """
        Return a string representation of ``v mod p!``, with ``p=self._prec()``

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._subset_label_text(4)
            '4 mod 6'
            sage: graph.view.prec += 3
            sage: graph._subset_label_text(5)
            '5 mod 120'
            sage: graph.view.prec += 3
            sage: graph._subset_label_text(6)
            '6 mod 8!'
        """
        if self._prec() <= 6:
            return "{} mod {}".format(v, factorial(self._prec()))
        return "{} mod {}!".format(v, self._prec())

    def _draw_subset_labels(self):
        """
        Draw the open subset labels (e.g. 0 mod 6, 2 mod 6, ...)

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: len(graph.canvas.find_all())
            274
            sage: graph._draw_subset_labels()
            sage: len(graph.canvas.find_all())
            298
            sage: graph.view.prec += 1
            sage: graph._draw_subset_labels()
            sage: len(graph.canvas.find_all())
            310
            sage: graph.view.prec += 1
            sage: graph._draw_subset_labels()
            sage: len(graph.canvas.find_all())
            326
        """
        n_blocks = prod(range(self.view.prec+1, self._prec()+1))
        y_up = self.border/ZZ(3)
        y_down = self.border + self.height + 2*self.border/ZZ(3)
        modulus = factorial(self._prec())
        view_modulus = factorial(self.view.prec)
        font = self._label_font()
        for i in range(n_blocks):
            a = Zhat(self.view.x + i*view_modulus, modulus)
            s, t = self._pixel_coordinates(a, is_x=True)
            x = self.border + (s+t)//2
            text = self._subset_label_text(a.value())
            self.canvas.create_text(x, y_up, text=text, font=font)
            self.canvas.create_text(x, y_down, text=text, font=font)
        x_left = self.border/ZZ(3)
        x_right = self.border + self.width + 2*self.border/ZZ(3)
        for i in range(n_blocks):
            a = Zhat(self.view.y + i*view_modulus, modulus)
            s, t = self._pixel_coordinates(a, is_x=False)
            y = self.border + self.height - (s+t)//2
            text = self._subset_label_text(a.value())
            self.canvas.create_text(x_left, y, text=text, font=font, angle=90)
            self.canvas.create_text(x_right, y, text=text, font=font, angle=-90)

    def _draw_footer(self):
        """
        Draw the footer with "currently viewing: (.. mod ..!) x (.. mod ..!)"

        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: len(graph.canvas.find_all())
            274
            sage: graph._draw_footer()
            sage: len(graph.canvas.find_all())
            275
        """
        text = "currently viewing: {}".format(self.view)
        x_middle = self.border+self.width/2
        y_bottom = 2*self.border + self.height + self.footer/2
        self.canvas.create_text(x_middle, y_bottom, text=text)

    def _draw(self):
        """
        Delete everything on the current canvas and redraw the whole graph
        
        TESTS::

            sage: graph = ProfiniteGraph(ProfiniteFibonacci())
            sage: graph._draw_footer()
            sage: len(graph.canvas.find_all())
            275
            sage: graph._draw()
            sage: len(graph.canvas.find_all())
            274
        """
        self.canvas.delete("all")
        self._draw_approximations()
        self._draw_identity_line()
        self._draw_axis()
        self._draw_axis_labels()
        self._draw_subset_labels()
        self._draw_footer()
        self.highlight_rect = self.canvas.create_rectangle(
            0, 0, 1, 1, width=0, fill=self.highlight_color, stipple="gray12"
        )


