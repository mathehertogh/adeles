"""
Load the SageMath package "Computing with adèles and idèles" at the GitHub
repository https://github.com/mathehertogh/adeles via the web.
This SageMath package is based on and part of [Her2021].

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

AUTHOR: Mathé Hertogh (2021-07)
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

files = [
	"profinite_integer.py",
	"profinite_number.py",
	"completion.py",
	"adele.py",
	"multiplicative_padic.py",
	"idele.py",
	"ray_class_group.py",
	"profinite_function.py",
	"profinite_graph.py",
	"matrix.py",
	"shimura.py",
	"modular.py"
]

repository = "https://raw.githubusercontent.com/mathehertogh/adeles/main/"

load([repository + file for file in files])
