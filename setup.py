# -*- coding: utf-8 -*-
#
# Setup script for the ACOSA Python module.
# Copyright (C) 2016 Malte Ziebarth
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Imports:
from setuptools import setup
from setuptools.extension import Extension

from Cython.Distutils import build_ext
from Cython.Compiler.Options import _directive_defaults

_directive_defaults['linetrace'] = True
_directive_defaults['binding'] = True

import numpy as np




# Extension:
extensions=[]

extensions.append(Extension('acosa',
	sources=['acosa/ACOSA.pyx',
	         'acosa/basic_types.cpp',
	         'acosa/fortunes_sphere.cpp',
	         'acosa/vdtesselation.cpp',
	         'acosa/beach.cpp',
	         'acosa/order_parameter.cpp',
	         'acosa/spherics.cpp',
	         'acosa/convexhull.cpp',
	         'acosa/circleevent.cpp',
	         'acosa/geometricgraph.cpp',
	         'acosa/alphaspectrum.cpp'],
	include_dirs=[np.get_include(),'acosa'],
	extra_compile_args=['-std=c++14'],
	language='c++'))

extensions[0].cython_c_in_temp = False

# Setup:

setup(
	name='acosa',
	version='1.0.0',
	description="A Compilation of Spherical Algorithms",
	long_description="Algorithsm to calculate the Voronoi tesselation, \
	                  Delaunay triangulation, convex hull , and alpha shape of \
	                  a set of nodes on a sphere.",
	author='Malte J. Ziebarth',
	author_email='contact@fmvkb.de',
	packages=['acosa'],
	cmdclass = {'build_ext': build_ext},
	ext_modules=extensions,
	provides=['acosa'],
	scripts=[],
	license='GPLv3',
	)
