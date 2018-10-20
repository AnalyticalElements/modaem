# Makefile v0.0.1 4/18/2001 VAK
#
# This is the makefile for ModAEM under Linux, using the Lahey F95 compiler
# NOTE: ModAEM makes use of the libraries BLAS and LAPACK, available at 
# Netlib (http://www.netlib.org)
#
# ModAEM 1.4pre1
# Copyright (c) 2001 WHPA Inc.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# Contact the author by e-mail at: vic@wittmanhydro.com
# Or by regular mail at:
# WHPA, Inc
# 320 W 8th St
# Bloomington, IN 47401

F95 = lf95
F95FLAGS = --ml cdecl --f95 -g --in --pca --tpp --wo
CC = gcc
CCFLAGS = -g
LD = lf95
LDFLAGS = --staticlink
AR = lf95
ARFLAGS = --nshared

INFLUENCE = i_well.o i_pond.o i_dipole.o i_linesink.o 
FUNCTIONS = f_well.o f_pond.o f_dipole.o
ELEMENTS = m_wl0.o m_pd0.o m_ls0.o m_ls1.o m_ls2.o m_hb0.o m_aqu.o m_in0.o m_aem.o m_inq.o m_wl1.o
UTILITY = u_constants.o u_io.o u_math.o u_matrix.o u_grid.o
ANALYSIS = a_grid.o a_tr0.o
LIBS = -lblas -llapack

all:	modaem.bin

modaem.bin : modaem.o ${INFLUENCE} ${FUNCTIONS} ${ELEMENTS} ${UTILITY} ${ANALYSIS}
	${LD} ${LDFLAGS} -o modaem.bin modaem.o ${INFLUENCE} ${FUNCTIONS} ${ELEMENTS} ${UTILITY} ${ANALYSIS} ${LIBS}

contour : ${CONTOUR}
	${LD} ${LDFLAGS} -o contour ${CONTOUR}

clean :
	rm -f *.o *.mod modaem.bin

distclean : clean
	rm -f *.a modtest

# Top-level ModAEM library
m_aem.o : m_aem.f90 u_constants.o u_io.o u_matrix.o f_well.o f_dipole.o f_pond.o \
          m_aqu.o m_wl0.o m_pd0.o m_ls0.o m_ls1.o m_ls2.o m_hb0.o m_wl1.o
	${F95} ${F95FLAGS} -c m_aem.f90
m_inq.o : m_inq.f90 u_constants.o u_io.o m_aem.o
	${F95} ${F95FLAGS} -c m_inq.f90

# Main programs
modaem.o : modaem.f90 u_constants.o u_io.o m_aqu.o m_wl0.o m_pd0.o m_ls0.o \
           m_ls1.o m_hb0.o m_wl1.o m_aem.o m_inq.o a_grid.o \
           a_tr0.o
	${F95} ${F95FLAGS} -c modaem.f90
tr0bln.o : tr0bln.f90 u_constants.o u_io.o
	${F95} ${F95FLAGS} -c tr0bln.f90

# Utility modules
u_constants.o : u_constants.f90
	${F95} ${F95FLAGS} -c u_constants.f90
u_io.o : u_io.f90 u_constants.o
	${F95} ${F95FLAGS} -c u_io.f90
u_math.o : u_math.f90 u_constants.o
	${F95} ${F95FLAGS} -c u_math.f90
u_matrix.o : u_constants.o u_io.o u_matrix.f90
	${F95} ${F95FLAGS} -c u_matrix.f90
u_grid.o : u_constants.o u_io.o u_grid.f90
	${F95} ${F95FLAGS} -c u_grid.f90
isvalidfp.o : isvalidfp.c
	${CC} ${CCFLAGS} -c isvalidfp.c

utility : ${UTILITY}

# Influence function modules
i_well.o : i_well.f90 u_constants.o
	${F95} ${F95FLAGS} -c i_well.f90
i_pond.o : i_pond.f90 u_constants.o
	${F95} ${F95FLAGS} -c i_pond.f90
i_linesink.o : i_linesink.f90 u_constants.o
	${F95} ${F95FLAGS} -c i_linesink.f90
i_dipole.o : i_dipole.f90 u_constants.o
	${F95} ${F95FLAGS} -c i_dipole.f90
i_rectarel.o : i_rectarel.f90 u_constants.o i_well.o i_dipole.o
	${F95} ${F95FLAGS} -c i_rectarel.f90
i_twl.o : i_twl.f90 u_constants.o u_math.o
	${F95} ${F95FLAGS} -c i_twl.f90

influence : ${INFLUENCE}

# Function collection modules
f_well.o : f_well.f90 u_constants.o u_io.o i_well.o
	${F95} ${F95FLAGS} -c f_well.f90
f_pond.o : f_pond.f90 u_constants.o u_io.o i_pond.o
	${F95} ${F95FLAGS} -c f_pond.f90
f_dipole.o : f_dipole.f90 u_constants.o u_io.o i_dipole.o i_linesink.o
	${F95} ${F95FLAGS} -c f_dipole.f90
functions : ${FUNCTIONS}

# Element modules
m_wl0.o : m_wl0.f90 f_well.o
	${F95} ${F95FLAGS} -c m_wl0.f90
m_wl1.o : m_wl1.f90 m_aqu.o f_well.o
	${F95} ${F95FLAGS} -c m_wl1.f90
m_pd0.o : m_pd0.f90 f_pond.o
	${F95} ${F95FLAGS} -c m_pd0.f90
m_ls0.o : m_ls0.f90 f_well.o f_dipole.o u_matrix.o
	${F95} ${F95FLAGS} -c m_ls0.f90
m_ls1.o : m_ls1.f90 f_well.o f_dipole.o m_aqu.o u_matrix.o 
	${F95} ${F95FLAGS} -c m_ls1.f90
m_ls2.o : m_ls2.f90 f_well.o f_dipole.o m_aqu.o u_matrix.o 
	${F95} ${F95FLAGS} -c m_ls2.f90
m_hb0.o : m_hb0.f90 f_dipole.o u_matrix.o
	${F95} ${F95FLAGS} -c m_hb0.f90
m_in0.o: m_in0.f90 f_dipole.o u_matrix.o
	${F95} ${F95FLAGS} -c m_in0.f90
m_aqu.o : m_aqu.f90 f_well.o f_dipole.o m_in0.o u_matrix.o
	${F95} ${F95FLAGS} -c m_aqu.f90
elements : ${ELEMENTS}

# Analysis modules
a_grid.o : a_grid.f90 m_aem.o
	${F95} ${F95FLAGS} -c a_grid.f90
a_tr0.o : a_tr0.f90 m_aem.o
	${F95} ${F95FLAGS} -c a_tr0.f90

