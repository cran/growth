/*
 *  growth : A Library of Special Functions for Repeated Measurements
 *  Copyright (C) 1998 J.K. Lindsey
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  SYNOPSIS
 *
 *    void flgamma_(double *x,double *y)
 *
 *  DESCRIPTION
 *
 *    This function allows a Fortran program to call the lgamma
 * function written in C.
 *
 */

#include <math.h>
#include "Fortran.h"
extern double lgamma(double x);
void F77_SYMBOL(flgamma)(double *x,double *y){
  *y=lgamma(*x);}
