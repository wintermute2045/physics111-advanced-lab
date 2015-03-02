/*******************************************************************************
 * This code is for mapping lists in C.
 * Copyright (C) 2013  Rachel Domagalski: idomagalski@berkeley.edu
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*------------------------------------------------------------------------------
 * Redefining math functions to work with mapping.
 *----------------------------------------------------------------------------*/
double map_log(double x, void *null){
    return log(x);
}

double map_log10(double x, void *null){
    return log10(x);
}

double map_pow(double x, void *p){
    double *y = (double *) p;
    return pow(x, *y);
}

double map_sqrt(double x, void *null){
    return sqrt(x);
}

double map_pow2(double x, double y, void *null){
    return pow(x, y);
}

/*------------------------------------------------------------------------------
 * End of redefining math functions to work with mapping.
 *----------------------------------------------------------------------------*/

/*******************************************************************************
 * map_d
 * Applies the function f to all elements in an array of doubles.
 * Input:
 *      dest:           Address of the destination array.
 *      src:            Address of the source array.
 *      size:           Size of the array.
 *      func:           Function that gets applied.
 *      params:         Parameters of the function that gets applied.
 ******************************************************************************/
void map_d(double *dest, double *src, long size,
        double (*func)(double, void *), void *params){
    long i;

    for (i=0; i<size; i++){
        dest[i] = (*func)(src[i], params);
    }
}

/*******************************************************************************
 * map2_d
 * Applies the function f to all elements in an array of doubles.
 * Input:
 *      dest:           Address of the destination array.
 *      src1:           Address of the first source array.
 *      src2:           Address of the second source array.
 *      size:           Size of the array.
 *      func:           Function that gets applied.
 *      params:         Parameters of the function that gets applied.
 ******************************************************************************/
void map2_d(double *dest, double *src1, double *src2, long size,
        double (*func)(double, double, void *), void *params){
    long i;

    for (i=0; i<size; i++){
        dest[i] = (*func)(src1[i], src2[i], params);
    }
}
