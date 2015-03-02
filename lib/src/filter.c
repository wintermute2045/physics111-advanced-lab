/*******************************************************************************
 * This code is for filterting lists in C.
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
#include <gsl/gsl_math.h>

/*******************************************************************************
 * finite_pairs
 * Determines whether or not two numbers are both finite.
 ******************************************************************************/
int finite_pairs(double num1, double num2, void *null){
    int finite = 0;

    if(gsl_finite(num1) && gsl_finite(num2)) finite = 1;

    return finite;
}

/*******************************************************************************
 * filter_d
 * Copies items of src to dest provide they satisfy a certain condition.
 * Input:
 *      dest:           Address of the destination array.
 *      src:            Address of the source array.
 *      size:           Size of the array.
 *      func:           Function that gets applied.
 *      params:         Parameters of the function that gets applied.
 * Output:
 *      npts:           Size of the new list.
 ******************************************************************************/
size_t filter_d(double *dest, double *src, size_t size,
        int (*func)(double, void *), void *params){
    size_t i, npts;

    for (i=npts=0; i<size; i++){
        if ((*func)(src[i], params)){
            dest[npts] = src[i];
            npts++;
        }
    }
    return npts;
}

/*******************************************************************************
 * filter_pairs_d
 * Input:
 *      dest1:          Address of the first destination array.
 *      dest2:          Address of the second destination array.
 *      src1:           Address of the first source array.
 *      src2:           Address of the second source array.
 *      size:           Size of the array.
 *      func:           Function that gets applied.
 *      params:         Parameters of the function that gets applied.
 ******************************************************************************/
size_t filter_pairs_d(double *dest1, double *dest2, double *src1, double *src2,
        size_t size, int (*func)(double, double, void *), void *params){
    size_t i, npts;

    for (i=npts=0; i<size; i++){
        if ((*func)(src1[i], src2[i], params)){
            dest1[npts] = src1[i];
            dest2[npts] = src2[i];
            npts++;
        }
    }

    return npts;
}
