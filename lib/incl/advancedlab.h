/*******************************************************************************
 * This code declares functions that might be commonly used throughout the labs
 * in the UC Berkeley Physics Advanced Lab course.
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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_histogram.h>

#ifndef _ADVLAB_H
#define _ADVLAB_H 1

/* In fitgaus.c */

double gaussian(double arg, double magnitude, double mean, double sigma);
int gaus_f(const gsl_vector *x, void *data, gsl_vector *f);
int gaus_df(const gsl_vector *x, void *data, gsl_matrix *J);
int gaus_fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J);
gsl_vector *fit_gaussian(gsl_histogram *hist,
        double *chisq, long *ndf, gsl_matrix *covar);

/* In fit_exp.c */

double expdecay(double arg, double magnitude, double mean, double sigma);
int exp_f(const gsl_vector *x, void *data, gsl_vector *f);
int exp_df(const gsl_vector *x, void *data, gsl_matrix *J);
int exp_fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J);
gsl_vector *fit_expdecay(gsl_histogram *hist,
        double *chisq, long *ndf, gsl_matrix *covar);

/* In map.c */

double map_log(double x, void *null);
double map_log10(double x, void *null);
double map_pow(double x, double y, void *null);
void map_d(double *dest, double *src, long size,
        double (*func)(double, void *), void *params);
void map2_d(double *dest, double *src1, double *src2, long size,
        double (*func)(double, double, void *), void *params);

/* In filter.c */

int finite_pairs(double num1, double num2, void *null);
size_t filter_d(double *dest, double *src, size_t size,
        int (*func)(double, void *), void *params);
size_t filter_pairs_d(double *dest1, double *dest2, double *src1, double *src2,
        size_t size, int (*func)(double, double, void *), void *params);

/* In misc.c */
void init_prog(char *progname);
void exit_prog(char *progname);

#endif
