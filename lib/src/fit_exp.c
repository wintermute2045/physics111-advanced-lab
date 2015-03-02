/*******************************************************************************
 * This code is for fitting gaussians using GSL.
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
#include <fitsio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_multifit_nlin.h>

/*******************************************************************************
 * gaussian
 * Computes the gaussian for given parameters
 * Input:
 *      arg:            Point to compute the gaussian at
 *      magnitude:      Height of the function
 *      lbda:           Constant in the exponent
 *      offset:         limit as t -> oo
 * Output:
 *      result:         Value of the gaussian at a given point.
 ******************************************************************************/
double expdecay(double arg, double magnitude, double lbda, double offset){
    double result;
    result = magnitude * exp(-lbda * arg);
    result += offset;
    return result;
}

/*******************************************************************************
 * gaus_f
 * Retrieves the components of the fit from x and uses them to compute values
 * to be stored in vector f. This function is required by the solver.
 * Input:
 *      x:              Vector with fit parameters
 *      data:           Data to fit to
 *      f:              set of computations of the data
 ******************************************************************************/
int exp_f(const gsl_vector *x, void *data, gsl_vector *f){
    long i, nbins, nonzero;
    double *xdata, *ydata, xpoint, bin_width, err;
    double magnitude, lbda, offset, ythr;
    gsl_histogram *hist;

    /* Get previous parameters */
    magnitude = gsl_vector_get(x, 0);
    lbda      = gsl_vector_get(x, 1);
    offset    = gsl_vector_get(x, 2);

    /* Histogram information */
    hist = (gsl_histogram *) data;
    nbins = hist -> n;
    xdata = hist -> range;
    ydata = hist -> bin;
    bin_width = (xdata[nbins] - xdata[0]) / nbins; /* for shifting to center */

    /* Set values of f */
    nonzero = 0;
    for (i=0; i<nbins; i++){
        if (ydata[i]){
            xpoint = xdata[i] + bin_width/2;
            ythr = expdecay(xpoint, magnitude, lbda, offset);
            err = sqrt(ydata[i]);
            gsl_vector_set(f, nonzero, (ythr - ydata[i]) / err);
            nonzero++;
        }
    }
    return GSL_SUCCESS;
}

/*******************************************************************************
 * gaus_df
 * Sets the jacobian of the gaussian
 * Input:
 *      x:              Vector with fit parameters
 *      data:           Data to fit to
 *      J:              Jacobian of the gaussian
 ******************************************************************************/
int exp_df(const gsl_vector *x, void *data, gsl_matrix *J){
    long i, nbins, nonzero;
    double *xdata, *ydata, xpoint, bin_width, err, jacobian;
    double magnitude, lbda;
    gsl_histogram *hist;

    /* Get previous parameters */
    magnitude = gsl_vector_get(x, 0);
    lbda      = gsl_vector_get(x, 1);

    /* Histogram information */
    hist = (gsl_histogram *) data;
    nbins = hist -> n;
    xdata = hist -> range;
    ydata = hist -> bin;
    bin_width = (xdata[nbins] - xdata[0]) / nbins; /* for shifting to center */

    /* Set values of J */
    nonzero = 0;
    for (i=0; i<nbins; i++){
        if (ydata[i]){
            xpoint = xdata[i] + bin_width/2;
            err = sqrt(ydata[i]);
            jacobian = exp(-lbda * xpoint);
            gsl_matrix_set(J, nonzero, 0, jacobian / err);
            jacobian *= -magnitude * xpoint;
            gsl_matrix_set(J, nonzero, 1, jacobian / err);
            jacobian = 1;
            gsl_matrix_set(J, nonzero, 2, jacobian / err);
            nonzero++;
        }
    }
    return GSL_SUCCESS;
}

/*******************************************************************************
 * gaus_fdf
 * Set up the system...
 *      x:              Vector with fit parameters
 *      data:           Data to fit to
 *      f:              set of computations of the data
 *      J:              Jacobian of the gaussian
 ******************************************************************************/
int exp_fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J){
    exp_f(x, data, f);
    exp_df(x, data, J);
    return GSL_SUCCESS;
}

/*******************************************************************************
 * fit_gaussian
 * Fit data to a guassian and return the results. Ideally, this should give the
 * same results as scipy.optimize.curve_fit.
 * Input:
 *      hist:           Histogram to fit the gaussian to
 * Output:
 *      chisq:          Chi^2 of the histogram
 *      ndf:            Number of degrees of freedom of the fit
 *      fit_params:     Fit parameters
 ******************************************************************************/
gsl_vector *fit_expdecay(gsl_histogram *hist,
        double *chisq, long *ndf, gsl_matrix *covar){
    double tol;
    double *hbin, *hrange, bin_width, xdata, min, max;
    double magnitude, lbda, offset;
    double error, ythr;
    int status;
    long epars, nonzero, nbins;
    long i;
    gsl_vector *pars, *fit_params;
    gsl_multifit_fdfsolver *efit;
    gsl_multifit_function_fdf exp_func;
    const gsl_multifit_fdfsolver_type *ftype;

    /* Allowed relative error is what scipy uses */
    /* tol = 1.49012e-8; scipy least squares default */
    tol = 1e-14;

    /* get number of bins containing data */
    nbins = hist -> n;
    hbin = hist -> bin;
    hrange = hist -> range;
    nonzero = 0;
    for (i=0; i<nbins; i++){
        if (hbin[i]) nonzero++;
    }

    /* Set the function */
    exp_func.f = &exp_f;
    exp_func.df = &exp_df;
    exp_func.fdf = &exp_fdf;
    exp_func.n = nonzero;
    exp_func.p = 3;
    exp_func.params = hist;

    /* Initialize the solver */
    epars = 3;
    pars = gsl_vector_alloc(epars);
    gsl_vector_set_all(pars, 1.0);
    ftype = gsl_multifit_fdfsolver_lmsder;
    efit = gsl_multifit_fdfsolver_alloc(ftype, nonzero, epars);
    gsl_multifit_fdfsolver_set(efit, &exp_func, pars);

    /* loop the solver and solve this thing */
    do {
        status = gsl_multifit_fdfsolver_iterate(efit);
        status = gsl_multifit_test_delta(efit -> dx, efit -> x, 0, tol);
    } while (status == GSL_CONTINUE);

    magnitude = gsl_vector_get(efit -> x, 0);
    lbda      = gsl_vector_get(efit -> x, 1);
    offset    = gsl_vector_get(efit -> x, 2);

    /* Compute the chi^2 */
    min = hrange[0];
    max = hrange[nbins];
    bin_width = (max - min) / nbins;
    *chisq = 0;
    for (i = 0; i<nbins; i++){
        if (hbin[i]){
            xdata = hrange[i] + bin_width/2.0;
            error = sqrt(hbin[i]);
            ythr = expdecay(xdata, magnitude, lbda, offset);
            *chisq += pow((hbin[i] - ythr)/error, 2);
        }
    }
    *ndf = nonzero - epars;

    /* Copy results to return vector */
    fit_params = gsl_vector_alloc(epars);
    gsl_vector_memcpy(fit_params, efit -> x);

    /* Compute the covariance matrix */
    gsl_multifit_covar(efit -> J, 0.0, covar);

    /* Free the solver's memory */
    gsl_vector_free(pars);
    gsl_multifit_fdfsolver_free(efit);

    /* Return the results of the fit */
    return fit_params;
}
