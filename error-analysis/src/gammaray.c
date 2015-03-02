/*******************************************************************************
 * This code is for problem 4 of the EAX homework for Physics 111 Advanced Lab
 * at UC Berkeley.
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
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <fitsio.h>
#include <gsl/gsl_cdf.h>
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
 *      mean:           mean of the gaussian
 *      sigma:          standard deviation of the gaussian
 * Output:
 *      result:         Value of the gaussian at a given point.
 ******************************************************************************/
double gaussian(double arg, double magnitude, double mean, double sigma){
    double result;
    result = -pow((arg - mean)/sigma, 2) / 2.0;
    result = magnitude * exp(result);
    return result;
}

/*******************************************************************************
 * double *gamma_stats
 * Get the mean, standard deviations, and uncertainties of those quantities.
 * Input:
 *      fptr:           Fits file containing data.
 ******************************************************************************/
void gamma_stats(fitsfile *fptr){
    double mean, sigma, err_mean, err_sigma;
    double data;
    int anynul, status = 0;
    long i, num_samples;

    fits_get_num_rows(fptr, &num_samples, &status);
    if (status){
        printf("Error retrieving info from fits table.\n");
        exit(status);
    }

    /* Compute mean */
    mean = 0;
    for (i=0; i<num_samples; i++){
        fits_read_col(fptr, TDOUBLE, 1, i+1, 1L, 1L, NULL,
                &data, &anynul, &status);
        if (status){
            printf("Error reading fits table.\n");
            exit(status);
        }
        mean += data;
    }
    mean /= num_samples;

    /* Compute standard deviation */
    sigma = 0;
    for (i=0; i<num_samples; i++){
        fits_read_col(fptr, TDOUBLE, 1, i+1, 1L, 1L, NULL,
                &data, &anynul, &status);
        if (status){
            printf("Error reading fits table.\n");
            exit(status);
        }
        sigma += pow(mean - data, 2);
    }
    sigma /= num_samples - 1;
    sigma = sqrt(sigma);

    err_mean = sigma / sqrt(num_samples);
    err_sigma = sigma / sqrt(2 * num_samples);

    /* Print the results */
    printf("Mean: %g\n", mean);
    printf("Error on the mean: %g\n", err_mean);
    printf("Standard deviation: %g\n", sigma);
    printf("Error on the standard deviation: %g\n\n", err_sigma);

    fits_write_key(fptr, TDOUBLE, "MEAN", &mean, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "SIGMA", &sigma, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "ERRMEAN", &err_mean, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "ERRSIGMA", &err_sigma, NULL, &status);
    if (status){
        printf("Error writing to fits header.\n");
        exit(status);
    }
}

/*******************************************************************************
 * gsl_histogram *fill_hist
 * Fit a histogram to a gaussian and write the fitted mean and sigma to the fits
 * header. Based on the gsl documentation, the range member of the gsl_histogram
 * struct seems to be consistent with the bins return value of the pyplot hist
 * function.
 * Intput:
 *      fptr:           Data used to fill the histogram
 *      nbins:          Number of bins in the histogram
 *      min:            Lower bound on histogram
 *      max:            Upper bound on histogram
 * Output:
 *      hist:           Histogram of the recorded data.
 ******************************************************************************/
gsl_histogram *fill_hist(fitsfile *fptr, long nbins, double min, double max){
    double data;
    int status = 0, anynul;
    long i, nrows;
    gsl_histogram *hist;

    fits_get_num_rows(fptr, &nrows, &status);
    if (status){
        printf("Error retrieving info from fits table.\n");
        exit(status);
    }

    /* Allocate and set ranges */
    hist = gsl_histogram_alloc(nbins);
    if (hist == NULL){
        printf("Failed to allocate histogram.\n");
        exit(1);
    }
    gsl_histogram_set_ranges_uniform(hist, min, max);

    /* Fill histogram */
    for (i=0; i<nrows; i++){
        fits_read_col(fptr, TDOUBLE, 1, i+1, 1L, 1L, NULL,
                &data, &anynul, &status);
        if (status){
            printf("Error reading fits table.\n");
            exit(status);
        }
        gsl_histogram_increment(hist, data);
    }

    return hist;
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
int gaus_f(const gsl_vector *x, void *data, gsl_vector *f){
    long i, nbins, nonzero;
    double *xdata, *ydata, xpoint, bin_width, err;
    double magnitude, mean, sigma, gaus;
    gsl_histogram *hist;

    /* Get previous parameters */
    magnitude = gsl_vector_get(x, 0);
    mean      = gsl_vector_get(x, 1);
    sigma     = gsl_vector_get(x, 2);

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
            gaus = gaussian(xpoint, magnitude, mean, sigma);
            err = sqrt(ydata[i]);
            gsl_vector_set(f, nonzero, (gaus - ydata[i]) / err);
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
int gaus_df(const gsl_vector *x, void *data, gsl_matrix *J){
    long i, nbins, nonzero;
    double *xdata, *ydata, xpoint, bin_width, err, jacobian;
    double magnitude, mean, sigma, gaus;
    gsl_histogram *hist;

    /* Get previous parameters */
    magnitude = gsl_vector_get(x, 0);
    mean      = gsl_vector_get(x, 1);
    sigma     = gsl_vector_get(x, 2);

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
            gaus = gaussian(xpoint, 1, mean, sigma);
            gsl_matrix_set(J, nonzero, 0, gaus / err);
            gaus = gaussian(xpoint, magnitude, mean, sigma);
            jacobian = (xpoint - mean) * gaus / pow(sigma, 2);
            gsl_matrix_set(J, nonzero, 1, jacobian / err);
            jacobian = pow(xpoint - mean, 2) * gaus / pow(sigma, 3);
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
int gaus_fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J){
    gaus_f(x, data, f);
    gaus_df(x, data, J);
    return GSL_SUCCESS;
}

/*******************************************************************************
 * fit_gaussian
 * Fit data to a guassian and write the results to a fits file. Ideally, this
 * should give the same results as scipy.optimize.curve_fit.
 * Input:
 *      fptr:           Fits file to be written to
 *      xdata:          X axis of the data
 *      ydata:          Y axis of the data
 *      npts:           Number of data points
 * Output:
 *      fit_params:     Fit parameters
 *      chisq_ndf:      chi^2 per degree of freedom
 ******************************************************************************/
void fit_gaussian(fitsfile *fptr, gsl_histogram *hist){
    double tol;
    double *hbin, *hrange, bin_width, xdata, min, max;
    double magnitude, mean, sigma;
    double error, chisq, chisq_ndf, ythr, chisq_prob;
    int status;
    long gpars, nonzero, nbins, ndf;
    long i;
    gsl_vector *pars;
    gsl_multifit_fdfsolver *gfit;
    gsl_multifit_function_fdf gaus;
    const gsl_multifit_fdfsolver_type *ftype;

    printf("Fitting the histogram using nonlinear least squares.\n");

    /* Allowed relative error is what scipy uses */
    tol = 1.49012e-8;

    /* get number of bins containing data */
    nbins = hist -> n;
    hbin = hist -> bin;
    hrange = hist -> range;
    nonzero = 0;
    for (i=0; i<nbins; i++){
        if (hbin[i]) nonzero++;
    }

    /* Set the function */
    gaus.f = &gaus_f;
    gaus.df = &gaus_df;
    gaus.fdf = &gaus_fdf;
    gaus.n = nonzero;
    gaus.p = 3;
    gaus.params = hist;

    /* Initialize the solver */
    gpars = 3;
    pars = gsl_vector_alloc(gpars);
    gsl_vector_set_all(pars, 1.0);
    ftype = gsl_multifit_fdfsolver_lmsder;
    gfit = gsl_multifit_fdfsolver_alloc(ftype, nonzero, gpars);
    gsl_multifit_fdfsolver_set(gfit, &gaus, pars);

    /* loop the solver and solve this thing */
    do {
        status = gsl_multifit_fdfsolver_iterate(gfit);
        status = gsl_multifit_test_delta(gfit -> dx, gfit -> x, 0, tol);
    } while (status == GSL_CONTINUE);

    magnitude = gsl_vector_get(gfit -> x, 0);
    mean = gsl_vector_get(gfit -> x, 1);
    /* The fitted sigma might be negative, but it is squared when computing the
     * gaussian, so taking the absolute value of sigma is ok */
    sigma = fabs(gsl_vector_get(gfit -> x, 2));

    /* Compute the chi^2 */
    min = hrange[0];
    max = hrange[nbins];
    bin_width = (max - min) / nbins;
    chisq = 0;
    for (i = 0; i<nbins; i++){
        if (hbin[i]){
            xdata = hrange[i] + bin_width/2.0;
            error = sqrt(hbin[i]);
            ythr = gaussian(xdata, magnitude, mean, sigma);
            chisq += pow((hbin[i] - ythr)/error, 2);
        }
    }

    ndf = nonzero - gpars;
    chisq_ndf = chisq / ndf;
    chisq_prob = gsl_cdf_chisq_P(chisq, ndf);

    /* Display some results */
    printf("Mean from the fit: %g\n", mean);
    printf("Sigma from the fit: %g\n", sigma);
    printf("Chi^2 / ndf = %g / %ld = %g\n", chisq, ndf, chisq_ndf);
    printf("Probability of getting a lower chi^2: %g\n\n", chisq_prob);

    /* Write results to a fits file */
    status = 0; /* reset from fit */
    fits_write_key(fptr, TDOUBLE, "HMAG", &magnitude, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "HMEAN", &mean, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "HSIGMA", &sigma, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "CHISQ", &chisq, NULL, &status);
    fits_write_key(fptr, TLONG, "NDF", &ndf, NULL, &status);
    if (status){
        printf("Error writing to fits header.\n");
        exit(status);
    }

    /* Free the solver's memory */
    gsl_vector_free(pars);
    gsl_multifit_fdfsolver_free(gfit);
}

/*******************************************************************************
 * gammaray
 * Reads from a data file supposedly containing energies from some gamma ray
 * experiment. The data is assumed to be on the order of magnitude of 1. This
 * function produces a histogram, and fits it to a gaussian. The fit parameters
 * are also compared to the mean and sigma that are computed directly.
 ******************************************************************************/
int main(int argc, char **argv){
    int status = 0;
    double energy, max_energy, min_energy;
    long num_samples = 0;
    /* Variables for reading in options */
    int opt;
    FILE *infile;
    fitsfile *outfile;
    char *infilename, *outfilename;
    /* Variables for fitsio */
    char *ttype[] = {"energy"};
    char *tform[] = {"D"};
    char *tunit[] = {NULL};
    /* Variables for fitting a histogram to a gaussian */
    long nbins;
    gsl_histogram *hist;

    /* Parse options */
    infilename = NULL;
    outfilename = NULL;
    while ((opt = getopt(argc, argv, "i:o:")) != -1){
        switch (opt){
            case 'i':
                infilename = (char *) malloc((strlen(optarg)+1)*sizeof(char));
                /* Prepend ! to filename to force overwriting it. */
                strcpy(infilename, optarg);
                break;
            case 'o':
                outfilename = (char *) malloc((strlen(optarg)+2)*sizeof(char));
                /* Prepend ! to filename to force overwriting it. */
                outfilename[0] = '!';
                strcat(outfilename, optarg);
                break;
            case '?':
                if (optopt == 'i'){
                    printf("Error: No input file specified.\n");
                    exit(1);
                } else if (optopt == 'o'){
                    printf("Error: No output file specified.\n");
                    exit(1);
                } else {
                    printf("Unknown option -%c. Exiting.\n", optopt);
                    exit(1);
                }
                break;
            default:
                exit(1);
        }
    }
    if (infilename == NULL){
        printf("Error: No input file specified.\n");
        exit(1);
    }
    if (outfilename == NULL){
        printf("Error: No output file specified.\n");
        exit(1);
    }

    /* Open input file for reading */
    infile = fopen(infilename, "r");
    if (infile == NULL){
        printf("Error opening file %s for reading.\n", infilename);
        exit(1);
    }

    /* Create file and generate results for it. */
    printf("Creating output file %s.\n", outfilename + 1);
    fits_create_file(&outfile, outfilename, &status);
    if (status){
        printf("Error creating file %s.\n", outfilename + 1);
        exit(status);
    }
    fits_open_file(&outfile, outfilename+1, READWRITE, &status);
    if (status){
        printf("Error opening file %s.\n", outfilename + 1);
        exit(status);
    }

    /* Create and fill fits table */
    fits_create_tbl(outfile, BINARY_TBL, 0, 1,
            ttype, tform, tunit, NULL, &status);
    if (status){
        printf("Error creating fits table in %s.\n", outfilename + 1);
        exit(status);
    }

    printf("Reading input from %s to %s.\n\n", infilename, outfilename + 1);
    num_samples = 0;
    while (fscanf(infile, "%lf", &energy) != EOF){
        /* Test for min an max */
        if (num_samples){
            if (energy > max_energy) max_energy = energy;
            if (energy < min_energy) min_energy = energy;
        } else {
            min_energy = max_energy = energy;
        }

        /* Write to fits file */
        num_samples++;
        fits_write_col(outfile, TDOUBLE, 1L, num_samples,
                1L, 1L, &energy, &status);
        if (status){
            printf("Error writing to fits table in %s.\n", outfilename + 1);
            exit(status);
        }
    }
    status = fclose(infile);
    if (status){
        printf("Error closing file %s.\n", infilename);
        exit(status);
    }

    /* Get mean and standard deviation */
    printf("Computing statistics.\n");
    gamma_stats(outfile);

    /* Fill histogram to fit gaussian to */
    min_energy = floor(min_energy * 1e1) * 1e-1;
    max_energy = (floor(max_energy * 1e1) + 1) * 1e-1;
    nbins = (long) (100 * (max_energy - min_energy));
    hist = fill_hist(outfile, nbins, min_energy, max_energy);
    fits_write_key(outfile, TLONG, "HBINS", &nbins, NULL, &status);
    fits_write_key(outfile, TDOUBLE, "HMIN", &min_energy, NULL, &status);
    fits_write_key(outfile, TDOUBLE, "HMAX", &max_energy, NULL, &status);
    if (status){
        printf("Error writing to fits header.\n");
        exit(status);
    }

    /* Fit gaussian */
    fit_gaussian(outfile, hist);

    /* Close file and exit */
    printf("Closing file %s.\n\n", outfilename + 1);
    fits_close_file(outfile, &status);
    if (status){
        printf("Error closing file %s.\n", outfilename + 1);
        exit(status);
    }
    free(infilename);
    free(outfilename);
    gsl_histogram_free(hist);
    return 0;
}
