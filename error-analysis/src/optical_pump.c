/*******************************************************************************
 * This code is for problem 5 of the EAX homework for Physics 111 Advanced Lab
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
#include <gsl/gsl_fit.h>

/*******************************************************************************
 * least_sq
 * Calculates the least square fit for a series of data.
 * Input:
 *      xdata:          Data along the x axis
 *      ydata:          Data along the y axis
 *      weights:        Weights of the data, NULL for unweighted fit
 *      npts:           Number of data points
 * Output:
 *      intercept:      Fitted y intercept of the data
 *      err_b:          Error on the intercept
 *      slope:          Fitted slope of the data
 *      err_m:          Error on the slope
 *      chisq:          Chi^2 of the fit.
 ******************************************************************************/
void least_sq(double *xdata,
              double *ydata,
              double *weights,
              long npts,
              double *intercept,
              double *err_b,
              double *slope,
              double *err_m,
              double *chisq)
{
    double cov00 = 0, cov01 = 0, cov11 = 0, weighted_mean, sum_weights;
    long i;

    if (weights == NULL){
        gsl_fit_linear(xdata, 1, ydata, 1, npts,
                intercept, slope, &cov00, &cov01, &cov11, chisq);
        *err_m = *err_b = -1; /* Errors for unweighted fit are irrelevant. */
    } else {
        gsl_fit_wlinear(xdata, 1, weights, 1, ydata, 1, npts,
                intercept, slope, &cov00, &cov01, &cov11, chisq);

        /* compute weighted mean of xdata */
        weighted_mean = 0;
        sum_weights = 0;
        for (i=0; i<npts; i++) {
            weighted_mean += weights[i] * xdata[i];
            sum_weights += weights[i];
        }
        weighted_mean /= sum_weights;

        /* get errors on coefficients */
        *err_b = *err_m = 0;
        for (i=0; i<npts; i++){
            *err_b += weights[i];
            *err_m += weights[i] * pow(xdata[i] - weighted_mean, 2);
        }
        *err_b = 1.0 / sqrt(*err_b);
        *err_m = 1.0 / sqrt(*err_m);
    }
}

/*******************************************************************************
 * fits_least_sq
 * This function computes the least square fit for data and then writes the
 * results/input data to a fits file.
 * Input:
 *      infile:         Input file for column info
 *      outfile:        Output file to be written to
 *      xdata:          The x-axis data
 *      ydata:          The y-axis data
 *      errors:         Errors on the ydata, NULL for unweighted fit.
 *      npts:           Number of data points.
 ******************************************************************************/
void fits_least_sq(fitsfile *infile, fitsfile *outfile,
                   double *xdata,
                   double *ydata,
                   double *errors,
                   long npts)
{
    double intercept = 0, err_b = 0, slope = 0, err_m = 0, chisq = 0;
    double chisq_prob;
    double *weights;
    int status = 0;
    long i;

    /* Create data table in outfile */
    fits_create_tbl(outfile, BINARY_TBL, 0, 0, NULL, NULL, NULL, NULL, &status);
    fits_copy_col(infile, outfile, 1, 1, 1, &status);
    fits_copy_col(infile, outfile, 2, 2, 1, &status);
    fits_write_col(outfile, TDOUBLE, 1, 1L, 1L, npts, xdata, &status);
    fits_write_col(outfile, TDOUBLE, 2, 1L, 1L, npts, ydata, &status);
    if (status){
        printf("Error copying data from infile to outfile.\n");
        exit(status);
    }

    /* calculate weights */
    if (errors == NULL)
        weights = NULL;
    else {
        fits_insert_col(outfile, 3, "sigma_f", "D", &status);
        fits_write_col(outfile, TDOUBLE, 3, 1L, 1L, npts, errors, &status);
        if (status){
            printf("Error writing sigma to fits table.\n");
            exit(status);
        }
        weights = (double *) malloc(npts * sizeof(double));
        for (i=0; i<npts; i++) weights[i] = 1.0 / pow(errors[i], 2);
    }

    /* compute least squares fit */
    least_sq(xdata, ydata, weights, npts,
            &intercept, &err_b, &slope, &err_m, &chisq);
    chisq_prob = gsl_cdf_chisq_P(chisq, npts-2);

    /* Display some results */
    printf("Data fit to the curve: y = %g * x + %g\n", slope, intercept);
    if (errors != NULL){
        printf("Error on the slope: %g\n", err_m);
        printf("Error on the intercept: %g\n", err_b);
    }
    printf("Chi^2 / ndf = %g / %ld = %g\n", chisq, npts - 2, chisq / (npts-2));
    printf("Probability for lower chi^2 %g for %ld degrees of freedom: %g\n",
            chisq, npts - 2, chisq_prob);

    /* Write results to fits table */
    fits_write_key(outfile, TDOUBLE, "M", &slope, "slope", &status);
    fits_write_key(outfile, TDOUBLE, "B", &intercept, "intercept", &status);
    if (errors != NULL){
        fits_write_key(outfile, TDOUBLE, "ERRM", &err_m, "error of slope",
                &status);
        fits_write_key(outfile, TDOUBLE, "ERRM", &err_b, "error of intercept",
                &status);
    }
    fits_write_key(outfile, TDOUBLE, "CHISQ", &chisq, "chi squared", &status);
    if (status){
        printf("Unable to write to fits header.\n");
        exit(status);
    }

    free(weights);
}

/*******************************************************************************
 * optimized_sigma
 * Given the slope and intercept of a linear fit, compute the optimal sigma and
 * use it to calculate the error on the slope and intercept.
 * Input:
 *      fptr:           Fits file containing fit information
 *      xdata:          x axis data
 *      ydata:          y axis data
 *      npts:           number of points
 ******************************************************************************/
void optimized_sigma(fitsfile *fptr, double *xdata, double *ydata, long npts){
    char *comment;
    double slope = 0, intercept = 0;
    double sigma, err_m, err_b;
    double mean;
    int status = 0;
    long i;

    /* Get slope and intercept from fits file */
    comment = (char *) malloc(81 * sizeof(char));
    fits_read_key(fptr, TDOUBLE, "M", &slope, comment, &status);
    fits_read_key(fptr, TDOUBLE, "B", &intercept, comment, &status);
    free(comment);
    if (status){
        printf("Unable to read fits header.\n");
        exit(status);
    }

    /* Compute optimal sigma */
    sigma = 0;
    for (i=0; i<npts; i++){
        sigma += pow(ydata[i] - slope * xdata[i] - intercept, 2);
    }
    sigma = sqrt(sigma / (npts - 2));

    /* compute mean of xdata */
    mean = 0;
    for (i=0; i<npts; i++) mean += xdata[i];
    mean /= npts;

    /* Compute errors on fitted coefficients */
    err_m = err_b = 0;
    for (i=0; i<npts; i++){
        err_b += pow(1.0 / sigma, 2);
        err_m += pow((xdata[i] - mean) / sigma, 2);
    }
    err_b = 1.0 / sqrt(err_b);
    err_m = 1.0 / sqrt(err_m);

    /* Display the results */
    printf("Computing error on data using chi^2 / ndf == 1 on the fit:\n");
    printf("y = %g * x + %g\n", slope, intercept);
    printf("Error on the y data: %g\n", sigma);
    printf("Error on the slope: %g\n", err_m);
    printf("Error in the intercept: %g\n", err_b);

    /* Write to the fits file */
    fits_write_key(fptr, TDOUBLE, "OPTERR", &sigma, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "OPTERRM", &err_m, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "OPTERRB", &err_b, NULL, &status);
    if (status){
        printf("Unable to write to fits header.\n");
        exit(status);
    }
}

/*******************************************************************************
 * optical_pump
 * This file reads data from some optical pump experiment and computes several
 * linear fits with certain errors. 
 ******************************************************************************/
int main(int argc, char **argv){
    double *current, *frequency, *errors, sigma;
    int status = 0, anynul = 0;
    long npts = 0;
    long i;
    /* Variables for reading in options */
    int opt, hdutype;
    fitsfile *infile, *outfile;
    char *infilename, *outfilename;

    /* Parse options */
    infilename = NULL;
    outfilename = NULL;
    while ((opt = getopt(argc, argv, "i:o:")) != -1){
        switch (opt){
            case 'i':
                infilename = (char *) malloc((strlen(optarg)+1)*sizeof(char));
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
    printf("Reading input from %s.\n", infilename);
    fits_open_table(&infile, infilename, READONLY, &status);
    if (status){
        printf("Error opening file %s.\n", infilename);
        exit(status);
    }
    fits_get_num_rows(infile, &npts, &status);
    if (status){
        printf("Error reading fits table.\n");
        exit(status);
    }
    current = (double *) malloc(npts * sizeof(double));
    frequency = (double *) malloc(npts * sizeof(double));
    errors = (double *) malloc(npts * sizeof(double)); /* assign later */
    fits_read_col(infile, TDOUBLE, 1, 1L, 1L, npts, NULL,
            current, &anynul, &status);
    fits_read_col(infile, TDOUBLE, 2, 1L, 1L, npts, NULL,
            frequency, &anynul, &status);
    if (status){
        printf("Error reading data from fits table.\n");
        exit(status);
    }

    /* Create output file for stroring results */
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

     /* Least square fit for unweighted data. */
    printf("\nComputing unwieighted least squares fit.\n");
    fits_least_sq(infile, outfile, current, frequency, NULL, npts);

    /* Least square fit for data with sigma = 0.01 MHz */
    sigma = 0.01;
    printf("\nComputing least squares fit with sigma_f = %g.\n", sigma);
    for (i=0; i<npts; i++) errors[i] = sigma;
    fits_least_sq(infile, outfile, current, frequency, errors, npts);

    /* Least square fit for data with sigma = 1 MHz */
    sigma = 1.0;
    printf("\nComputing least squares fit with sigma_f = %g.\n", sigma);
    for (i=0; i<npts; i++) errors[i] = sigma;
    fits_least_sq(infile, outfile, current, frequency, errors, npts);
    
    /* Least square fit for data with sigma = 0.3 + 0.3*f */
    sigma = 0.03;
    printf("\nComputing least squares fit with sigma_f = %g + %g * frequency.\n",
            sigma, sigma);
    for (i=0; i<npts; i++) errors[i] = sigma + sigma * frequency[i];
    fits_least_sq(infile, outfile, current, frequency, errors, npts);
    printf("\n");

    /* compute errors assuming first fit is optimal */
    fits_movabs_hdu(outfile, 2, &hdutype, &status);
    optimized_sigma(outfile, current, frequency, npts);

    /* Close files and exit */
    printf("\nClosing file %s.\n", infilename);
    fits_close_file(infile, &status);
    if (status){
        printf("Error closing file %s.\n", infilename);
        exit(status);
    }
    printf("Closing file %s.\n\n", outfilename + 1);
    fits_close_file(outfile, &status);
    if (status){
        printf("Error closing file %s.\n", outfilename + 1);
        exit(status);
    }
    free(infilename);
    free(outfilename);
    free(current);
    free(frequency);
    free(errors);
    return 0;
}
