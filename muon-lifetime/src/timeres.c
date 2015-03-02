/*******************************************************************************
 * This code is part of the analysis for the muon lifetime lab.
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
#include <fitsio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_histogram.h>
#include <advancedlab.h>
#include <muonlifetime.h>

/*******************************************************************************
 * timeres
 * Compute the time resultion of the digitizer in the lab computers. The
 * resolution is the standard deviation of the gaussian of time measurements for
 * passing muons.
 ******************************************************************************/
int main(int argc, char **argv){
    init_prog("timeres");
    char *infilename;
    double chisq = 0;
    double hmax; /* Max time delay allowed in histogram */
    double magnitude, mean, sigma;
    int opt, status = 0, hdutype = 0;
    long ndf = 0, nbins;
    cutmaker *cut;
    fitsfile *infile;
    gsl_matrix *covar;
    gsl_vector *fitpars;
    gsl_histogram *hist;

    /* Parse options */
    infilename = NULL;
    while ((opt = getopt(argc, argv, "i:")) != -1){
        switch (opt){
            case 'i':
                infilename = (char *) malloc((strlen(optarg)+1)*sizeof(char));
                strcpy(infilename, optarg);
                break;
            case '?':
                if (optopt == 'i'){
                    printf("Error: No input file specified.\n");
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

    /* Open the input file */
    printf("Calculating the time resolution from %s.\n\n", infilename);
    fits_open_file(&infile, infilename, READWRITE, &status);
    fits_movabs_hdu(infile, 2, &hdutype, &status);
    if (status || hdutype != BINARY_TBL){
        printf("Error opening file %s.\n", infilename);
        exit(status);
    }

    /* Define generic cutmaker struct */
    cut              = (cutmaker *) malloc(sizeof(cutmaker));
    cut -> vmin      = -1;
    cut -> vmax      = -1;
    cut -> pulse     = IGNORE;
    cut -> fix_time  = NULL;
    cut -> params    = NULL;
    cut -> apply_fix = 0;

    /* Fill histogram and fit gaussian */
    hmax = 1e-7;
    hmax *= 1e6;
    nbins = 25;
    covar = gsl_matrix_alloc(3, 3);
    hist = fill_hist(infile, nbins, hmax, 1, 7, cut);
    fitpars = fit_gaussian(hist, &chisq, &ndf, covar);
    magnitude = gsl_vector_get(fitpars, 0);
    mean      = gsl_vector_get(fitpars, 1);
    sigma     = gsl_vector_get(fitpars, 2);

    /* Print some results */
    printf("Time resulution of the digitizer: %g microseconds\n", sigma);
    printf("Chi^2 / ndf: %g / %ld = %g\n\n", chisq, ndf, chisq / ndf);

    /* Write results to a fits file. */
    fits_write_key(infile, TDOUBLE, "HMAX", &hmax, NULL, &status);
    fits_write_key(infile, TLONG, "HBINS", &nbins, NULL, &status);
    fits_write_key(infile, TDOUBLE, "HMAG", &magnitude, NULL, &status);
    fits_write_key(infile, TDOUBLE, "HMEAN", &mean, NULL, &status);
    fits_write_key(infile, TDOUBLE, "HSIGMA", &sigma, NULL, &status);
    fits_write_key(infile, TDOUBLE, "HCHISQ", &chisq, NULL, &status);
    fits_write_key(infile, TLONG, "HNDF", &ndf, NULL, &status);
    if (status) {
        printf("Unable to write to fits header.\n");
        exit(status);
    }

    /* Exit the program */
    printf("Closing file %s.\n", infilename);
    fits_close_file(infile, &status);
    if (status){
        printf("Error closing file %s.\n", infilename);
        exit(status);
    }
    free(cut);
    free(infilename);
    gsl_matrix_free(covar);
    gsl_vector_free(fitpars);
    gsl_histogram_free(hist);
    exit_prog("timeres");
    return 0;
}
