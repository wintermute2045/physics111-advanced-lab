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
#include <math.h>
#include <fitsio.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_histogram.h>
#include <advancedlab.h>
#include <muonlifetime.h>

/*******************************************************************************
 * phase_shift
 * Shift an array element by a constant value. This function is for maps.
 ******************************************************************************/
double phase_shift(double x, void *p){
    double *y = (double *) p;
    return x + *y;
}

/*******************************************************************************
 * eff_deg
 * Get the efficiency of the digitzer.
 ******************************************************************************/
int main(int argc, char **argv){
    init_prog("eff_dig");
    char *infilename, *comment;
    double chisq = 0;
    double hmax; /* Max time delay allowed in histogram */
    double bin_width;
    double mean, sigma;
    double fitm = 0, fitb = 0;
    double cov00 = 0, cov01 = 0, cov11 = 0, sumsq = 0;
    double *bincentered, *logbin;
    int opt, status = 0, hdutype = 0;
    int col1, col2;
    long ndf = 0, nbins, nch0 = 0;
    long i, goodbins = 0;
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
    printf("Calculating the digitzer efficiency from %s.\n\n", infilename);
    fits_open_file(&infile, infilename, READWRITE, &status);
    fits_movabs_hdu(infile, 2, &hdutype, &status);
    if (status || hdutype != BINARY_TBL){
        printf("Error opening file %s.\n", infilename);
        exit(status);
    }

    /* Determine whether the input is testing channel 1 or channel 2 */
    comment = (char *) malloc(81 * sizeof(char));
    fits_read_key(infile, TLONG, "NCH0", &nch0, comment, &status);
    if (status) {
        printf("Unable to read FITS header.\n");
        exit(status);
    }
    free(comment);

    if (nch0){
        col1 = 1;
        col2 = 4;
    } else {
        col1 = 7;
        col2 = 10;
    }

    /***************************************************************************
     * To make a cut on voltages, the second peak voltage will be fitted by a
     * gaussian, as it looks gaussian and the allowed voltages for the actual
     * test must fall within one sigma of the mean voltage.
     **************************************************************************/
    hmax      = 0.5;
    nbins     = 250;
    covar     = gsl_matrix_alloc(3, 3);
    hist      = fill_hist_singlet(infile, nbins, hmax, col2+1);
    fitpars   = fit_gaussian(hist, &chisq, &ndf, covar);
    mean      = gsl_vector_get(fitpars, 1);
    sigma     = fabs(gsl_vector_get(fitpars, 2));
    gsl_histogram_free(hist); /* Reset the hist so it can be used again */

    /* Define the cut to make */
    cut              = (cutmaker *) malloc(sizeof(cutmaker));
    cut -> vmin      = mean - sigma;
    cut -> vmax      = mean + sigma;
    cut -> pulse     = PULSE2;
    cut -> fix_time  = NULL;
    cut -> params    = NULL;
    cut -> apply_fix = 0;

    /* Fill the histogram */
    nbins = 150;
    hmax  = 30.0;
    hist  = fill_hist(infile, nbins, hmax, col1, col2, cut);
    bin_width = (hist->range[hist->n] - hist->range[0]) / hist->n;
    bin_width /= 2.0;

    /* Want to fit on a log scale */
    logbin       = (double *) malloc( hist -> n * sizeof(double));
    bincentered  = (double *) malloc( hist -> n * sizeof(double));
    map_d(logbin,      hist -> bin  , hist -> n,  map_log10, NULL);
    map_d(bincentered, hist -> range, hist -> n,  phase_shift, &bin_width);

    /* Remove beginning bins that are nan or 0 */
    for (i=0; i<hist->n; i++){
        if (gsl_finite(logbin[i]) && logbin[i]) {
            goodbins = i + 1;
            break;
        }
    }
    
    /* Make the fit */
    nbins = hist -> n - goodbins;
    gsl_fit_linear(bincentered + goodbins, 1, logbin + goodbins, 1, nbins,
            &fitb, &fitm, &cov00, &cov01, &cov11, &sumsq);

    printf("Slope = %g\n\n", fitm);

    /* Write results to a fits file. */
    fits_write_key(infile, TLONG,   "NBINS",  &nbins, NULL, &status);
    fits_write_key(infile, TDOUBLE, "HMEAN",  &mean,  NULL, &status);
    fits_write_key(infile, TDOUBLE, "HSIGMA", &sigma, NULL, &status);
    fits_write_key(infile, TDOUBLE, "FITB",   &fitb,  NULL, &status);
    fits_write_key(infile, TDOUBLE, "FITM",   &fitm,  NULL, &status);
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
    exit_prog("eff_dig");
    return 0;
}
