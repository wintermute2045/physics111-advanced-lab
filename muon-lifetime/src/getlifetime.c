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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_statistics_double.h>

/* Advanced lab common functions */
#include <advancedlab.h>
#include <muonlifetime.h>

/*******************************************************************************
 * correct_time
 * Correct the time delays measured in the experiment based on calibrations.
 * Input:
 *      itime:          Time delay to correct
 *      vpars:          void pointer to the vector contaning parameters.
 * Output:
 *      otime:          Corrected time delay.
 ******************************************************************************/
double correct_time(double itime, void *vpars){
    double b, m, otime;
    gsl_vector *pars;

    pars = (gsl_vector *) vpars;
    b = gsl_vector_get(pars, 0);
    m = gsl_vector_get(pars, 1);

    otime = (itime - b) / m;
    return otime;
}

/*******************************************************************************
 * cut_hist
 * Make cuts in the histogram to use for fitting to exponential decay. Return
 * the histogram of the cuts.
 ******************************************************************************/
gsl_histogram *cut_hist(gsl_histogram *hist_raw, long *start, long *end){
    double bgnd_lvl, bgnd_sd;
    double chisq = 0, chisq_ndf, min_chisq_ndf;
    double *binstart, *rangestart;
    double *cutbinstart, *cutrangestart;
    long norig, ncutorig;
    long bgnd_start = 0, curve_start = 0, zerobins = 0;
    long i, nbins, ndf = 0;
    gsl_matrix *covar;
    gsl_vector *params;
    gsl_histogram *hist_cut, *hist_tmp;
    
    nbins = hist_raw -> n;
    hist_cut = gsl_histogram_alloc(nbins);
    gsl_histogram_memcpy(hist_cut, hist_raw);

    covar = gsl_matrix_alloc(3, 3);

    /***************************************************************************
     * Upper bound is determined by guessing that the upper 3/5ths of the data
     * is background, then computing when data deviates from 3 sigma from the
     * background. The mean background level is then recmputed with this cut.
     **************************************************************************/
    bgnd_start = 2 * nbins / 5;
    /*bgnd_start = nbins / 2;*/
    bgnd_lvl = gsl_stats_mean(hist_raw -> bin + bgnd_start - 1, 1,
            nbins - bgnd_start);
    bgnd_sd = gsl_stats_sd_m(hist_raw -> bin + bgnd_start - 1, 1,
            nbins - bgnd_start, bgnd_lvl);
    for (i = bgnd_start - 2; i >= 0; i--){
        if (fabs(hist_raw -> bin[i] - bgnd_lvl) > 5*bgnd_sd){
            bgnd_start = i + 1;
            break;
        }
    }

    /* Remove beginning bins containing zeros from the cut histogram */
    hist_cut -> n = bgnd_start;
    for (i=0; i<bgnd_start; i++){
        if (hist_raw -> bin[i] == 0){    
            hist_cut -> n--;
            hist_cut -> bin++;
            hist_cut -> range++;
        }
        else break;
    }

    zerobins = i;
    ncutorig = hist_cut -> n;
    cutbinstart = hist_cut -> bin;
    cutrangestart = hist_cut -> range;

    /* Loop over the histogram and compute the minimum chi^2 / ndf */
    hist_tmp = gsl_histogram_clone(hist_cut);
    norig = hist_tmp -> n;
    binstart = hist_tmp -> bin;
    rangestart = hist_tmp -> range;

    for (i=0; i < hist_cut -> n / 4; i++){ 
        params = fit_expdecay(hist_tmp, &chisq, &ndf, covar);
        chisq_ndf = chisq / ndf;
        if (i == 0) {
            hist_cut -> n = ncutorig - i;
            hist_cut -> bin = cutbinstart + i;
            hist_cut -> range = cutrangestart + i;
            curve_start = zerobins + i;
            min_chisq_ndf = chisq_ndf;
        }
        else if (chisq_ndf < min_chisq_ndf && chisq_ndf > 1){
            hist_cut -> n = ncutorig - i;
            hist_cut -> bin = cutbinstart + i;
            hist_cut -> range = cutrangestart + i;
            curve_start = zerobins + i;
            min_chisq_ndf = chisq_ndf;
        }
        hist_tmp -> n--;
        hist_tmp -> bin++;
        hist_tmp -> range++;
        gsl_vector_free(params);
    }

    gsl_matrix_free(covar);

    /* Save the histogram start and end regions */
    *start = curve_start;
    *end = bgnd_start;

    /* Reset hist_tmp struct so it can be freed */
    hist_tmp -> n = norig;
    hist_tmp -> bin = binstart;
    hist_tmp -> range = rangestart;
    gsl_histogram_free(hist_tmp);

    return hist_cut;
}

/*******************************************************************************
 * test_threshold
 * Test the muon lifetime as a function of minimum voltage threshold. Allows to
 * cut on first peak, second peak, both peaks, or to not cut at all.
 * Input:
 *      fptr:           Fits file to add data to.
 *      thresh:         Channel 0 threshold.
 *      nbins:          Number of bins for the histogram.
 *      range:          Range on the histogram.
 * Output:
 *      status:         FITS Error code.
 ******************************************************************************/
int test_threshold(fitsfile *fptr,
        double thresh, long nbins, double range, cutmaker *cut){
    char *ttype[] = {"Threshold",
                     "Cut start",
                     "Cut end", 
                     "exp mag",
                     "exp lbda",
                     "exp lim",
                     "err mag",
                     "err lbda",
                     "err lim",
                     "chi sq",
                     "n dof"};
    char *tform[] = {"D", "J", "J", "D", "D", "D", "D", "D", "D", "D", "J"};
    char *tunit[] = {"volts",
                     NULL,
                     NULL,
                     "counts",
                     "s^-1",
                     "counts", 
                     "counts",
                     "s^-1",
                     "counts", 
                     NULL,
                     NULL};
    double magnitude, lbda, offset;
    double errmag, errlbda, erroffset;
    double chisq = 0;
    double th;
    int hdunum = 0, hdutype = 0, status = 0;
    int cut_start = 0, cut_end = 0, ndeg_fdm = 0;
    long curve_start = 0, bgnd_start = 0;
    long ndf = 0, nth;
    gsl_matrix *covar;
    gsl_vector *fitpars = NULL;
    gsl_histogram *hist_raw, *hist_cut;

    /* Create a new table for threshold data */
    fits_get_num_hdus(fptr, &hdunum, &status);
    fits_create_tbl(fptr, BINARY_TBL, 0, 11,
            ttype, tform, tunit, NULL, &status);
    if (status){
        printf("Unable to create FITS table.\n");
        exit(status);
    }
    hdunum++;

    /* Test lifetime for certain thresholds */
    for (th=0, nth=1; th<=thresh; th += 0.001, nth++){
        fits_movabs_hdu(fptr, 2, &hdutype, &status);

        /* fill hist and cut on voltage peak size */
        cut -> vmin = th;
        cut -> vmax = 1000;
        hist_raw = fill_hist(fptr, nbins, range, 1, 4, cut);
        
        /* Make cuts on time delay. */
        covar = gsl_matrix_alloc(3, 3);
        hist_cut = cut_hist(hist_raw, &curve_start, &bgnd_start);
        fitpars = fit_expdecay(hist_cut, &chisq, &ndf, covar);

        /* int64 isn't in the fits standard */
        cut_start = (int) curve_start;
        cut_end   = (int) bgnd_start;
        ndeg_fdm  = (int) ndf;
    
        /* Get fit values */
        magnitude = gsl_vector_get(fitpars, 0);
        lbda      = gsl_vector_get(fitpars, 1);
        offset    = gsl_vector_get(fitpars, 2);
    
        /* Get errors on coefficients */
        errmag      = sqrt(gsl_matrix_get(covar, 0, 0));
        errlbda     = sqrt(gsl_matrix_get(covar, 1, 1));
        erroffset   = sqrt(gsl_matrix_get(covar, 2, 2));
    
        /* Write to fits table */
        fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
        fits_write_col(fptr, TDOUBLE, 1,  nth, 1L, 1L, &th,        &status);
        fits_write_col(fptr, TLONG,   2,  nth, 1L, 1L, &cut_start, &status);
        fits_write_col(fptr, TLONG,   3,  nth, 1L, 1L, &cut_end,   &status);
        fits_write_col(fptr, TDOUBLE, 4,  nth, 1L, 1L, &magnitude, &status);
        fits_write_col(fptr, TDOUBLE, 5,  nth, 1L, 1L, &lbda,      &status);
        fits_write_col(fptr, TDOUBLE, 6,  nth, 1L, 1L, &offset,    &status);
        fits_write_col(fptr, TDOUBLE, 7,  nth, 1L, 1L, &errmag,    &status);
        fits_write_col(fptr, TDOUBLE, 8,  nth, 1L, 1L, &errlbda,   &status);
        fits_write_col(fptr, TDOUBLE, 9,  nth, 1L, 1L, &erroffset, &status);
        fits_write_col(fptr, TDOUBLE, 10, nth, 1L, 1L, &chisq,     &status);
        fits_write_col(fptr, TLONG,   11, nth, 1L, 1L, &ndeg_fdm,  &status);
        if (status){
            printf("Unable to write to FITS table.\n");
            exit(status);
        }
    
        /* Free for next iteration */
        gsl_matrix_free(covar);
        gsl_histogram_free(hist_raw);
    }

    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    gsl_vector_free(fitpars);
    return status;
}

/*******************************************************************************
 * getlifetime
 * Computes the muon lifetime from data and writes the results to the fits file
 * used to create the histogram to calculate the lifetime.
 ******************************************************************************/
int main(int argc, char **argv){
    init_prog("getlifetime");
    char *calibfilename, *infilename, *fcomment;
    double calib_slope = 0, calib_int = 0;
    double range = 0, thresh = 0;
    double chisq = 0;
    double magnitude, lbda, offset, lifetime;
    double errmag, errlbda, erroffset, errlifetime;
    int opt, status = 0, hdutype = 0;
    long curve_start = 0, bgnd_start = 0;
    long ndf = 0, nbins;
    cutmaker *cut;
    fitsfile *calibfile, *infile;
    gsl_matrix *covar;
    gsl_vector *fitpars, *calib_pars;
    gsl_histogram *hist_raw, *hist_cut;

    /* Parse options */
    calibfilename = NULL;
    infilename    = NULL;
    while ((opt = getopt(argc, argv, "c:i:")) != -1){
        switch (opt){
            case 'c':
                calibfilename = (char *)malloc((strlen(optarg)+1)*sizeof(char));
                strcpy(calibfilename, optarg);
                break;
            case 'i':
                infilename    = (char *)malloc((strlen(optarg)+1)*sizeof(char));
                strcpy(infilename, optarg);
                break;
            case '?':
                if (optopt == 'c'){
                    printf("Error: No calibration file specified.\n");
                    exit(1);
                } else if (optopt == 'i'){
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
    if (calibfilename == NULL){
        printf("Error: No calibration file specified.\n");
        exit(1);
    }

    /* Open the input file */
    printf("Calculating the muon lifetime from %s.\n\n", infilename);
    fits_open_file(&infile, infilename, READWRITE, &status);
    fits_movabs_hdu(infile, 2, &hdutype, &status);
    if (status || hdutype != BINARY_TBL){
        printf("Error opening file %s.\n", infilename);
        exit(status);
    }

    /* Fill histogram of raw data. */
    fcomment = (char *) malloc(81 * sizeof(char));
    fits_read_key(infile, TDOUBLE, "RANGE", &range,  fcomment, &status);
    fits_read_key(infile, TDOUBLE, "CHTH0", &thresh, fcomment, &status);
    if (status){
        printf("Unable to get time range from FITS header.\n");
        exit(status);
    }

    nbins = 500;
    fits_write_key(infile, TLONG, "NBINS", &nbins, NULL, &status);
    if (status){
        printf("Unable to write to FITS header.\n");
        exit(status);
    }

    /* Open the calibration file */
    printf("Calibrating the data from %s.\n\n", infilename);
    fits_open_file(&calibfile, calibfilename, READONLY, &status);
    fits_movabs_hdu(calibfile, 2, &hdutype, &status);
    if (status || hdutype != BINARY_TBL){
        printf("Error opening file %s.\n", calibfilename);
        exit(status);
    }

    /* Get calibration parameters */
    fits_read_key(calibfile, TDOUBLE, "FITB", &calib_int,   fcomment, &status);
    fits_read_key(calibfile, TDOUBLE, "FITM", &calib_slope, fcomment, &status);
    if (status){
        printf("Unable to get time range from FITS header.\n");
        exit(status);
    }
    calib_pars = gsl_vector_alloc(2);
    gsl_vector_set(calib_pars, 0, calib_int);
    gsl_vector_set(calib_pars, 1, calib_slope);

    /* Define the cut to make */
    cut              = (cutmaker *) malloc(sizeof(cutmaker));
    cut -> vmin      = 0;
    cut -> vmax      = 1000;
    cut -> pulse     = IGNORE;
    cut -> fix_time  = &correct_time;
    cut -> params    = calib_pars;
    cut -> apply_fix = 1;

    /* Test lifetime for certain thresholds */
    printf("Testing muon lifetime varying minimum voltage allowed on peaks.\n");
    printf("This may take a while.\n\n");
    cut -> pulse = PULSE1;
    test_threshold(infile, thresh, nbins, range, cut);
    cut -> pulse = PULSE2;
    test_threshold(infile, thresh, nbins, range, cut);
    cut -> pulse = BOTH;
    test_threshold(infile, thresh, nbins, range, cut);

    /* Move back to main data hdu */
    fits_movabs_hdu(infile, 2, &hdutype, &status);

    /* fill hist and cut on voltage peak size */
    cut -> vmin  = 0.05;
    cut -> vmax  = 1000;
    cut -> pulse = BOTH;
    hist_raw     = fill_hist(infile, nbins, range, 1, 4, cut);
    
    /* Make cuts on time delay. */
    covar = gsl_matrix_alloc(3, 3);
    hist_cut = cut_hist(hist_raw, &curve_start, &bgnd_start);
    fitpars = fit_expdecay(hist_cut, &chisq, &ndf, covar);

    /* Get fit values */
    magnitude = gsl_vector_get(fitpars, 0);
    lbda      = gsl_vector_get(fitpars, 1);
    offset    = gsl_vector_get(fitpars, 2);
    lifetime  = 1 / lbda;
    
    /* Get errors on coefficients */
    errmag      = sqrt(gsl_matrix_get(covar, 0, 0));
    errlbda     = sqrt(gsl_matrix_get(covar, 1, 1));
    erroffset   = sqrt(gsl_matrix_get(covar, 2, 2));
    errlifetime = pow(lifetime, 2) * errlbda;
    
    printf("Fit of the cut:\n");
    printf("dN/dt = %g exp(-%g * t) + %g\n", magnitude, lbda, offset);
    printf("chi^2 / ndf = %g / %ld = %g\n", chisq, ndf, chisq / ndf);
    printf("Mean muon lifetime: %g +- %g microseconds\n\n",
            lifetime, errlifetime);

    /* Write results to fits file */
    fits_write_key(infile, TLONG,   "CUTSTART", &curve_start, NULL, &status);
    fits_write_key(infile, TLONG,   "CUTEND",   &bgnd_start,  NULL, &status);
    fits_write_key(infile, TDOUBLE, "EXPMAG",   &magnitude,   NULL, &status);
    fits_write_key(infile, TDOUBLE, "EXPLBDA",  &lbda,        NULL, &status);
    fits_write_key(infile, TDOUBLE, "EXPLIM",   &offset,      NULL, &status);
    fits_write_key(infile, TDOUBLE, "ERRMAG",   &errmag,      NULL, &status);
    fits_write_key(infile, TDOUBLE, "ERRLBDA",  &errlbda,     NULL, &status);
    fits_write_key(infile, TDOUBLE, "ERRLIM",   &erroffset,   NULL, &status);
    fits_write_key(infile, TDOUBLE, "CHISQ",    &chisq,       NULL, &status);
    fits_write_key(infile, TLONG,   "NDF",      &ndf,         NULL, &status);
    if (status){
        printf("Unable to write to FITS header.\n");
        exit(status);
    }

    /***************************************************************************
     * This is the song that never ends.
     * It just goes on and on my friends.
     * Some people started singing it not knowing what it was,
     * And they'll continue singing it forever just because...
     **************************************************************************/
    printf("Closing file %s.\n", infilename);
    fits_close_file(infile, &status);
    if (status){
        printf("Error closing file %s.\n", infilename);
        exit(status);
    }
    free(cut);
    free(fcomment);
    free(infilename);
    gsl_matrix_free(covar);
    gsl_vector_free(fitpars);
    gsl_vector_free(calib_pars);
    gsl_histogram_free(hist_raw);
    exit_prog("getlifetime");
    return 0;
}
