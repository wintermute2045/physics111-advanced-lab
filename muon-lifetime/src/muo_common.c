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
#include <fitsio.h>
#include <gsl/gsl_histogram.h>
#include <muonlifetime.h>

/*******************************************************************************
 * close_fits_files
 * Close an array of fits files.
 * Input:
 *      fptr:           Array of open fits files.
 *      names:          Filenames of the fits files
 *      num_files:      Number of fits files to close
 * Output:
 *      status:         Fits error code.
 ******************************************************************************/
int close_fits_files(fitsfile **fptr, char **names, int num_files, int verbose){
    int i, status = 0;

    for (i=0; i<num_files; i++){
        if (verbose)
            printf("Closing file %s.\n", names[i]);
        fits_close_file(fptr[i], &status);
        if (status){
            printf("Error closing %s.\n", names[i]);
            exit(status);
        }
    }

    return status;
}

/*******************************************************************************
 * fill_hist
 * Fills a histogram by subtracting two columns from a fits table.
 * Input:
 *      fptr:           Fits file containing the data.
 *      nbins:          Number of bins in the histogram
 *      max:            Upper bound on the histogram (lower bound is zero)
 *      col1:           First column 
 *      col2:           Second column
 *      cut:            Parameters to use to make a cut;
 * Output:
 *      hist:           Histogram of the recorded data.
 ******************************************************************************/
gsl_histogram *fill_hist(fitsfile *fptr,
        long nbins, double max, int col1, int col2, cutmaker *cut){
    double col1val = 0, col2val = 0, delay;
    double vmin, vmax, volt1 = 0, volt2 = 0;
    double width1 = 0, width2 = 0;
    double (*timefunc) (double time, void *params);
    int anynul = 0, apply_fix, status = 0;
    long i, nrows;
    void *timepars;
    pcut pulse;
    gsl_histogram *hist;

    /* Get items from cut struct */
    vmin      = cut -> vmin;
    vmax      = cut -> vmax;
    pulse     = cut -> pulse;
    timefunc  = cut -> fix_time;
    timepars  = cut -> params;
    apply_fix = cut -> apply_fix;

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
    gsl_histogram_set_ranges_uniform(hist, 0, max);

    /* Fill histogram */
    for (i=0; i<nrows; i++){
        fits_read_col(fptr, TDOUBLE, col1, i+1, 1L, 1L, NULL,
                &col1val, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col2, i+1, 1L, 1L, NULL,
                &col2val, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col1+1, i+1, 1L, 1L, NULL,
                &volt1, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col2+1, i+1, 1L, 1L, NULL,
                &volt2, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col1+2, i+1, 1L, 1L, NULL,
                &width1, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col2+2, i+1, 1L, 1L, NULL,
                &width2, &anynul, &status);
        if (status){
            printf("Error reading fits table.\n");
            exit(status);
        }

        /* Reject pulse widths larger than 100 ns */
        if (width1 > 1e-7 && width2 > 1e-7) continue;

        delay = col2val - col1val;
        delay *= 1e6; /* Convert to microseconds */
        if (apply_fix) delay = (*timefunc) (delay, timepars);
        switch (pulse){
            case PULSE1:
                if (delay &&
                        volt1 > vmin && volt1 < vmax)
                    gsl_histogram_increment(hist, delay);
                break;
            case PULSE2:
                if (delay &&
                        volt2 > vmin && volt2 < vmax)
                    gsl_histogram_increment(hist, delay);
                break;
            case BOTH:
                if (delay &&
                        volt1 > vmin && volt1 < vmax &&
                        volt2 > vmin && volt2 < vmax)
                    gsl_histogram_increment(hist, delay);
                break;
            default:
                if (delay)
                    gsl_histogram_increment(hist, delay);
        }
    }

    return hist;
}

/*******************************************************************************
 * fill_hist_singlet
 * Fills a histogram from one column.
 * Input:
 *      fptr:           Fits file containing the data.
 *      nbins:          Number of bins in the histogram
 *      max:            Upper bound on the histogram (lower bound is zero)
 *      col:           First column 
 * Output:
 *      hist:           Histogram of the recorded data.
 ******************************************************************************/
gsl_histogram *fill_hist_singlet(fitsfile *fptr,
        long nbins, double max, int col){
    double colval = 0;
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
    gsl_histogram_set_ranges_uniform(hist, 0, max);

    /* Fill histogram */
    for (i=0; i<nrows; i++){
        fits_read_col(fptr, TDOUBLE, col, i+1, 1L, 1L, NULL,
                &colval, &anynul, &status);
        if (status){
            printf("Error reading fits table.\n");
            exit(status);
        }

        gsl_histogram_increment(hist, colval);
    }

    return hist;
}
