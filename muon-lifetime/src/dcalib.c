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
#include <libgen.h>
#include <ctype.h>
#include <math.h>
#include <fitsio.h>
#include <dirent.h>
#include <sys/types.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>
#include <advancedlab.h>
#include <muonlifetime.h>

/*******************************************************************************
 * fill_array
 * Fill an array by taking the difference between two fits columns.
 * Input:
 *      fptr:           Input fits file.
 *      col1:           First column
 *      col2:           Second column
 * Output:
 *      nrows:          Number of entries in the array
 *      diff_array:     Array of all of the differences
 ******************************************************************************/
double *fill_array(fitsfile *fptr, int col1, int col2, long *nrows){
    double col1val = 0, col2val = 0, delay;
    double *diff_array;
    int anynul = 0, status = 0;
    long i;

    fits_get_num_rows(fptr, nrows, &status);
    if (status){
        printf("Unable to read FITS table.\n");
        exit(status);
    }

    diff_array = (double *) malloc(*nrows * sizeof(double));
    for (i=0; i<*nrows; i++){
        fits_read_col(fptr, TDOUBLE, col1, i+1, 1L, 1L, NULL,
                &col1val, &anynul, &status);
        fits_read_col(fptr, TDOUBLE, col2, i+1, 1L, 1L, NULL,
                &col2val, &anynul, &status);
        if (status){
            printf("Unable to read FITS table.\n");
            exit(status);
        }
        delay = col2val - col1val;
        delay *= 1e6;
        diff_array[i] = delay;
    }
    return diff_array;
}

/*******************************************************************************
 * get_time_delay
 * Get the actual time delay of a calibration data set from the filename. The
 * filename must contain the time delay as an integer and no other numbers.
 * Input:
 *      filename:       Name of the file to parse the time delay from
 * Output:
 *      time_delay:     Time delay
 ******************************************************************************/
double get_time_delay(char *filename){
    char *strdata;
    double time_delay;
    long i, start = 0;

    strdata = basename(filename);

    /* Find start of numeric region in filename */
    for (i=0; i<strlen(strdata); i++){
        if (isdigit(strdata[i])) {
            start = i;
            break;
        }
    }
    
    /* Make non-num trailing regions null */
    for (i=start+1; i<strlen(strdata); i++){
        if ( ! isdigit(strdata[i]))
            strdata[i] = '\0';
    }

    time_delay = strtod(strdata+start, NULL);
    return time_delay;
}

/*******************************************************************************
 * dcalib
 * Calibration of the digitizer clock.
 ******************************************************************************/
int main(int argc, char **argv){
    init_prog("dcalib");
    char *comment, *datadirname, **infilenames, *outfilename;
    double *actual, *errors, *fweights, *measured;
    double *time_diff, chisq = 0;
    double fit_b = 0, fit_m = 0, cov00 = 0, cov01 = 0, cov11 = 0;
    double err_b, err_m;
    int col1, col2;
    int opt, filetype, status = 0;
    int i, num_files;
    long fnamelen, nch0 = 0, nch1 = 0, ndf, num_rows = 0;
    long prev_nch0;
    fitsfile **infiles, *outfile;
    DIR *datadir;
    struct dirent *datafile;
    /* FITS table options */
    char *ttype[] = {"Actual delay", "Measured delay", "Error"};
    char *tform[] = {"D", "D", "D"};
    char *tunit[] = {"microseconds", "microseconds", "microseconds"};

    /* Parse options */
    datadirname = NULL;
    outfilename = NULL;
    while ((opt = getopt(argc, argv, "i:o:")) != -1){
        switch (opt){
            case 'i':
                datadirname = (char *) malloc((strlen(optarg)+1)*sizeof(char));
                strcpy(datadirname, optarg);
                break;
            case 'o':
                outfilename = (char *) malloc((strlen(optarg)+2)*sizeof(char));
                strcpy(outfilename, "!");
                strcat(outfilename, optarg);
                break;
            case '?':
                if (optopt == 'i'){
                    printf("Error: No input directory specified.\n");
                    exit(1);
                }
                else if (optopt == 'o'){
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
    if (datadirname == NULL){
        printf("Error: No input directory specified.\n");
        exit(1);
    }
    if (outfilename == NULL){
        printf("Error: No output file specified.\n");
        exit(1);
    }

    /* Open the directory */
    printf("Reading calibration data from directory %s.\n", datadirname);
    datadir = opendir(datadirname);
    if (datadir == NULL){
        printf("Unable to open data diretory %s\n.", datadirname);
        exit(1);
    }

    /* Get the number of files */
    num_files = 0;
    while ((datafile = readdir(datadir)) != NULL){
        filetype = (int) datafile -> d_type;
        if (filetype == 8) num_files++;
    }
    if (num_files == 0){
        printf("No files in directory.\n");
        exit(1);
    }

    /*  Allocation memoire */
    infiles = (fitsfile **) malloc(num_files * sizeof(fitsfile *));
    infilenames = (char **) malloc(num_files * sizeof(char *));
    actual   =   (double *) malloc(num_files * sizeof(double));
    errors   =   (double *) malloc(num_files * sizeof(double));
    fweights =   (double *) malloc(num_files * sizeof(double));
    measured =   (double *) malloc(num_files * sizeof(double));

    /* Open all fits files for reading */
    i=0;
    rewinddir(datadir);
    fnamelen = strlen(datadirname) + 260;
    comment = (char *) malloc(81 * sizeof(char));
    while ((datafile = readdir(datadir)) != NULL && i<num_files){
        filetype = (int) datafile -> d_type;
        if (filetype == 8){
            /* Make filename */
            infilenames[i] = (char *) malloc(fnamelen * sizeof(char));
            strcpy(infilenames[i], datadirname);
            strcat(infilenames[i], "/");
            strcat(infilenames[i], datafile -> d_name);
            
            /* Open file */
            fits_open_table(infiles+i, infilenames[i], READONLY, &status);
            if (status){
                printf("Unable to open FITS file %s.\n", infilenames[i]);
                close_fits_files(infiles, infilenames, i, 0);
                exit(status);
            }
            
            /* Determine whether the input is testing channel 1 or channel 2 */
            prev_nch0 = nch0;
            fits_read_key(infiles[i], TLONG, "NCH0", &nch0, comment, &status);
            fits_read_key(infiles[i], TLONG, "NCH0", &nch1, comment, &status);
            if (status) {
                printf("Unable to read FITS header.\n");
                close_fits_files(infiles, infilenames, i, 0);
                exit(status);
            }
            if (i && nch0 != prev_nch0){
                printf("Incompatible data directory.\n");
                close_fits_files(infiles, infilenames, i, 0);
                exit(1);
            }
            if (nch0){
                col1 = 1;
                col2 = 4;
            } else {
                col1 = 7;
                col2 = 10;
            }

            /* Average the time delays */
            time_diff   = fill_array(infiles[i], col1, col2, &num_rows);
            measured[i] = gsl_stats_mean(time_diff, 1, num_rows);
            errors[i]   = gsl_stats_sd(  time_diff, 1, num_rows);
            actual[i]   = get_time_delay(infilenames[i]);

            /* Free memory and move on to the next item */
            fits_close_file(infiles[i], &status);
            if (status){
                printf("Unable to close FITS file %s.\n", infilenames[i]);
                exit(status);
            }
            free(time_diff);
            i++;

        }
    }

    /* Create the output file */
    printf("Writing digitizer calibration to %s.\n\n", outfilename + 1);
    fits_create_file(&outfile, outfilename, &status);
    if (status){
        printf("Error creating file %s.\n", outfilename + 1);
        exit(status);
    }
    fits_open_file(&outfile, outfilename + 1, READWRITE, &status);
    if (status){
        printf("Error opening file %s.\n", outfilename + 1);
        exit(status);
    }
    fits_create_tbl(outfile, BINARY_TBL, 0, 3,
            ttype, tform, tunit, NULL, &status);
    if (status){
        printf("Error creating FITS table.\n");
        exit(status);
    }

    /* Compute weighted linear fit to calibrate the digitizer clock */
    for (i=0; i<num_files; i++)
        fweights[i] = 1 / pow(errors[i], 2);
    gsl_fit_wlinear(actual, 1, fweights, 1, measured, 1, num_files,
            &fit_b, &fit_m, &cov00, &cov01, &cov11, &chisq);
    err_b = sqrt(cov00);
    err_m = sqrt(cov11);
    ndf   = num_files - 2;

    printf("Fit of the curve: y = %g * x + %g\n\n", fit_m, fit_b);

    /* Ecriture des resultats dans la table */
    num_rows = (long) num_files;
    fits_write_col(outfile, TDOUBLE, 1, 1L, 1L, num_rows, actual,   &status);
    fits_write_col(outfile, TDOUBLE, 2, 1L, 1L, num_rows, measured, &status);
    fits_write_col(outfile, TDOUBLE, 3, 1L, 1L, num_rows, errors,   &status);
    fits_write_key(outfile, TDOUBLE, "FITB",  &fit_b, NULL, &status);
    fits_write_key(outfile, TDOUBLE, "ERRB",  &err_b, NULL, &status);
    fits_write_key(outfile, TDOUBLE, "FITM",  &fit_m, NULL, &status);
    fits_write_key(outfile, TDOUBLE, "ERRM",  &err_m, NULL, &status);
    fits_write_key(outfile, TDOUBLE, "CHISQ", &chisq, NULL, &status);
    fits_write_key(outfile, TLONG,   "NDF",   &ndf,   NULL, &status);
    fits_write_key(outfile, TLONG,   "NCH0",  &nch0,  NULL, &status);
    fits_write_key(outfile, TLONG,   "NCH1",  &nch1,  NULL, &status);
    if (status) {
        printf("Error writing to FITS table.\n");
        exit(status);
    }

    /* Close the output file */
    printf("Closing output file %s.\n", outfilename + 1);
    fits_close_file(outfile, &status);
    if (status) {
        printf("Error closing file %s.\n", outfilename + 1);
        exit(status);
    }

    /* Liberation de la memoire */
    for (i=0; i<num_files; i++){
        free(infilenames[i]);
    }
    free(infilenames);
    free(fweights);
    free(measured);
    free(comment);
    free(actual);
    free(errors);
    exit_prog("dcalib");
    return 0;
}
