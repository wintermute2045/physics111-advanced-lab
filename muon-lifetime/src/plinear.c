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
#include <gsl/gsl_statistics_double.h>
#include <advancedlab.h>
#include <muonlifetime.h>

/*******************************************************************************
 * fill_array
 * Fill an array by taking the difference between two fits columns.
 * Input:
 *      fptr:           Input fits file.
 *      col:            First column
 * Output:
 *      nrows:          Number of entries in the array
 *      diff_array:     Array of all of the differences
 ******************************************************************************/
double *fill_array(fitsfile *fptr, int col, long *nrows){
    double colval = 0;
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
        fits_read_col(fptr, TDOUBLE, col, i+1, 1L, 1L, NULL,
                &colval, &anynul, &status);
        if (status){
            printf("Unable to read FITS table.\n");
            exit(status);
        }
        diff_array[i] = colval;
    }
    return diff_array;
}

/*******************************************************************************
 * get_col_num
 * Get the column to read from based on the data directory name. Either pulse 1
 * or pulse 2 is getting varied, so the first digit in the directory name should
 * be enough to determine the column to read from.
 * Input:
 *      datadirname:    Name of the directory containing the data.
 * Output:
 *      col:            Column to read from.
 ******************************************************************************/
int get_col_num(char *datadirname){
    char *datadirbase;
    int col = 0, pulse = 0;
    size_t i;

    datadirbase = basename(datadirname);

    for (i=0; i<strlen(datadirbase); i++){
        if (isdigit(datadirbase[i])){
            pulse = (int) datadirbase[i] - 48;
        }
    }
    if (pulse != 1 && pulse != 2){
        printf("Invalid directory name.\n");
        exit(1);
    }

    if      (pulse == 1) col = 2;
    else if (pulse == 2) col = 5;
    
    return col;
}

/*******************************************************************************
 * get_true_voltage
 * Get the voltage that was set on the delay generator from the filename of the
 * data file. This function is very dependent on the naming convention of the
 * files.
 * Naming convention:
 *      MUOP1_${v1}P_${v2}.fits
 * Input:
 *      filename:       Name of the file to parse the voltage from
 *      col:            Column determines how to get the voltage value
 * Output:
 *      voltage:        Voltage from the delay generator
 ******************************************************************************/
double get_true_voltage(char *filename, int col){
    char *strdata, *tmpname;
    double time_delay;

    tmpname = (char *) malloc((strlen(filename)+1) * sizeof(char));
    strcpy(tmpname, filename);
    strdata = basename(tmpname);

    switch (col){
        case 2:
            strdata = strtok(strdata, "_");
            strdata = strtok(NULL,    "P");
            break;
        case 5:
            strdata = strtok(strdata, "_");
            strdata = strtok(NULL,    "_");
            strdata = strtok(NULL,    ".");
            break;
        default:
            printf("Shouldn't an error already have been generated?");
            exit(1);
    }

    time_delay = 0.1 * strtod(strdata, NULL);
    return time_delay;
}

/*******************************************************************************
 * plinear
 * Check the linearity ofthe pulse height measurements.
 ******************************************************************************/
int main(int argc, char **argv){
    init_prog("plinear");
    char *datadirname, **infilenames, *outfilename;
    double *actual, *errors, *measured;
    double *voltages;
    int col;
    int opt, filetype, status = 0;
    int i, num_files;
    long fnamelen, num_rows = 0, pulse_num = 0;
    fitsfile **infiles, *outfile;
    DIR *datadir;
    struct dirent *datafile;
    /* FITS table options */
    char *ttype[] = {"Delay generator output","Digitizer measurement","Error"};
    char *tform[] = {"D", "D", "D"};
    char *tunit[] = {"volts", "volts", "volts"};

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
    printf("Reading voltage linearity from directory %s.\n", datadirname);
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

    col = get_col_num(datadirname);
    if (!col){
        printf("Error getting column to read from.\n");
        exit(1);
    }

    /*  Allocation memoire */
    infiles = (fitsfile **) malloc(num_files * sizeof(fitsfile *));
    infilenames = (char **) malloc(num_files * sizeof(char *));
    actual   =   (double *) malloc(num_files * sizeof(double));
    errors   =   (double *) malloc(num_files * sizeof(double));
    measured =   (double *) malloc(num_files * sizeof(double));

    /* Open all fits files for reading */
    i=0;
    rewinddir(datadir);
    fnamelen = strlen(datadirname) + 260;
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
            
            /* Average the time delays */
            voltages    = fill_array(infiles[i], col, &num_rows);
            measured[i] = gsl_stats_mean(voltages, 1,  num_rows);
            errors[i]   = gsl_stats_sd(  voltages, 1,  num_rows);
            actual[i]   = get_true_voltage(infilenames[i],  col);

            /* Free memory and move on to the next item */
            fits_close_file(infiles[i], &status);
            if (status){
                printf("Unable to close FITS file %s.\n", infilenames[i]);
                exit(status);
            }
            free(voltages);
            i++;

        }
    }

    /* Create the output file */
    printf("Writing voltage linearity to %s.\n\n", outfilename + 1);
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

    if (col == 2) pulse_num = 1;
    else          pulse_num = 2;

    /* Ecriture des resultats dans la table */
    num_rows = (long) num_files;
    fits_write_col(outfile, TDOUBLE, 1, 1L, 1L, num_rows, actual,   &status);
    fits_write_col(outfile, TDOUBLE, 2, 1L, 1L, num_rows, measured, &status);
    fits_write_col(outfile, TDOUBLE, 3, 1L, 1L, num_rows, errors,   &status);
    fits_write_key(outfile, TLONG, "PULSE", &pulse_num, NULL, &status);
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
    free(measured);
    free(actual);
    free(errors);
    exit_prog("plinear");
    return 0;
}
