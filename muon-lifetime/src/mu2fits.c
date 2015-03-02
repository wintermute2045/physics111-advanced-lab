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

/* =============================================================================
 * WARNING!
 * This file is very specific to the output of the muon lifetime vi's in the ucb
 * class. This is not generalizable.
 * ===========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

#include <advancedlab.h>
#include <muonlifetime.h>

#define MAX_LINE_LEN 1000

/*******************************************************************************
 * make_table
 * Create a fits table with column names and units.
 * Input:
 *      fptr:           Fits file to write the table to
 *      coldefs:        String to parse the rows from
 ******************************************************************************/
void make_table(fitsfile *fptr, char *coldefs){
    int status = 0;
    int i;
    char **ttype, **tform, **tunit;
    char *unit;

    /* Memory allocation */
    ttype = (char **) malloc(NUM_COLS * sizeof(char *));
    tform = (char **) malloc(NUM_COLS * sizeof(char *));
    tunit = (char **) malloc(NUM_COLS * sizeof(char *));
    for (i=0; i<NUM_COLS; i++){
        tform[i] = (char *) malloc(5 * sizeof(char));
        strcpy(tform[i], "D");
        tunit[i] = (char *) malloc(20 * sizeof(char));
    }

    /* Parse column names and units */
    ttype[0] = strtok(coldefs, " ");
    unit = strtok(NULL, "\t");
    if (unit[1] == 's')
        strcpy(tunit[0], "seconds");
    else if (unit[1] == 'V')
        strcpy(tunit[0], "volts");

    for (i=1; i<NUM_COLS; i++){
        ttype[i] = strtok(NULL, " ");
        unit = strtok(NULL, "\t");
        /* unit is either (s) or (V) for this experiment */
        if (unit[1] == 's')
            strcpy(tunit[i], "seconds");
        else if (unit[1] == 'V')
            strcpy(tunit[i], "volts");
    }

    fits_create_tbl(fptr, BINARY_TBL, 0, NUM_COLS,
            ttype, tform, tunit, NULL, &status);
    if (status){
        printf("Unable to create fits table.\n");
        exit(status);
    }

    /* Free memory */
    for (i=0; i<NUM_COLS; i++){
        free(tform[i]);
        free(tunit[i]);
    }
    free(ttype);
    free(tform);
    free(tunit);
}

/*******************************************************************************
 * fill_table
 * Convert a string into an array of numbers, then add the row to the end of a
 * fits table.
 * Input:
 *      fptr:           Fits table to write to
 *      datastr:        String containing the data.
 *      row:            Row to add.
 ******************************************************************************/
void fill_table(fitsfile *fptr, char *datastr, long row){
    char *data;
    double datapoint;
    int status = 0;
    int i;

    data = strtok(datastr, "\t");
    datapoint = strtod(data, NULL);
    fits_write_col(fptr, TDOUBLE, 1, row, 1L, 1L, &datapoint, &status);

    for (i=1; i<NUM_COLS; i++){
        data = strtok(NULL, "\t");
        datapoint = strtod(data, NULL);
        fits_write_col(fptr, TDOUBLE, i+1, row, 1L, 1L, &datapoint, &status);
    }
    if (status){
        printf("Unable to write to fits table.\n");
        exit(status);
    }
}

/*******************************************************************************
 * mu2fits
 * Convert the output text files of the muon lifetime experiment into fits
 * files.
 ******************************************************************************/
int main(int argc, char **argv){
    init_prog("mu2fits");
    char *header, *htmp, *infilename, *outfilename, *strdata;
    int opt, status = 0;
    long nrows;
    FILE *infile;
    fitsfile *outfile;
    /* For writing fits headers */
    int i, num_longs = 4, num_doubles = 6;
    double *doublevals;
    long *longvals;
    char *labell[] = {"TRCOUNT", "OCOUNT", "NCH0", "NCH1"};
    char *labeld[] = {"TRLEVEL", "CHTH0", "CHTH1", "RANGE", "DELAY", "LOWPASS"};
    char *comment[] = {"Trigger count",
                       "Output count",
                       "Number of pulses in Channel 0",
                       "Number of pulses in Channel 1",
                       "Trigger level (V)",
                       "Channel 0 threshold (V)",
                       "Channel 1 threshhold (V)",
                       "Range (microseconds)",
                       "Pulse 2 delay (microseconds)",
                       "Lower cut-off frequency (Hz)"};

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
                strcpy(outfilename, "!");
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
    infile = fopen(infilename, "r");
    if (infile == NULL){
        printf("Error opening file %s for reading.\n", infilename);
        exit(1);
    }

    /* Create the output file */
    printf("Creating output file %s.\n\n", outfilename + 1);
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

    /* memory allocation */
    longvals = (long *) malloc(num_longs * sizeof(long));
    doublevals = (double *) malloc(num_doubles * sizeof(double));

    /* Read header line into string */
    header = (char *) malloc(MAX_LINE_LEN * sizeof(char));
    htmp   = (char *) malloc(MAX_LINE_LEN * sizeof(char));
    fgets(header, MAX_LINE_LEN, infile);

    /* Get header values that are longs */
    strdata = strtok(header, "=");
    strdata = strtok(NULL, "\t");
    longvals[0] = strtol(strdata, NULL, 10);
    for (i=1; i < num_longs; i++){
        if (i == 2) strdata = strtok(NULL, "="); /* skip over elapsed time */
        strdata = strtok(NULL, "=");
        strdata = strtok(NULL, "\t");
        longvals[i] = strtol(strdata, NULL, 10);
    }

    /* Get header values that are doubles */
    /* Sometimes the first value is an '=' character. this causes seg faults,
     * but is easily fixed */
    fgets(htmp, MAX_LINE_LEN, infile);
    if (htmp[0] == '=')
        header[0] = ')';
    else
        header[0] = '\0';
    strcat(header, htmp);
    free(htmp);
    strdata = strtok(header, "=");
    strdata = strtok(NULL, "=");
    for (i=0; i<num_doubles; i++){
        strdata = strtok(NULL, "=");
        strdata = strtok(NULL, "\t");
        doublevals[i] = strtod(strdata, NULL);
    }

    /* Print header values */
    printf("Information on this data set:\n");
    for(i=0; i<num_longs; i++)
        printf("%s: %ld\n", comment[i], longvals[i]);
    for (i=0; i<num_doubles; i++)
        printf("%s: %g\n", comment[num_longs+i], doublevals[i]);
    printf("\n");
    
    /* Get column names and units */
    fgets(header, MAX_LINE_LEN, infile);
    fgets(header, MAX_LINE_LEN, infile);

    /* create arrays of column names and units */
    make_table(outfile, header);

    /* fill table */
    printf("Writing data to table.\n");
    nrows = 0;
    while (fgets(header, MAX_LINE_LEN, infile) != NULL){
        nrows++;
        fill_table(outfile, header, nrows);
        printf("Wrote %ld rows to table.\r", nrows);
    }
    printf("Wrote %ld rows to table.\n", nrows);

    /* Write to fits header */
    printf("Writing keys/values to header.\n");
    for (i=0; i<num_longs; i++){
        fits_write_key(outfile, TLONG, labell[i], longvals+i,
                comment[i], &status);
    }
    for (i=0; i<num_doubles; i++){
        fits_write_key(outfile, TDOUBLE, labeld[i], doublevals+i,
                comment[num_longs+i], &status);
    }
    if (status){
        printf("Unable to write to fits header.\n");
        exit(status);
    }

    /* Time to exit the program */
    printf("Closing files.\n");
    fclose(infile);
    fits_close_file(outfile, &status);
    if (status){
        printf("Error closing file %s.\n", outfilename + 1);
        exit(status);
    }
    free(header);
    free(infilename);
    free(outfilename);
    free(longvals);
    free(doublevals);
    exit_prog("mu2fits");
    return 0;
}
