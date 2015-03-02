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

/*******************************************************************************
 * make_table
 * This function creates a table using the columns from an input file, then
 * copies data into the output.
 * Input:
 *      infptr:         File to get the columns from.
 *      outfile:        File to copy the columns to.
 *      names:          Filenames of the input files (for error handling).
 *      num_files       Number of input files.
 * Output:
 *      status:         Fits error status
 ******************************************************************************/
int make_table(fitsfile **infptr, fitsfile *outfptr, 
        char **names, int num_files){
    double datapt = 0;
    int i, k, anynul = 0, ncols, status = 0;
    long j, nrows, total_rows;

    /* Create table and copy the columns */
    fits_create_tbl(outfptr, BINARY_TBL, 0, 0, NULL, NULL, NULL, NULL, &status);
    if (status){
        printf("Error creating FITS table.\n");
        close_fits_files(infptr, names, num_files, 0);
        exit(status);
    }
    for (i=0; i<NUM_COLS; i++){
        fits_copy_col(*infptr, outfptr, i+1, i+1, 1, &status);
        if (status){
            printf("Unable to copy FITS columns.\n");
            close_fits_files(infptr, names, num_files, 0);
            exit(status);
        }
    }

    /* Copy data to output file */
    ncols      = 0;
    nrows      = 0;
    total_rows = 0;
    for (i=0; i<num_files; i++){
        printf("Copying data from %s.\n", names[i]);
        /* Get number of rows to copy */
        fits_get_num_cols(infptr[i], &ncols, &status);
        fits_get_num_rows(infptr[i], &nrows, &status);
        if (status){
            printf("Error retriving info about FITS table.\n");
            close_fits_files(infptr, names, num_files, 0);
            exit(status);
        }

        /* Copy data */
        for (j=0; j<nrows; j++){
            for (k=0; k<ncols; k++){
                fits_read_col(infptr[i], TDOUBLE, k+1, j+1, 1L, 1L, NULL,
                        &datapt, &anynul, &status);
                fits_write_col(outfptr, TDOUBLE, k+1, total_rows+j+1, 1L, 1L,
                        &datapt, &status);
                if (status){
                    printf("Unable to copy data.\n");
                    close_fits_files(infptr, names, num_files, 0);
                    exit(status);
                }
            }
        }
        total_rows += nrows;
    }
    printf("\n");

    return status;
}

/*******************************************************************************
 * trigger_output_count
 * Gets the total trigger and output counts for all files.
 * Input:
 *      fptr:           Array of fitsfile pointers
 *      names:          Filenames (for error handling)
 *      num_files:      Number of fits files
 * Output:
 *      trigger_count:  Total trigger count
 *      output_count:   Total output count
 *      status:         Fits error status
 ******************************************************************************/
int trigger_output_count(fitsfile **fptr, char **names, int num_files,
        long *trigger_count, long *output_count){
    char *comment;
    int i, status = 0;
    long trigger = 0, output = 0;

    *trigger_count = 0;
    *output_count  = 0;
    comment = (char *) malloc(81 * sizeof(char));
    for (i=0; i<num_files; i++){
        fits_read_key(fptr[i], TLONG, "TRCOUNT", &trigger, comment, &status);
        fits_read_key(fptr[i], TLONG, "OCOUNT",  &output,  comment, &status);
        if (status){
            printf("Unable to read from FITS headers.\n");
            close_fits_files(fptr, names, num_files, 0);
            exit(status);
        }
        *trigger_count += trigger;
        *output_count  += output;
    }

    free(comment);
    return status;
}

/*******************************************************************************
 * catmufits
 * Combine two or more fits files into a single file. The fits file that will
 * fill the output file first determine the headers to use to determine whether
 * or not the rest of the files are compatible and can be combined.
 ******************************************************************************/
int main(int argc, char **argv){
    init_prog("catmufits");
    char **infilenames, *outfilename, *comment;
    int opt, status = 0;
    fitsfile **infiles, *outfile;
    /* For writing fits headers */
    int i, j, num_longs = 4, num_args;
    double checkdblkey, *doublevals;
    long arglen, checklngkey, *longvals, trigger_count = 0, output_count = 0;
    long ndbl = 6, nlng = 2;
    char *labell[] = {"NCH0", "NCH1"};
    char *labeld[] = {"TRLEVEL", "CHTH0", "CHTH1", "RANGE", "DELAY", "LOWPASS"};
    char *allcomments[] = {"Trigger count",
                           "Output count",
                           "Number of pulses in Channel 0",
                           "Number of pulses in Channel 1",
                           "Trigger level (V)",
                           "Channel 0 threshold (V)",
                           "Channel 1 threshhold (V)",
                           "Range (microseconds)",
                           "Pulse 2 delay (microseconds)",
                           "Lower cut-off frequency (Hz)"};

    /* Storage for fits headers being read */
    comment    = (char *)   malloc(81   * sizeof(char));
    longvals   = (long *)   malloc(nlng * sizeof(long));
    doublevals = (double *) malloc(ndbl * sizeof(double));

    /* Parse options */
    outfilename = NULL;
    while ((opt = getopt(argc, argv, "o:")) != -1){
        switch (opt){
            case 'o':
                outfilename = (char *) malloc((strlen(optarg)+2)*sizeof(char));
                /* Prepend ! to filename to force overwriting it. */
                strcpy(outfilename, "!");
                strcat(outfilename, optarg);
                break;
            case '?':
                if (optopt == 'o'){
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
    if (outfilename == NULL){
        printf("Error: No output file specified.\n");
        exit(1);
    }

    /* Determine how many files to concatenate */
    num_args = argc - optind;
    if (!num_args){
        printf("Error: No input files specified.\n");
        exit(1);
    }
    infiles = (fitsfile **) malloc(num_args * sizeof(fitsfile *));
    infilenames = (char **) malloc(num_args * sizeof(char *));

    /* Open files */
    for (i=0; i<num_args; i++){
        arglen = 1 + strlen(argv[optind+i]);
        infilenames[i] = (char *) malloc(arglen * sizeof(char));
        strcpy(infilenames[i], argv[optind+i]);
        
        printf("Opening file %s.\n", infilenames[i]);

        fits_open_table(infiles+i, infilenames[i], READONLY, &status);
        if (status){
            printf("Error opening %s.\n", infilenames[i]);
            close_fits_files(infiles, infilenames, i, 0);
            exit(status);
        }
    }

    /* Get fits headers/keys from first infile */
    printf("\nGetting keys to match from %s.\n\n", *infilenames);
    for (i=0; i<nlng; i++){
        fits_read_key(*infiles, TLONG, labell[i], longvals+i, comment, &status);
        if (status){
            printf("Error reading keys from %s.\n", *infilenames);
            close_fits_files(infiles, infilenames, num_args, 0);
            exit(status);
        }
    }
    for (i=0; i<ndbl; i++){
        fits_read_key(*infiles, TDOUBLE, labeld[i],
                doublevals+i, comment, &status);
        if (status){
            printf("Error reading keys from %s.\n", *infilenames);
            close_fits_files(infiles, infilenames, num_args, 0);
            exit(status);
        }
    }

    /* Check the rest of the header keys to see if they match the first file */
    printf("Checking that all files are compatible.\n");
    for (i=1; i<num_args; i++){
        for (j=0; j<nlng; j++){
            fits_read_key(infiles[i], TLONG, labell[j],
                    &checklngkey, comment, &status);
            if (status){
                printf("Error reading keys from %s.\n", infilenames[i]);
                close_fits_files(infiles, infilenames, num_args, 0);
                exit(status);
            }

            if (checklngkey != longvals[j]){
                printf("Header key for %s is invalid.\n", labell[j]);
                close_fits_files(infiles, infilenames, num_args, 0);
                exit(1);
            }
        }
        for (j=0; j<ndbl; j++){
            fits_read_key(infiles[i], TDOUBLE, labeld[j],
                    &checkdblkey, comment, &status);
            if (status){
                printf("Error reading keys from %s.\n", infilenames[i]);
                close_fits_files(infiles, infilenames, num_args, 0);
                exit(status);
            }
            
            if (checkdblkey != doublevals[j]){
                printf("Header key for %s is invalid.\n", labell[j]);
                close_fits_files(infiles, infilenames, num_args, 0);
                exit(1);
            }
        }
    }

    /* Create the output file */
    printf("Creating output file %s.\n\n", outfilename + 1);
    fits_create_file(&outfile, outfilename, &status);
    if (status){
        printf("Error creating file %s.\n", outfilename + 1);
        close_fits_files(infiles, infilenames, num_args, 0);
        exit(status);
    }
    fits_open_file(&outfile, outfilename + 1, READWRITE, &status);
    if (status){
        printf("Error opening file %s.\n", outfilename + 1);
        close_fits_files(infiles, infilenames, num_args, 0);
        exit(status);
    }

    /* Create table in the output file and copy data */
    make_table(infiles, outfile, infilenames, num_args);

    /* Write to fits header */
    printf("Copying FITS header.\n\n");
    trigger_output_count(infiles, infilenames, num_args,
            &trigger_count, &output_count);
    fits_write_key(outfile, TLONG, "TRCOUNT", &trigger_count,
            allcomments[i], &status);
    fits_write_key(outfile, TLONG, "OCOUNT", &output_count,
            allcomments[i], &status);
    for (i=0; i<nlng; i++)
        fits_write_key(outfile, TLONG, labell[i], longvals + i,
                allcomments[i+2], &status);
    for (i=0; i<ndbl; i++)
        fits_write_key(outfile, TDOUBLE, labeld[i], doublevals + i,
                allcomments[num_longs+i], &status);
    if (status){
        printf("Unable to write to FITS header.\n");
        close_fits_files(infiles, infilenames, num_args, 0);
        exit(status);
    }

    /* Time to exit the program */
    close_fits_files(infiles, infilenames, num_args, 1);
    printf("Closing file %s.\n\n", outfilename + 1);
    fits_close_file(outfile, &status);
    if (status){
        printf("Error closing file %s.\n", outfilename + 1);
        close_fits_files(infiles, infilenames, num_args, 0);
        exit(status);
    }
    free(outfilename);
    exit_prog("catmufits");
    return 0;
}
