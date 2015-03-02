/*******************************************************************************
 * This code is for problem 3 of the EAX homework for Physics 111 Advanced Lab
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
#include <time.h>
#include <fitsio.h>

/* For random number generation. */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Global random number generator */
gsl_rng *randnumgen;

/*******************************************************************************
 * double *sample_gauss
 * Get the mean, standard deviation, and error of N samples from a gaussian.
 * Input:
 *      num_samples:    Number of samples.
 *      true_mean:      Actual mean of the distribution
 *      true_sigma:     Actual standard deviation
 * Output:
 *      result:         Array containing computed mean, sigma, and error
 ******************************************************************************/
double *sample_gauss(long num_samples, double true_mean, double true_sigma){
    double *result, *samples, mean = 0, sigma = 0, error;
    long i;

    /* Exit if there is only one sample; */
    if (num_samples == 1) return NULL;
    
    /* Memory allocation */
    result = (double *) malloc(3*sizeof(double));
    samples = (double *) malloc(num_samples*sizeof(double));

    for (i=0; i<num_samples; i++){
        samples[i] = true_mean + gsl_ran_gaussian(randnumgen, true_sigma);
        mean += samples[i];
    }
    mean /= num_samples;

    for (i=0; i<num_samples; i++){
        sigma += pow(mean - samples[i], 2);
    }
    free(samples); /* This memory is no longer needed. */
    sigma /= num_samples - 1;
    sigma = sqrt(sigma);
    error = sigma / sqrt(num_samples);
    
    result[0] = mean;
    result[1] = sigma;
    result[2] = error;
    return result;
}

/*******************************************************************************
 * gaussian_test
 * This function creates a fits table, then populates it with the mean, sigma,
 * and error on the mean from experiments sampling from a gaussian num_samples
 * number of times.
 * Input:
 *      fptr:           Fits file to write the table to
 *      num_iter:       Number of experiments to run
 *      num_samples:    Number of samples for each experiment
 *      mean:           Mean for the distibution in the experiments
 *      sigma:          Sigma for the distribution in the experiments
 * Output:
 *      status:         Fits error status
 ******************************************************************************/
void gaussian_test(fitsfile *fptr, long num_iter, long num_samples, double mean,
        double sigma, int *status){
    long nrows, i, num_sigma;
    double *stats, ratio_sigma, ave_mean, mean_sd, *all_means;
    /* Fits file columns */
    char *ttype[] = {"mean", "sigma", "error"};
    char *tform[] = {"D", "D", "D"};
    char *tunit[] = {NULL, NULL, NULL};
    char *key[] = {"NSAMPLES", "MEAN", "SIGMA", "MEANAVE", "MEANSD"};
    char *comment[] = {"Number of samples per iteration",
        "True mean of the gaussian",
        "True standard deviation of the gaussian",
        "Average mean",
        "Standard deviation of all means"};

    /* Create the fits table */
    fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, tunit, NULL, status);
    if (*status){
        printf("Error creating fits table in guassian_test.\n");
        exit(*status);
    }

    /* Get mean, sigma, and error, then save to a fits file */
    all_means = (double *) malloc(num_iter * sizeof(double));
    ave_mean = 0;
    for (nrows = 0; nrows < num_iter; nrows++){
        stats = sample_gauss(num_samples, mean, sigma);
        ave_mean += stats[0];
        all_means[nrows] = stats[0];
        /* Write to file */
        for (i=0; i < 3; i++){
            fits_write_col(fptr, TDOUBLE, i+1, nrows+1, 1L,1L, stats+i, status);
            if (*status){
                printf("Error writing to fits table in gaussian_test.\n");
                exit(1);
            }
        }
        free(stats); /* free memory for next iteration */
    }
    ave_mean /= num_iter;

    /* Get standard deviation of all means */
    mean_sd = 0;
    for (i=0; i<num_iter; i++)
        mean_sd += pow(ave_mean - all_means[i], 2);
    mean_sd /= num_iter - 1;
    mean_sd = sqrt(mean_sd);

    /* Get percentage of experiments within 2 sigma */
    num_sigma = 0;
    for (i=0; i<num_iter; i++){
        if (fabs(all_means[i] - ave_mean) < 2 * mean_sd)
            num_sigma++;
    }
    ratio_sigma = (double) num_sigma / num_iter;
    free(all_means);

    /* print some results */
    printf("Average mean: %g\n", ave_mean);
    printf("Sigma of all means: %g\n", mean_sd);
    printf("Percentage of experiments within 2 sigma: ");
    printf("%g%c\n", 100 * ratio_sigma, (char) 37);

    /* Write the number of samples to the fits table */
    fits_write_key(fptr, TLONG, key[0], &num_samples, comment[0], status);
    fits_write_key(fptr, TDOUBLE, key[1], &mean, comment[1], status);
    fits_write_key(fptr, TDOUBLE, key[2], &num_sigma, comment[2], status);
    fits_write_key(fptr, TDOUBLE, key[3], &ave_mean, comment[3], status);
    fits_write_key(fptr, TDOUBLE, key[4], &mean_sd, comment[4], status);
    if (*status){
        printf("Error writing to fits header.\n");
        exit(*status);
    }
}

/*******************************************************************************
 * gauss_stats
 * This function ggenerates sampling experiments with gaussians and collects the
 * statistics on them.
 ******************************************************************************/
int main(int argc, char **argv){
    double *stats;
    double mean = 0, sigma = 1;
    int i, status = 0;
    long num_iter;
    long num_samples[5] = {5, 10, 50, 100, 1000};
    /* Variables for reading in options */
    int opt;
    fitsfile *outfile;
    char *filename;
    
    /***************************************************************************
     * Initialize random number generator using time for the seed. The RNG
     * algorithm is Mersenne Twister.
     **************************************************************************/
    randnumgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(randnumgen, time(NULL));

    /* Parse options */
    filename = NULL;
    while ((opt = getopt(argc, argv, "o:")) != -1){
        switch (opt){
            case 'o':
                filename = (char *) malloc((strlen(optarg)+2)*sizeof(char));
                /* Prepend ! to filename to force overwriting it. */
                filename[0] = '!';
                strcat(filename, optarg);
                printf("Writing output to %s.\n\n", filename + 1);
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
    if (filename == NULL){
        printf("Error: No output file specified.\n");
        exit(1);
    }

   /* Create file and generate results for it. */
   fits_create_file(&outfile, filename, &status);
   if (status){
       printf("Error creating file %s.\n", filename + 1);
       exit(status);
   }
   fits_open_file(&outfile, filename+1, READWRITE, &status);
   if (status){
       printf("Error opening file %s.\n", filename + 1);
       exit(status);
   }
    
    /* Loop over number of samples */
    num_iter = 1000;
    for (i = 0; i < 5; i++){
        stats = sample_gauss(num_samples[i], mean, sigma);
    
        /* Print example list. */
        printf("Example list for %ld number of samples.\n", num_samples[i]);
        printf("Mean: %g\n", stats[0]);
        printf("Standard deviation: %g\n", stats[1]);
        printf("Error on the mean: %g\n", stats[2]);
        printf("Expected error on the mean: %g\n", 1 / sqrt(num_samples[i]));
        printf("\n");
        free(stats);
    
        /* Sample guassian and save to file */
        printf("Sampling from %ld gaussians %ld times each.\n",
                num_iter, num_samples[i]);
        gaussian_test(outfile, num_iter, num_samples[i], mean, sigma, &status);
        if (status){
            printf("Error: unable to write to file.\n");
            exit(status);
        }
        printf("\nWritten to file %s.\n\n", filename + 1);
    }

    printf("Closing file %s.\n\n", filename + 1);
    fits_close_file(outfile, &status);
    if (status){
        printf("Error closing file %s.\n", filename + 1);
        exit(status);
    }
    free(stats);
    free(filename);
    return 0;
}
