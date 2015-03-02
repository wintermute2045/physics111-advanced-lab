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

#ifndef _MUO_H
#define _MUO_H 1

#define NUM_COLS 12

/* Pulses to cut on */
typedef enum {IGNORE, PULSE1, PULSE2, BOTH} pcut;

/* Structure of cuts and corrections to make */
typedef struct {
    double  vmin;       /* Minimum voltage allowed */
    double  vmax;       /* Maximum voltage allowed */
    pcut    pulse;      /* The pulses to make the cut on */
    double  (*fix_time) (double time, void *params); /* Apply time correction */
    void    *params;    /* Parameters to correction function */
    int     apply_fix;  /* Whether or not to apply a time correction */
} cutmaker;

/* In muo_common.c */

int close_fits_files(fitsfile **fptr, char **names, int num_files, int verbose);
gsl_histogram *fill_hist(fitsfile *fptr,
        long nbins, double max, int col1, int col2, cutmaker *cut);
gsl_histogram *fill_hist_singlet(fitsfile *fptr,
        long nbins, double max, int col);

#endif
