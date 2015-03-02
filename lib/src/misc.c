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
#include <time.h>

/*******************************************************************************
 * init_prog
 * Initialize a program.
 * Input:
 *      progname        The name of the program being run
 ******************************************************************************/
void init_prog(char *progname){
    time_t     current;
    struct tm *date;

    time(&current);
    date = localtime(&current);
    printf("Initializing program %s on %s", progname, asctime(date));
    printf("Author: Rachel Domagalski\n\n");
}

/*******************************************************************************
 * init_prog
 * Initialize a program.
 * Input:
 *      progname        The name of the program being run
 ******************************************************************************/
void exit_prog(char *progname){
    time_t     current;
    struct tm *date;

    time(&current);
    date = localtime(&current);
    printf("Exiting program %s on %s", progname, asctime(date));
}
