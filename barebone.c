#include <stdio.h>
#include <math.h>

#define SIZE 200

//the code can be modified for time steps 250 and 1000 where the former will show a forward travelling wave and the latter will show a reflected wave due to the PMC boundary condition.

int main()
{
    FILE *barebones;
    // Open a file in writing mode
    barebones = fopen("barebone1.csv", "w");
    if (barebones == NULL) {
        perror("Error opening file");
        return 1;
    }

    double ez[SIZE] = {0.}, hy[SIZE] = {0.}, imp0 = 377.0; //terminates the end node of the magnetic field with a 0 field.
    int qTime, maxTime = 250, mm;

    /* do time stepping */
    for (qTime = 0; qTime < maxTime; qTime++) {
        /* update magnetic field */
        for (mm = 0; mm < SIZE - 1; mm++) {
            hy[mm] = hy[mm] + (ez[mm + 1] - ez[mm]) / imp0;
        }

        /* update electric field */
        for (mm = 1; mm < SIZE; mm++) {
            ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0;
        }

        /* hardwire a source node */
        ez[0] = exp(-(qTime - 30.) * (qTime - 30.) / 100.);

        /* write ez values as comma-separated row */
        fprintf(barebones, "%.6f", ez[50]);
        fprintf(barebones, "\n");  // new line after each time step
    }

    // Close the file
    fclose(barebones);
    return 0;
}
