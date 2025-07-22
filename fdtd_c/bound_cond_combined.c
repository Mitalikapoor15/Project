#include <stdio.h>
#include <math.h>

#define SIZE 200

//applying an additive source at node ez[50] which will propagate in both directions from that node.

int main()
{
    double ez[SIZE] = {0.}, hy[SIZE] = {0.}, imp0 = 377.0; //terminates the end node of the magnetic field with a 0 field.
    int qTime, maxTime = 200, mm;

    char filename[100] = "sim.dat";
    FILE *snapshots;

    snapshots = fopen(filename, "w");

    /* do time stepping */
    for (qTime = 0; qTime < maxTime; qTime++) {
        /* update magnetic field */
        for (mm = 0; mm < SIZE - 1; mm++) {
            hy[mm] = hy[mm] + (ez[mm + 1] - ez[mm]) / imp0;
        }

        /* update electric field */
        for (mm = 1; mm < SIZE; mm++) 
            ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0;
        
      //instead of the hardwire source we add the additive source

        ez[50] += exp(-(qTime - 30.) * (qTime - 30.) / 100.);

        //write snapshot if time is a multiple of 2
        if (qTime % 2 == 0) {
            for (mm = 0; mm < SIZE; mm++)
                fprintf(snapshots, "%g ", ez[mm]);
            fprintf(snapshots, "\n");
  
        }
    }// end of time-stepping
    fclose(snapshots);
    return 0;
}
