//this is a copy of the absorbing boundary conditions 1d fdtd.
//designing the coupler
//two lines - waveguide and ring
//define four arrays - electric fields before and after passing through the port for both waveguide and ring.
//ez_r_i ring section before the signal encounters the port
//ez_r_f ring section after it has encountered the port.
//same goes for the waveguide sections as well.

#include <stdio.h>
#include <math.h>

#define SIZE 200
#define SIZE_i 51

//applying an additive source at node ez[50] which will propagate in both directions from that node.

int main()
{
    double p[SIZE] = {0.}; //ports
    double ez_g_i[SIZE_i] = {0.}, hy_g_i[SIZE_i] = {0.}, ez_g_f[SIZE] = {0.}, hy_g_f[SIZE] = {0.}; //i and f stand for initial and final
    double ez_r_i[SIZE_i] = {0.}, hy_r_i[SIZE_i] = {0.}, ez_r_f[SIZE] = {0.}, hy_r_f[SIZE] = {0.}, imp0 = 377.0; //divided by 2 in order to define the length of the 4 branches, the signal is injected at 50.
    int qTime, maxTime = 200, mm;
    int t = 0.5; //arbitrary (will change later)
    int k = 0.3;

    char filename[100] = "port_sim.dat";
    FILE *snapshots;

    snapshots = fopen(filename, "w");

    /* do time stepping */
    for (qTime = 0; qTime < maxTime; qTime++) {

       //hy_r_i[SIZE_i - 1] = hy_r_i[SIZE_i-2];  //boundary conditions to be implemented later to avoid confusion.
        /* update magnetic field */
        for (mm = 0; mm < SIZE_i - 1; mm++) {
            hy_g_i[mm] = hy_g_i[mm] + (ez_g_i[mm + 1] - ez_g_i[mm]) / imp0;
        }

        //ez_r_i[0] = ez_r_i[1];   //implementing boundary conditions [later]

        /* update electric field */
        for (mm = 1; mm < SIZE_i-1; mm++) 
            ez_g_i[mm] = ez_g_i[mm] + (hy_g_i[mm] - hy_g_i[mm - 1]) * imp0;

        
      //instead of the hardwire source we add the additive source

        ez_g_i[50] += exp(-(qTime - 30.) * (qTime - 30.) / 100.); //injecting a Gaussian pulse to the waveguide??

        //defining the 4 ports
        p[0] = ez_g_i[50];
        p[1] = t * p[0]; //through
        p[3] = k * p[0]; //drop

        ez_g_f[0] = p[1]; //check!

        for (mm = 0; mm < SIZE - 1; mm++) {  //will still count from 0 - 199 instead of 50 - 199
            hy_g_f[mm] = hy_g_f[mm] + (ez_g_f[mm + 1] - ez_g_f[mm]) / imp0;
        }

        for (mm = 1; mm < SIZE  -1; mm++) 
            ez_g_f[mm] = ez_g_f[mm] + (hy_g_f[mm] - hy_g_f[mm - 1]) * imp0;

        //through port definition finished ___________

        ez_r_f[0] = p[3]; //do i need to define the first half of the ring branch as well? which is before the port.

        for (mm = 0; mm < SIZE - 1; mm++) {  //will still count from 0 - 199 instead of 50 - 199
            hy_r_f[mm] = hy_r_f[mm] + (ez_r_f[mm + 1] - ez_r_f[mm]) / imp0;
        }

        for (mm = 1; mm < SIZE - 1; mm++) 
            ez_r_f[mm] = ez_r_f[mm] + (hy_r_f[mm] - hy_r_f[mm - 1]) * imp0;

        
        

        //write snapshot if time is a multiple of 2
        if (qTime % 2 == 0) {
            for (mm = 0; mm < SIZE; mm++)
                fprintf(snapshots, "%g ", ez_g_i[mm]);
            fprintf(snapshots, "\n");
  
        }
    }// end of time-stepping
    fclose(snapshots);
    return 0;
}
