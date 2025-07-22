//this is a copy of the absorbing boundary conditions 1d fdtd.
//designing the coupler
//two lines - waveguide and ring
//define four arrays - electric fields before and after passing through the port for both waveguide and ring.
//ez_r_i ring section before the signal encounters the port
//ez_r_f ring section after it has encountered the port.
//same goes for the waveguide sections as well.
//boundary conditions have been applied.

#include <stdio.h>
#include <math.h>

#define SIZE 100
#define NUM_PORTS 4

//applying an additive source at node ez[50] which will propagate in both directions from that node.

int main()
{
    double ez[NUM_PORTS][SIZE] = {0.}; // ports 0 and 1 for guide; 2 and 3 for ring
    double hy[NUM_PORTS][SIZE] = {0.};
    // double ez_g_i[SIZE_i] = {0.}, hy_g_i[SIZE_i] = {0.}, ez_g_f[SIZE] = {0.}, hy_g_f[SIZE] = {0.}; //i and f stand for initial and final
    // double ez_r_i[SIZE_i] = {0.}, hy_r_i[SIZE_i] = {0.}, ez_r_f[SIZE] = {0.}, hy_r_f[SIZE] = {0.}, 
    double p[SIZE] = {0.}; //ports
    double imp0 = 377.0; //divided by 2 in order to define the length of the 4 branches, the signal is injected at 50.
    int qTime, maxTime = 200, j;
    int t = 3; //arbitrary (will change later)
    int k = 4;
    int connect = 50;  //grid value where the port is made

    const int num_ports = 4;
    char *filenames[] = {
        "portData\\port1.dat",
        "portData\\port2.dat",
        "portData\\port3.dat",
        "portData\\port4.dat"
    };
    
    FILE *snapshots[4];

    // Open all files
    for (int i = 0; i < num_ports; i++) {
        snapshots[i] = fopen(filenames[i], "w");
        if (snapshots[i] == NULL) {
            printf("Error opening file %s\n", filenames[i]);
            return 1;
        }
    }


    char *filename = "portData\\p_vals.dat";  //always use forward slash in addresses. backslash is interpreted as an escape sequence in C.

    FILE *snapshot = fopen(filename, "w");

    for (qTime = 0; qTime < maxTime; qTime++) {

         /* do time stepping */

           //UPDATE equations for the waveguide written separately
           
           hy[0][SIZE - 1] = hy[0][SIZE-2];  //boundary conditions to be implemented later to avoid confusion.

            /* update magnetic field */
            for (j = 0; j < SIZE - 1; j++) {
                hy[0][j] = hy[0][j] + (ez[0][j + 1] - ez[0][j]) / imp0;
            }

            ez[0][0] = ez[0][1];   //implementing boundary conditions [later]

            /* update electric field */
            for (j = 1; j < SIZE - 1; j++) 
                ez[0][j] = ez[0][j] + (hy[0][j] - hy[0][j - 1]) * imp0;


          //adding a source to the waveguide

            ez[0][10] += exp(-(qTime - 10.) * (qTime - 10.) / 100.); //injecting a Gaussian pulse at 10th node of the i/p waveguide

            //defining the 4 ports
            p[0] = ez[0][50];  //port is fixed
            p[1] = t * p[0]; //through
            p[2] = k * p[0]; //drop

            // for (int i = 0; i < 3; i++) {
            // printf("p[%d] = %.6f\n", i, p[i]);
            // }

            for (int i = 0; i < (NUM_PORTS - 1); i++) {
                fprintf(snapshot, "%g ", p[i]);
            }
            fprintf(snapshot, "\n");    
             //stored the data of the three defined ports (so far). Observed that the ports start having the signals around the 40th time-step.
             //The file is called p_vals.dat
            

            for (int i = 1; i < (NUM_PORTS - 1); i++) {  //excludes 3rd port for now
             
                ez[i][0] = p[i];  //transferring signal across ports
               
                hy[0][SIZE - 1] = hy[0][SIZE-2]; //boundary condition

                /* update magnetic field */
                for (j = 0; j < (SIZE - 1); j++) {
                    hy[i][j] = hy[i][j] + (ez[i][j + 1] - ez[i][j]) / imp0;
                }
              //implementing boundary condition
                ez[i][0] = ez[i][1]; 
                /* update electric field */
                for (j = 1; j < (SIZE - 1); j++) 
                    ez[i][j] = ez[i][j] + (hy[i][j] - hy[i][j - 1]) * imp0;
            }
        //write snapshot if time is a multiple of 2
        for (int i = 0; i < num_ports; i++) {
            if (qTime % 2 == 0) {
                for (j = 0; j < SIZE; j++)
                    fprintf(snapshots[i], "%g ", ez[i][j]);
                }
            fprintf(snapshots[i], "\n");
        }
    } // end of time-stepping
    
    for (int i = 0; i < num_ports; i++) {
    fclose(snapshots[i]);
    }
    fclose(snapshot);
return 0;
}
