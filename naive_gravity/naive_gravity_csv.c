#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_LINE 1024
#define G 1.0                    // Gravitational constant
#define EPSILON 1.0e-9           // Softening factor to prevent division by zero

typedef struct {
    double mass;
    double pos[3];   // x, y, z
    double vel[3];   // vx, vy, vz (not used here)
} Particle;

int N = 100;           // Number of particles

void compute_gravity(Particle particles[], double forces[][3]){
    for(int i = 0; i < N; i++){
        forces[i][0] = 0.0;
        forces[i][1] = 0.0;
        forces[i][2] = 0.0;

        for(int j = 0; j < N; j++){
            if(i == j) continue;

            double dx = particles[j].pos[0] - particles[i].pos[0];
            double dy = particles[j].pos[1] - particles[i].pos[1];
            double dz = particles[j].pos[2] - particles[i].pos[2];

            double dist_sq = dx*dx + dy*dy + dz*dz + EPSILON;
            double dist = sqrt(dist_sq);
            double force_mag = G * particles[i].mass * particles[j].mass / dist_sq;

            // Directional force components
            forces[i][0] += force_mag * dx / dist;
            forces[i][1] += force_mag * dy / dist;
            forces[i][2] += force_mag * dz / dist;
        }
    }
}

int main() {
    Particle particles[N];

    FILE *fptr;
    fptr = fopen("particles.csv", "r");
    if(!fptr){
        perror("Error opening file!");
        exit(1);
    }

    int count = 0;
    char line[MAX_LINE];
    while(fgets(line, sizeof(line), fptr)){
        line[strcspn(line, "\n")] = '\0';

        char* token = strtok(line, ",");
        particles[count].mass = atof(token);
        for(int i = 0; i < 3; i++){
            token = strtok(NULL, ",");
            particles[count].pos[i] = atof(token);
        }
        for(int i = 0; i < 3; i++){
            token = strtok(NULL, ",");
            particles[count].vel[i] = atof(token);
        }
        count++;
    }
    fclose(fptr);

    double forces[N][3];

    compute_gravity(particles, forces);

    // Output the force vectors
    for (int i = 0; i < N; i++) {
        printf("Particle %d: Fx = %.5e, Fy = %.5e, Fz = %.5e\n",
               i, forces[i][0], forces[i][1], forces[i][2]);
    }

    return 0;
}
