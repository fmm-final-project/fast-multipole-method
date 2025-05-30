#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define G 1.0                    // Gravitational constant
#define EPSILON 1.0e-9           // Softening factor to prevent division by zero

typedef struct {
    double mass;
    double pos[3];   // x, y, z
    double vel[3];   // vx, vy, vz
} Particle;

int N;           // Number of particles
unsigned long long count = 0;

void compute_gravity(Particle particles[], double forces[][3]){
    for(int i = 0; i < N; i++){
        forces[i][0] = 0.0;
        forces[i][1] = 0.0;
        forces[i][2] = 0.0;

        for(int j = 0; j < N; j++){
            if(i == j) continue;
            count++;

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

    printf("Start datafile input\n");
    FILE *fptr;
    fptr = fopen("particles.bin", "rb");
    if(!fptr){
        printf("Error opening file!\n");
        exit(1);
    }

    fseek(fptr, 0, SEEK_END);
    int file_size = ftell(fptr);
    rewind(fptr);

    int num_of_doubles = file_size / sizeof(double);
    if(num_of_doubles % 7 != 0){
        printf("File is not a 7 doubles for each particles!\n");
        fclose(fptr);
        exit(1);
    }

    N = num_of_doubles / 7;

    double* buffer = (double*)malloc(file_size);
    if(!buffer){
        printf("Memory allocate failed\n");
        fclose(fptr);
        exit(1);
    }

    if(fread(buffer, sizeof(double), num_of_doubles, fptr)){};

    Particle particles[N];
    for(int i = 0; i < N; i++){
        particles[i].mass = buffer[i * 7];
        particles[i].pos[0] = buffer[i * 7 + 1];
        particles[i].pos[1] = buffer[i * 7 + 2];
        particles[i].pos[2] = buffer[i * 7 + 3];
        particles[i].vel[0] = buffer[i * 7 + 4];
        particles[i].vel[1] = buffer[i * 7 + 5];
        particles[i].vel[2] = buffer[i * 7 + 6];
    }
    fclose(fptr);
    printf("Finish datafile input\n\n");

    double forces[N][3];

    double start, end;
    printf("Start force evaluation\n");
    start = omp_get_wtime();
    compute_gravity(particles, forces);
    end = omp_get_wtime();
    printf("Finish force evaluation\n\n");
    printf("Evaluation count: %llu\n", count);
    double total_time = (end - start) * 1000;
    printf("Total time: %.3lf ms\n", total_time);

    // Output the force vectors
    FILE* fcsv;
    fcsv = fopen("force.csv", "w");
    if(!fcsv){
        printf("Failed to open output file!\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(fcsv, "%.10e,%.10e,%.10e\n", forces[i][0], forces[i][1], forces[i][2]);
    }
    fclose(fcsv);

    FILE* fbin;
    fbin = fopen("force.bin", "wb");
    if(!fbin){
        printf("Failed to open output file!\n");
        exit(1);
    }
    size_t written = fwrite(forces, sizeof(double), 3 * N, fbin);
    if(written != 3 * N){
        printf("Failed to write all data!\n");
    }
    fclose(fbin);
    printf("Finish output\n");

    return 0;
}
