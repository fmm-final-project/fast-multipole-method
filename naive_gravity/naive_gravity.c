#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define EPSILON 1.0e-12          // Softening factor to prevent division by zero
#define MAX_LINE_LEN 1024

typedef struct {
    double mass;
    double pos[3];   // x, y, z
    double vel[3];   // vx, vy, vz
} Particle;

// Input parameters
double G;
double THETA;
int MAX_PARTICLES_PER_CELL;
char datafile[MAX_LINE_LEN];
char outfile[MAX_LINE_LEN];

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

    FILE *fp = fopen("tree.in", "r");
    if(!fp){
        printf("Failed to open file\n");
        exit(1);
    }
    char line[MAX_LINE_LEN];
    int lcount = 0;
    while(fgets(line, sizeof(line), fp)){
        if(line[0] == '#' || line[0] == '\n'){
            continue;
        }
        lcount++;
        if(lcount == 1){
            if(sscanf(line, "%lf", &G)){}
        }
        else if(lcount == 2){
            if(sscanf(line, "%lf", &THETA)){}
        }
        else if(lcount == 3){
            if(sscanf(line, "%d", &MAX_PARTICLES_PER_CELL)){}
        }
        else if(lcount == 4){
            if(sscanf(line, "%s", datafile)){}
        }
        else if(lcount == 5){
            if(sscanf(line, "%s", outfile)){}
        }
    }

    fclose(fp);

    printf("Start datafile input\n");
    FILE *fptr;
    fptr = fopen(datafile, "rb");
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
    Particle* particles = (Particle*)malloc(sizeof(Particle) * N);
    if (!particles) {
        printf("Memory allocate failed\n");
        free(buffer);
        fclose(fptr);
        exit(1);
    }

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

    double (*forces)[3] = malloc(sizeof(double) * N * 3);
    if (!forces) {
        printf("Memory allocate failed\n");
        free(buffer);
        free(particles);
        fclose(fptr);
        exit(1);
    }

    double start, end;
    printf("Start force evaluation\n");
    start = omp_get_wtime();
    compute_gravity(particles, forces);
    end = omp_get_wtime();
    printf("Finish force evaluation\n\n");
    printf("Evaluation count: %llu\n", count);
    double total_time = (end - start) * 1000;
    printf("Total execution time: %.5lf ms\n", total_time);

    // Output the force vectors
    printf("Start output\n");
    /*
    FILE* fcsv;
    fcsv = fopen(outfile, "w");
    if(!fcsv){
        printf("Failed to open output file!\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(fcsv, "%.10e,%.10e,%.10e\n", forces[i][0], forces[i][1], forces[i][2]);
    }
    fclose(fcsv);
    */

    FILE* fbin = fopen(outfile, "wb");
    if(!fbin){
        perror("Failed to open file");
        return 1;
    }

    // Write the array to file
    size_t written = fwrite(forces[0], sizeof(double), 3 * N, fbin);
    if(written != 3 * N){
        perror("Failed to write data");
    }
    // Close the file
    fclose(fbin);
    printf("Finish output\n");

    free(forces);
    free(buffer);
    free(particles);

    return 0;
}