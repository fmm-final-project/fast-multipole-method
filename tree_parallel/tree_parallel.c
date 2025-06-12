#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define MAX_LINE_LEN 1024

typedef struct{
    double x[3];
    double mass;
    double force[3];
} Particle;

typedef struct Cell{
    double center[3];     // cell geometric center
    double size;          // cell edge length
    double mass;          // total mass (monopole)
    double massCenter[3]; // center of mass
    // dipole vanished by choosing CoM as center
    double quad[3][3];    // quadrupole tensor
    int nParticles;
    Particle** particles;
    struct Cell* children[8];
} Cell;

// Input parameters
double G;
double THETA;
int MAX_PARTICLES_PER_CELL;
char datafile[MAX_LINE_LEN];
char outfile[MAX_LINE_LEN];
int NUM_OF_THREADS;

// Performance
unsigned long long direct_count = 0;
double start, end;

// Create a new cell, compute its center and halfSize from particle extents
Cell* createCell(Particle** particles, int n){
    Cell* cell = (Cell*)malloc(sizeof(Cell));
    cell->nParticles = n;
    cell->particles = (Particle**)malloc(n * sizeof(Particle*));
    for(int i = 0; i < n; i++){
        cell->particles[i] = particles[i];
    }
    // Compute bounding box
    double min[3], max[3];
    for(int d = 0; d < 3; d++){
        min[d] = max[d] = particles[0]->x[d];
    }
    for(int i = 1; i < n; i++){
        for(int d = 0; d < 3; d++){
            double v = particles[i]->x[d];
            if(v < min[d]) min[d] = v;
            if(v > max[d]) max[d] = v;
        }
    }
    // Set center and halfSize
    for(int d = 0; d < 3; d++){
        cell->center[d] = 0.5 * (min[d] + max[d]);
    }
    double dx = max[0] - min[0];
    double dy = max[1] - min[1];
    double dz = max[2] - min[2];
    cell->size = fmax(fmax(dx, dy), dz);
    // Initialize children
    for(int c = 0; c < 8; c++) cell->children[c] = NULL;
    return cell;
}

void subdivideCell(Cell *cell){
    if(cell->nParticles <= MAX_PARTICLES_PER_CELL) return;

    // Distribute particles into octants
    Particle** plist = cell->particles;
    int counts[8] = {0};
    for(int i = 0; i < cell->nParticles; i++){
        Particle* p = plist[i];
        int idx = (p->x[0] > cell->center[0])
                + 2 * (p->x[1] > cell->center[1])
                + 4 * (p->x[2] > cell->center[2]);
        counts[idx]++;
    }
    Particle*** lists = (Particle***)malloc(8 * sizeof(Particle**));
    for(int c = 0; c < 8; c++){
        lists[c] = (Particle**)malloc(counts[c] * sizeof(Particle*));
        counts[c] = 0;
    }
    for(int i = 0; i < cell->nParticles; i++){
        Particle* p = plist[i];
        int idx = (p->x[0] > cell->center[0])
                + 2 * (p->x[1] > cell->center[1])
                + 4 * (p->x[2] > cell->center[2]);
        lists[idx][counts[idx]++] = p;
    }
    free(cell->particles);
    for(int c = 0; c < 8; c++){
        if(counts[c] > 0){
            cell->children[c] = createCell(lists[c], counts[c]);
            subdivideCell(cell->children[c]);
        }
        free(lists[c]);
    }
    free(lists);
    cell->nParticles = 0;
    cell->particles = NULL;
}

void computeMassDistribution(Cell *cell){
    // Initialize multipoles
    cell->mass = 0.0;
    for(int i = 0; i < 3; i++){
        cell->massCenter[i] = 0.0;
        for(int j = 0; j < 3; j++) cell->quad[i][j] = 0.0;
    }
    if(cell->nParticles > 0){
        // Leaf: compute multipoles directly
        for(int i = 0; i < cell->nParticles; i++){
            Particle* p = cell->particles[i];
            cell->mass += p->mass;
            for(int d = 0; d < 3; d++) cell->massCenter[d] += p->mass * p->x[d];
        }
        for(int d = 0; d < 3; d++) cell->massCenter[d] /= cell->mass;
        // Quadrupole
        for(int i = 0; i < cell->nParticles; i++){
            Particle* p = cell->particles[i];
            double dx[3], r2 = 0;
            for (int d = 0; d < 3; d++){
                dx[d] = p->x[d] - cell->massCenter[d];
                r2 += dx[d] * dx[d];
            }
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    cell->quad[a][b] += p->mass * (3*dx[a]*dx[b] - (a==b ? r2 : 0.0));
                }
            }
        }
    }
    else{
        // Internal: aggregate from children
        for(int c = 0; c < 8; c++){
            if (!cell->children[c]) continue;
            Cell* ch = cell->children[c];
            computeMassDistribution(ch);
            cell->mass += ch->mass;
            for(int d = 0; d < 3; d++) cell->massCenter[d] += ch->mass * ch->massCenter[d];
        }
        for(int d = 0; d < 3; d++) cell->massCenter[d] /= cell->mass;
        for(int c = 0; c < 8; c++){
            if(!cell->children[c]) continue;
            Cell* ch = cell->children[c];
            double dC[3], r2 = 0;
            for(int d = 0; d < 3; d++){
                dC[d] = ch->massCenter[d] - cell->massCenter[d];
                r2 += dC[d] * dC[d];
            }
            for(int a = 0; a < 3; a++){
                for(int b = 0; b < 3; b++){
                    cell->quad[a][b] += ch->quad[a][b] + ch->mass * (3*dC[a]*dC[b] - (a==b ? r2 : 0.0));
                }
            }
        }
    }
}

unsigned long long evaluateForceOnParticle(Particle* p, Cell* cell){
    unsigned long long local_count = 0;
    if(cell->nParticles > 0){
        // Direct summation
        for(int i = 0; i < cell->nParticles; i++){
            Particle* q = cell->particles[i];
            if(p == q) continue;
            local_count++;
            double dx[3], r2 = 1e-12;
            for(int d = 0; d <3; d++){
                dx[d] = p->x[d] - q->x[d];
                r2 += dx[d] * dx[d];
            }
            double inv3 = 1.0 / (r2 * sqrt(r2));
            double f = -G * p->mass * q->mass * inv3;
            for(int d = 0; d < 3; d++){
                p->force[d] += f * dx[d];
            }
        }
    }
    else{
        // Multipole approximation
        double dx[3], r2 = 1e-12;
        for(int d = 0; d < 3; d++){
            dx[d] = p->x[d] - cell->massCenter[d];
            r2 += dx[d] * dx[d];
        }
        double r = sqrt(r2);
        if(cell->size < r * THETA){
            double r3 = r2 * r;
            double r5 = r2 * r3;
            double r7 = r2 * r5;
            // Monopole
            double monoF = -G * p->mass * cell->mass / r3;
            for(int d = 0; d < 3; d++){
                p->force[d] += monoF * dx[d];
            }
            // Quadrupole
            double c1 = -3.0 * G * p->mass / r5;
            double c2 = 7.5 * G * p->mass / r7;
            double Q_r[3] = {0.0, 0.0, 0.0};
            for(int i = 0; i < 3; i++){
                Q_r[i] = cell->quad[i][0] * dx[0] + cell->quad[i][1] * dx[1] + cell->quad[i][2] * dx[2];
            }
            double r_Q_r = dx[0] * Q_r[0] + dx[1] * Q_r[1] + dx[2] * Q_r[2];
            for(int k = 0; k < 3; k++){
                p->force[k] += c1 * Q_r[k] + c2 * r_Q_r * dx[k];
            }
            
        }
        else{
            // Too close: recurse downward
            for(int c = 0; c < 8; c++){
                if(cell->children[c]){
                    local_count += evaluateForceOnParticle(p, cell->children[c]);
                }
            }
        }
    }
    return local_count;
}

void freeTree(Cell* cell){
    if(!cell) return;
    if(cell->particles) free(cell->particles);
    for(int c = 0; c < 8; c++) freeTree(cell->children[c]);
    free(cell);
}

int main(){

    printf("Start Input\n");
    fflush(stdout);
    start = omp_get_wtime();
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
        else if(lcount == 6){
            if(sscanf(line, "%d", &NUM_OF_THREADS));
        }
    }
    fclose(fp);
    omp_set_num_threads(NUM_OF_THREADS);

    // Input particle data
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
    double* buffer = (double*)malloc(file_size);
    if(!buffer){
        printf("Memory allocate failed\n");
        fclose(fptr);
        exit(1);
    }
    if(fread(buffer, sizeof(double), num_of_doubles, fptr)){};

    int N = num_of_doubles / 7;
    Particle** plist = (Particle**)malloc(N * sizeof(Particle*));
    Particle* particles = (Particle*)malloc(N * sizeof(Particle));
    for(int i = 0; i < N; i++){
        particles[i].mass = buffer[i * 7];
        particles[i].x[0] = buffer[i * 7 + 1];
        particles[i].x[1] = buffer[i * 7 + 2];
        particles[i].x[2] = buffer[i * 7 + 3];
        particles[i].force[0] = particles[i].force[1] = particles[i].force[2] = 0;
    }
    for(int i = 0; i < N; i++) plist[i] = &particles[i];
    end = omp_get_wtime();
    printf("Finish Input\n");
    double input_time = (end - start) * 1000;
    printf("Input time: %.3lf ms\n\n", input_time);

    // Build Octree
    printf("Start Building Octree\n");
    fflush(stdout);
    start = omp_get_wtime();
    Cell* root = createCell(plist, N);
    subdivideCell(root);
    end = omp_get_wtime();
    printf("Finish Building Octree\n");
    double build_time = (end - start) * 1000;
    printf("Building time: %.3lf ms\n\n", build_time);

    // Compute Multipole Expansion for each cell
    printf("Start Multipole Expansion\n");
    fflush(stdout);
    start = omp_get_wtime();
    computeMassDistribution(root);
    end = omp_get_wtime();
    printf("Finish Multipole Expansion\n");
    double expansion_time = (end - start) * 1000;
    printf("Expansion time: %.3lf ms\n\n", expansion_time);

    // Evaluate Force on every particles
    printf("Start Force Evaluation\n");
    fflush(stdout);
    start = omp_get_wtime();
    #pragma omp parallel for schedule(guided, 16) reduction(+:direct_count)
    for(int i = 0; i < N; i++){
        direct_count += evaluateForceOnParticle(&particles[i], root);
    }
    end = omp_get_wtime();
    printf("Finish Force Evaluation\n");
    double eval_time = (end - start) * 1000;
    printf("Evaluation time: %.3lf ms\n\n", eval_time);

    printf("Direct eval: %llu\n", direct_count);
    printf("Direct ratio: %.5f %%\n\n", direct_count * 100.0 / N / (N - 1));

    double total_time = build_time + expansion_time + eval_time;
    printf("Total execution time: %.3lf ms\n\n", total_time);

    // Output the force vectors
    printf("Start Output\n");
    start = omp_get_wtime();
    FILE* fcsv;
    fcsv = fopen(outfile, "w");
    if(!fcsv){
        printf("Failed to open output file!\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(fcsv, "%.16e,%.16e,%.16e\n", particles[i].force[0], particles[i].force[1], particles[i].force[2]);
    }
    fclose(fcsv);

    free(plist);
    freeTree(root);
    free(particles);
    free(buffer);
    end = omp_get_wtime();
    printf("Finish Output\n");
    double output_time = (end - start) * 1000;
    printf("Output time: %.3lf ms\n\n", output_time);
    return 0;
}