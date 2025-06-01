#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define MAX_LINE_LEN 1024
#define EPSILON 1e-12

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
    double L1[3];         // 1st order local expansion
    double L2[3][3];      // 2nd order local expansion
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

// Performance
unsigned long long direct_count = 0;
double start, end;

// Create a new cell, compute its center and size from particle extents
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
        cell->L1[i] = 0.0;
        for(int j = 0; j < 3; j++){
            cell->quad[i][j] = 0.0;
            cell->L2[i][j] = 0.0;
        }
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
            double dx[3], r2 = EPSILON;
            for (int d = 0; d < 3; d++){
                dx[d] = p->x[d] - cell->massCenter[d];
                r2 += dx[d] * dx[d];
            }
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    cell->quad[a][b] += p->mass * (3 * dx[a] * dx[b] - (a==b ? r2 : 0.0));
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
            double dC[3], r2 = EPSILON;
            for(int d = 0; d < 3; d++){
                dC[d] = ch->massCenter[d] - cell->massCenter[d];
                r2 += dC[d] * dC[d];
            }
            for(int a = 0; a < 3; a++){
                for(int b = 0; b < 3; b++){
                    cell->quad[a][b] += ch->quad[a][b] + ch->mass * (3 * dC[a] * dC[b] - (a==b ? r2 : 0.0));
                }
            }
        }
    }
}

void dualTreeWalk(Cell* A, Cell* B){
    // Self-Interaction -> Subnode pairs
    if(A == B){
        if(A->nParticles > 0) return;
        for(int i = 0; i < 8; i++){
            for(int j = i; j < 8; j++){
                if(A->children[i] && A->children[j]){
                    dualTreeWalk(A->children[i], A->children[j]);
                }
            }
        }
        return;
    }

    double dx[3], r2 = EPSILON;
    for(int d = 0; d < 3; d++){
        dx[d] = B->massCenter[d] - A->massCenter[d];
        r2 += dx[d] * dx[d];
    }

    // Well-separated
    double S = A->size + B->size;
    if(S * S < r2 * THETA * THETA){
        double r = sqrt(r2);
        double r3 = r2 * r;
        double r5 = r2 * r3;
        double r7 = r2 * r5;
        double c1 = -3.0 / r5;
        double c2 = 7.5 / r7;

        // Compute L1
        for(int i = 0; i < 3; i++){
            B->L1[i] += -A->mass * dx[i] / r3;
            A->L1[i] += -B->mass * (-dx[i]) / r3;
        }

        // New 25 flops
        double QA_r[3] = {0.0, 0.0, 0.0};
        for(int i = 0; i < 3; i++){
            QA_r[i] = A->quad[i][0] * dx[0] + A->quad[i][1] * dx[1] + A->quad[i][2] * dx[2];
        }
        double r_QA_r = dx[0] * QA_r[0] + dx[1] * QA_r[1] + dx[2] * QA_r[2];
        for(int k = 0; k < 3; k++) B->L1[k] += c1 * QA_r[k] + c2 * r_QA_r * dx[k];

        double QB_r[3] = {0.0, 0.0, 0.0};
        for(int i = 0; i < 3; i++){
            QB_r[i] = B->quad[i][0] * (-dx[0]) + B->quad[i][1] * (-dx[1]) + B->quad[i][2] * (-dx[2]);
        }
        double r_QB_r = (-dx[0]) * QB_r[0] + (-dx[1]) * QB_r[1] + (-dx[2]) * QB_r[2];
        for(int k = 0; k < 3; k++) A->L1[k] += c1 * QB_r[k] + c2 * r_QB_r * (-dx[k]);

        // Compute L2
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                B->L2[i][j] += -A->mass * (3 * dx[i] * dx[j] / r5 - (i == j ? 1.0/r3 : 0));
                A->L2[i][j] += -B->mass * (3 * (-dx[i]) * (-dx[j]) / r5 - (i == j ? 1.0/r3 : 0));
            }
        }
        /*
        double r4 = r2 * r2;
        // Quadrupole contribution to L2
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int m = 0; m < 3; m++) {
                    for (int n = 0; n < 3; n++) {
                        // Coefficients for B's L2
                        double coeffB = ((i == m && j == n) + (i == n && j == m)
                                        - 5.0 * dx[i] * dx[j] * dx[m] * dx[n] / r4) / r5;
                        B->L2[i][j] += -1.5 * A->quad[m][n] * coeffB;

                        // Coefficients for A's L2 (reverse dx)
                        double coeffA = ((i == m && j == n) + (i == n && j == m)
                                        - 5.0 * (-dx[i]) * (-dx[j]) * (-dx[m]) * (-dx[n]) / r4) / r5;
                        A->L2[i][j] += -1.5 * B->quad[m][n] * coeffA;
                    }
                }
            }
        }
        */
        return;
    }

    // Split into subnodes for new pairs
    if(A->nParticles == 0 && B->nParticles == 0){
        // Open the larger one
        if(A->size >= B->size){
            for(int i = 0; i < 8; i++){
                if(!A->children[i]) continue;
                dualTreeWalk(A->children[i], B);
            }
        }
        else{
            for(int i = 0; i < 8; i++){
                if(!B->children[i]) continue;
                dualTreeWalk(A, B->children[i]);
            }
        }
    }
    else if(A->nParticles == 0){
        // Only A can open
        for(int i = 0; i < 8; i++) {
            if(!A->children[i]) continue;
            dualTreeWalk(A->children[i], B);
        }
    }
    else if(B->nParticles == 0){
        // Only B can open
        for(int j = 0; j < 8; j++){
            if(!B->children[j]) continue;
            dualTreeWalk(A, B->children[j]);
        }
    }
    else{
        // direct sum between every particle in A and every in B
        for(int i = 0; i < A->nParticles; i++){
            Particle *p = A->particles[i];
            for(int j = 0; j < B->nParticles; j++){
                Particle *q = B->particles[j];
                direct_count += 2;
                double dx[3], r2 = EPSILON;
                for(int d = 0; d < 3; d++){
                    dx[d] = p->x[d] - q->x[d];
                    r2 += dx[d] * dx[d];
                }
                double inv3 = 1.0 / (r2 * sqrt(r2));
                double f = -p->mass * q->mass * inv3;
                for(int d = 0; d < 3; d++){
                    p->force[d] += f * dx[d];
                    q->force[d] += -f * dx[d];
                }
            }
        }
    }
}

void localPassDown(Cell* cell){
    if(cell->nParticles == 0){
        for(int i = 0; i < 8; i++){
            if(cell->children[i]){
                Cell* ch = cell->children[i];
                double dx[3];
                for(int a = 0; a < 3; a++){
                    dx[a] = ch->massCenter[a] - cell->massCenter[a];
                    cell->children[i]->L1[a] += cell->L1[a];
                }

                for(int a = 0; a < 3; a++){
                    for(int b = 0; b < 3; b++){
                        ch->L1[a] += cell->L2[a][b] * dx[b];
                        ch->L2[a][b] += cell->L2[a][b];
                    }
                }
                localPassDown(ch);
            }
        }
    }
    else{
        for(int i = 0; i < cell->nParticles; i++){
            Particle* p = cell->particles[i];
            // Direct summation
            for(int j = i + 1; j < cell->nParticles; j++){
                Particle* q = cell->particles[j];
                direct_count += 2;
                double dx[3], r2 = EPSILON;
                for(int d = 0; d < 3; d++){
                    dx[d] = p->x[d] - q->x[d];
                    r2 += dx[d] * dx[d];
                }
                double inv3 = 1.0 / (r2 * sqrt(r2));
                double f = -p->mass * q->mass * inv3;
                for(int d = 0; d < 3; d++){
                    p->force[d] += f * dx[d];
                    q->force[d] -= f * dx[d];
                }
            }

            // From local expansion
            double dx[3];
            for(int d = 0; d < 3; d++){
                dx[d] = p->x[d] - cell->massCenter[d];
            }
            for(int d = 0; d < 3; d++){
                p->force[d] += p->mass * cell->L1[d];
                for(int b = 0; b < 3; b++){
                    p->force[d] += p->mass * cell->L2[d][b] * dx[b];
                }
            }
        }
    }
}

void freeTree(Cell* cell){
    if(!cell) return;
    if(cell->particles) free(cell->particles);
    for(int c = 0; c < 8; c++) freeTree(cell->children[c]);
    free(cell);
}

int main(){
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

    int N = num_of_doubles / 7;

    double* buffer = (double*)malloc(file_size);
    if(!buffer){
        printf("Memory allocate failed\n");
        fclose(fptr);
        exit(1);
    }
    if(fread(buffer, sizeof(double), num_of_doubles, fptr)){};

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

    // Compute local expansion using dual tree walk
    printf("Start Dual Tree Walk\n");
    fflush(stdout);
    start = omp_get_wtime();
    dualTreeWalk(root, root);
    end = omp_get_wtime();
    printf("Finish Dual Tree Walk\n");
    double dual_time = (end - start) * 1000;
    printf("Dual Tree Walk time: %.3lf ms\n\n", dual_time);

    // Compute L2L and evaluation
    printf("Start Local Passing Down\n");
    fflush(stdout);
    start = omp_get_wtime();
    localPassDown(root);
    for(int i = 0; i < N; i++){
        for(int d = 0; d < 3; d++){
            particles[i].force[d] *= G;
        }
    }
    end = omp_get_wtime();
    printf("Finish Local Passing Down\n");
    double local_time = (end - start) * 1000;
    printf("Local Passing Down time: %.3lf ms\n\n", local_time);
    printf("Direct count %llu\n", direct_count);
    printf("Direct Ratio = %f %%\n\n", direct_count * 100.0 / N / (N - 1));

    double total_time = build_time + expansion_time + dual_time + local_time;
    printf("Total execution time: %.3lf ms\n\n", total_time);

    // Output the force vectors
    printf("Start Output\n");
    FILE* fcsv;
    fcsv = fopen(outfile, "w");
    if(!fcsv){
        printf("Failed to open output file!\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(fcsv, "%.10e,%.10e,%.10e\n", particles[i].force[0], particles[i].force[1], particles[i].force[2]);
    }
    fclose(fcsv);
    printf("Finish Output\n");

    free(plist);
    freeTree(root);
    free(particles);
    free(buffer);
    return 0;
}