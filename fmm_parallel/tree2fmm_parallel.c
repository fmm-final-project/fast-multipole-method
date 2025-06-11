#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define MAX_LINE_LEN 1024
#define EPSILON 1e-12

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
    int* pids;
    struct Cell* children[8];
    omp_lock_t lock;
} Cell;

// Input parameters
double G;
double THETA;
int MAX_PARTICLES_PER_CELL;
char datafile[MAX_LINE_LEN];
char outfile[MAX_LINE_LEN];
int THRESHOLD;

// Particle arrays
double* mass;
double** pos;
double** force;

// Parallelization
int NUM_OF_THREADS;
unsigned long long direct_count = 0;
unsigned long long* count_array;
double*** partial_force;
double start, end;

// Create a new cell, compute its center and size from particle extents
Cell* createCell(int* pids, int n){
    Cell* cell = (Cell*)malloc(sizeof(Cell));
    cell->nParticles = n;
    cell->pids = (int*)malloc(n * sizeof(int));
    for(int i = 0; i < n; i++){
        cell->pids[i] = pids[i];
    }
    // Compute bounding box
    double min[3], max[3];
    for(int d = 0; d < 3; d++){
        min[d] = max[d] = pos[pids[0]][d];
    }
    for(int i = 1; i < n; i++){
        for(int d = 0; d < 3; d++){
            double v = pos[pids[i]][d];
            if(v < min[d]) min[d] = v;
            if(v > max[d]) max[d] = v;
        }
    }
    // Set center and size
    for(int d = 0; d < 3; d++){
        cell->center[d] = 0.5 * (min[d] + max[d]);
    }
    double dx = max[0] - min[0];
    double dy = max[1] - min[1];
    double dz = max[2] - min[2];
    cell->size = fmax(fmax(dx, dy), dz);
    // Initialize children
    for(int c = 0; c < 8; c++) cell->children[c] = NULL;
    omp_init_lock(&cell->lock);
    return cell;
}

void subdivideCell(Cell *cell){
    if(cell->nParticles <= MAX_PARTICLES_PER_CELL) return;

    // Distribute particles into octants
    int* plist = cell->pids;
    int counts[8] = {0};
    for(int i = 0; i < cell->nParticles; i++){
        int p = plist[i];
        int idx = (pos[p][0] > cell->center[0])
                + 2 * (pos[p][1] > cell->center[1])
                + 4 * (pos[p][2] > cell->center[2]);
        counts[idx]++;
    }
    int** lists = (int**)malloc(8 * sizeof(int*));
    for(int c = 0; c < 8; c++){
        lists[c] = (int*)malloc(counts[c] * sizeof(int));
        counts[c] = 0;
    }
    for(int i = 0; i < cell->nParticles; i++){
        int p = plist[i];
        int idx = (pos[p][0] > cell->center[0])
                + 2 * (pos[p][1] > cell->center[1])
                + 4 * (pos[p][2] > cell->center[2]);
        lists[idx][counts[idx]++] = p;
    }
    free(cell->pids);
    for(int c = 0; c < 8; c++){
        if(counts[c] > 0){
            cell->children[c] = createCell(lists[c], counts[c]);
            subdivideCell(cell->children[c]);
        }
        free(lists[c]);
    }
    free(lists);
    cell->pids = NULL;
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
    if(cell->pids != NULL){
        // Leaf: compute multipoles directly
        for(int i = 0; i < cell->nParticles; i++){
            int p = cell->pids[i];
            cell->mass += mass[p];
            for(int d = 0; d < 3; d++) cell->massCenter[d] += mass[p] * pos[p][d];
        }
        for(int d = 0; d < 3; d++) cell->massCenter[d] /= cell->mass;
        // Quadrupole
        for(int i = 0; i < cell->nParticles; i++){
            int p = cell->pids[i];
            double dx[3], r2 = EPSILON;
            for(int d = 0; d < 3; d++){
                dx[d] = pos[p][d] - cell->massCenter[d];
                r2 += dx[d] * dx[d];
            }
            for(int a = 0; a < 3; a++){
                for(int b = 0; b < 3; b++){
                    cell->quad[a][b] += mass[p] * (3 * dx[a] * dx[b] - (a==b ? r2 : 0.0));
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
    int nParticles = A->nParticles + B->nParticles;
    // Self-Interaction -> Subnode pairs
    if(A == B){
        if(A->pids != NULL) return;
        for(int i = 0; i < 8; i++){
            for(int j = i; j < 8; j++){
                if(A->children[i] && A->children[j]){
                    #pragma omp task firstprivate(i, j) final(nParticles < THRESHOLD)
                    dualTreeWalk(A->children[i], A->children[j]);
                }
            }
        }
        #pragma omp taskwait
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
        //double r7 = r2 * r5;
        //double c1 = -3.0 / r5;
        //double c2 = 7.5 / r7;

        // New 25 flops
        /*
        double QA_r[3] = {0.0, 0.0, 0.0};
        for(int i = 0; i < 3; i++){
            QA_r[i] = A->quad[i][0] * dx[0] + A->quad[i][1] * dx[1] + A->quad[i][2] * dx[2];
        }
        double r_QA_r = dx[0] * QA_r[0] + dx[1] * QA_r[1] + dx[2] * QA_r[2];
        */
        omp_set_lock(&B->lock);
        for(int i = 0; i < 3; i++){
            // Compute L1
            B->L1[i] += -A->mass * dx[i] / r3;
            //B->L1[i] += c1 * QA_r[i] + c2 * r_QA_r * dx[i];
            // Compute L2
            for(int j = 0; j < 3; j++){
                B->L2[i][j] += -A->mass * (3 * dx[i] * dx[j] / r5 - (i == j ? 1.0/r3 : 0));
            }
        }
        omp_unset_lock(&B->lock);

        /*
        double QB_r[3] = {0.0, 0.0, 0.0};
        for(int i = 0; i < 3; i++){
            QB_r[i] = B->quad[i][0] * (-dx[0]) + B->quad[i][1] * (-dx[1]) + B->quad[i][2] * (-dx[2]);
        }
        double r_QB_r = (-dx[0]) * QB_r[0] + (-dx[1]) * QB_r[1] + (-dx[2]) * QB_r[2];
        */
        omp_set_lock(&A->lock);
        for(int i = 0; i < 3; i++){
            // Compute L1
            A->L1[i] += -B->mass * (-dx[i]) / r3;
            //A->L1[i] += c1 * QB_r[i] + c2 * r_QB_r * (-dx[i]);
            // Compute L2
            for(int j = 0; j < 3; j++){
                A->L2[i][j] += -B->mass * (3 * (-dx[i]) * (-dx[j]) / r5 - (i == j ? 1.0/r3 : 0));
            }
        }
        omp_unset_lock(&A->lock);

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
    if(A->pids == NULL && B->pids == NULL){
        // Open the larger one
        if(A->size >= B->size){
            for(int i = 0; i < 8; i++){
                if(!A->children[i]) continue;
                #pragma omp task firstprivate(i) final(nParticles < THRESHOLD)
                dualTreeWalk(A->children[i], B);
            }
        }
        else{
            for(int i = 0; i < 8; i++){
                if(!B->children[i]) continue;
                #pragma omp task firstprivate(i) final(nParticles < THRESHOLD)
                dualTreeWalk(A, B->children[i]);
            }
        }
        #pragma omp taskwait
    }
    else if(A->pids == NULL){
        // Only A can open
        for(int i = 0; i < 8; i++) {
            if(!A->children[i]) continue;
            #pragma omp task firstprivate(i) final(nParticles < THRESHOLD)
            dualTreeWalk(A->children[i], B);
        }
        #pragma omp taskwait
    }
    else if(B->pids == NULL){
        // Only B can open
        for(int i = 0; i < 8; i++){
            if(!B->children[i]) continue;
            #pragma omp task firstprivate(i) final(nParticles < THRESHOLD)
            dualTreeWalk(A, B->children[i]);
        }
        #pragma omp taskwait
    }
    else{
        // direct sum between every particle in A and every in B
        int tid = omp_get_thread_num();
        for(int i = 0; i < A->nParticles; i++){
            int p = A->pids[i];
            for(int j = 0; j < B->nParticles; j++){
                int q = B->pids[j];
                count_array[tid] += 2;
                double dx[3], r2 = EPSILON;
                for(int d = 0; d < 3; d++){
                    dx[d] = pos[p][d] - pos[q][d];
                    r2 += dx[d] * dx[d];
                }
                double inv3 = 1.0 / (r2 * sqrt(r2));
                double f = -mass[p] * mass[q] * inv3;
                for(int d = 0; d < 3; d++){
                    partial_force[tid][p][d] += f * dx[d];
                    partial_force[tid][q][d] += -f * dx[d];
                }
            }
        }
    }
}

void localPassDown(Cell* cell){
    if(cell->pids == NULL){
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
            int p = cell->pids[i];
            // Direct summation
            for(int j = i + 1; j < cell->nParticles; j++){
                int q = cell->pids[j];
                direct_count += 2;
                double dx[3], r2 = EPSILON;
                for(int d = 0; d < 3; d++){
                    dx[d] = pos[p][d] - pos[q][d];
                    r2 += dx[d] * dx[d];
                }
                double inv3 = 1.0 / (r2 * sqrt(r2));
                double f = -mass[p] * mass[q] * inv3;
                for(int d = 0; d < 3; d++){
                    force[p][d] += f * dx[d];
                    force[q][d] -= f * dx[d];
                }
            }
            // From partial sums of each threads
            for(int i = 0; i < NUM_OF_THREADS; i++){
                for(int d = 0; d < 3; d++){
                    force[p][d] += partial_force[i][p][d];
                }
            }

            // From local expansion
            double dx[3];
            for(int d = 0; d < 3; d++){
                dx[d] = pos[p][d] - cell->massCenter[d];
            }
            for(int d = 0; d < 3; d++){
                force[p][d] += mass[p] * cell->L1[d];
                for(int b = 0; b < 3; b++){
                    force[p][d] += mass[p] * cell->L2[d][b] * dx[b];
                }
            }
        }
    }
}

void freeTree(Cell* cell){
    if(!cell) return;
    if(cell->pids) free(cell->pids);
    omp_destroy_lock(&cell->lock);
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
    THRESHOLD = MAX_PARTICLES_PER_CELL * 1000;

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

    // Initialize particle arrays
    mass = (double*)calloc(N, sizeof(double));
    pos = (double**)malloc(N * sizeof(double*));
    pos[0] = (double*)calloc(N * 3, sizeof(double));
    for(int i = 1; i < N; i++){
        pos[i] = pos[0] + 3 * i;
    }
    force = (double**)malloc(N * sizeof(double*));
    force[0] = (double*)calloc(N * 3, sizeof(double));
    for(int i = 1; i < N; i++){
        force[i] = force[0] + 3 * i;
    }
    for(int i = 0; i < N; i++){
        mass[i] = buffer[i * 7];
        pos[i][0] = buffer[i * 7 + 1];
        pos[i][1] = buffer[i * 7 + 2];
        pos[i][2] = buffer[i * 7 + 3];
        force[i][0] = force[i][1] = force[i][2] = 0.0;
    }
    
    // Initialize parallelization arrays
    omp_set_num_threads(NUM_OF_THREADS);
    count_array = (unsigned long long*)calloc(NUM_OF_THREADS, sizeof(unsigned long long));
    partial_force = (double***)malloc(NUM_OF_THREADS * sizeof(double**));
    for(int i = 0; i < NUM_OF_THREADS; i++){
        partial_force[i] = (double**)malloc(N * sizeof(double*));
        partial_force[i][0] = (double*)calloc(N * 3, sizeof(double));
        for(int j = 1; j < N; j++){
            partial_force[i][j] = partial_force[i][0] + 3 * j;
        }
    }
    end = omp_get_wtime();
    printf("Finish Input\n");
    double input_time = (end - start) * 1000;
    printf("Input time: %.3lf ms\n\n", input_time);

    // Build Octree
    printf("Start Building Octree\n");
    fflush(stdout);
    start = omp_get_wtime();
    int* plist = (int*)malloc(N * sizeof(int));
    for(int i = 0; i < N; i++) plist[i] = i;
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
    #pragma omp parallel
    {
        #pragma omp single
        {
            dualTreeWalk(root, root);
        }
    }
    end = omp_get_wtime();
    printf("Finish Dual Tree Walk\n");
    double dual_time = (end - start) * 1000;
    printf("Dual Tree Walk time: %.3lf ms\n\n", dual_time);

    // Compute L2L and evaluation
    printf("Start Local Passing Down\n");
    fflush(stdout);
    start = omp_get_wtime();
    localPassDown(root);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < N; i++){
        for(int d = 0; d < 3; d++){
            force[i][d] *= G;
        }
    }
    end = omp_get_wtime();
    printf("Finish Local Passing Down\n");
    double local_time = (end - start) * 1000;
    printf("Local Passing Down time: %.3lf ms\n\n", local_time);

    for(int i = 0; i < NUM_OF_THREADS; i++){
        direct_count += count_array[i];
    }
    printf("Direct count %llu\n", direct_count);
    printf("Direct Ratio = %f %%\n\n", direct_count * 100.0 / N / (N - 1));

    double total_time = build_time + expansion_time + dual_time + local_time;
    printf("Number of Treads: %d\n", NUM_OF_THREADS);
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
        fprintf(fcsv, "%.10e,%.10e,%.10e\n", force[i][0], force[i][1], force[i][2]);
    }
    fclose(fcsv);

    free(plist);
    freeTree(root);
    free(buffer);
    free(count_array);
    for(int i = 0; i < NUM_OF_THREADS; i++){
        free(partial_force[i][0]);
        free(partial_force[i]);
    }
    free(mass);
    free(pos[0]);
    free(pos);
    free(force[0]);
    free(force);
    end = omp_get_wtime();
    printf("Finish Output\n");
    double output_time = (end - start) * 1000;
    printf("Output time: %.3lf ms\n", output_time);
    return 0;
}