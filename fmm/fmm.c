#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define G 1.0                    // Gravitational constant
#define EPSILON 1.0e-9           // Softening factor to prevent division by zero
#define INIT_LEN 32
#define MAX_PARTICLES_PER_LEAF 10
#define MAX_DEPTH 20
#define HALF_BOX_SIZE 5.0

typedef struct {
    double mass;
    double pos[3];   // x, y, z
} Particle;

typedef struct{
    Particle* arr;
    int size;
    int capacity;
} pList;

typedef struct octreenode{
    double center[3];
    double half_size;
    pList* particle_list;
    struct octreenode* children_list;
    bool is_leaf;
} OctreeNode;

int N;           // Number of particles

// Basic operations of particle lists
pList* pListInit(int init_len){
    pList* newList = (pList*)malloc(sizeof(pList));
    newList->capacity = init_len;
    newList->size = 0;
    newList->arr = (Particle*)malloc(init_len * sizeof(Particle));
    return newList;
}

void pListInsert(pList* list, Particle particle){
    list->arr[list->size] = particle;
    list->size++;
    if(list->size > list->capacity){
        list->capacity *= 2;
        list->arr = (Particle*)realloc(list->arr, list->capacity * sizeof(Particle));
    }
}

// Compute which octree child node the particle should go
int ComputeOctreeIndex(Particle particle, double* center){
    int x = particle.pos[0] > center[0] ? 1 : 0;
    int y = particle.pos[1] > center[1] ? 1 : 0;
    int z = particle.pos[2] > center[2] ? 1 : 0;
    return x + 2 * y + 4 * z;
}

// Compute the center of a child octree node from its parent's center
double* ComputeChildCenter(double* center, double half_size, int index){
    double offset = half_size / 2;
    int x = (index & 1) == 1 ? 1 : -1;
    int y = ((index / 2) & 1) == 1 ? 1 : -1;
    int z = ((index / 4) & 1) == 1 ? 1 : -1;
    double* child_center = (double*)malloc(3 * sizeof(double));
    child_center[0] = center[0] + offset * x;
    child_center[1] = center[1] + offset * y;
    child_center[2] = center[2] + offset * z;
    return child_center;
}

// Building the octree of given particles
OctreeNode* BuildOctree(pList* particle_list, double* center, double half_size, int depth){
    OctreeNode* node = (OctreeNode*)malloc(sizeof(OctreeNode));
    for(int i = 0; i < 3; i++){
        node->center[i] = center[i];
    }
    node->half_size = half_size;
    node->particle_list = particle_list;
    node->children_list = NULL;
    node->is_leaf = true;

    if(particle_list->size < MAX_PARTICLES_PER_LEAF || depth > MAX_DEPTH){
        // Is a leaf node
        printf("construct a leaf node with %d particles\n", particle_list->size);
        for(int i = 0; i < particle_list->size; i++){
            printf("%.10e\n", particle_list->arr[i].mass);
        }
        return node;
    }

    // Is an internal node, split children
    node->is_leaf = false;
    node->children_list = (OctreeNode*)malloc(8 * sizeof(OctreeNode));
    pList** children_particles = (pList**)malloc(8 * sizeof(pList*));
    for(int i = 0; i < 8; i++){
        children_particles[i] = pListInit(INIT_LEN);
    }

    for(int i = 0; i < particle_list->size; i++){
        int id = ComputeOctreeIndex(particle_list->arr[i], center);
        pListInsert(children_particles[id], particle_list->arr[i]);
    }

    for(int i = 0; i < 8; i++){
        double* child_center = ComputeChildCenter(center, half_size, i);
        double child_half_size = half_size / 2.0;
        BuildOctree(children_particles[i], child_center, child_half_size, depth + 1);
    }
    printf("construct an internal node with %d nodes\n", particle_list->size);
    return node;
}

void FMM_compute_gravity(Particle particles[], double forces[][3]){
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

    pList* particles = (pList*)malloc(sizeof(pList));
    particles->arr = (Particle*)malloc(N * sizeof(Particle));
    particles->capacity = N;
    particles->size = N;
    for(int i = 0; i < N; i++){
        particles->arr[i].mass = buffer[i * 7];
        particles->arr[i].pos[0] = buffer[i * 7 + 1];
        particles->arr[i].pos[1] = buffer[i * 7 + 2];
        particles->arr[i].pos[2] = buffer[i * 7 + 3];
    }
    fclose(fptr);

    double forces[N][3];

    double init_center[] = {0.0, 0.0, 0.0};
    OctreeNode* root = BuildOctree(particles, init_center, HALF_BOX_SIZE, 0);







    //FMM_compute_gravity(particles, forces);

    /*
    // Output the force vectors
    FILE* fcsv;
    fcsv = fopen("force.csv", "w");
    if(!fcsv){
        printf("Failed to open output file!\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fprintf(fcsv, "%.10e,%.10e,%.10e\n",forces[i][0], forces[i][1], forces[i][2]);
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
    */
    return 0;
}
