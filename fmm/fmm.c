#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define G 1.0                      // Gravitational constant
#define EPSILON 1.0e-9             // Softening factor to prevent division by zero

#define INIT_LEN 16
#define MAX_PARTICLES_PER_LEAF 10
#define MAX_DEPTH 20
#define HALF_BOX_SIZE 5.0
#define DIM 3
#define THETA 0.5

typedef struct {
    double mass;
    double pos[DIM];   // x, y, z
} Particle;

typedef struct{
    Particle* arr;
    int size;
    int capacity;
} pList;

typedef struct octreenode{
    unsigned long long id;
    double center[DIM];
    double half_size;
    pList* particle_list;
    struct octreenode** children_list;

    // Multipole expansions
    double monopole;
    double dipole[DIM];
    double quadrupole[DIM][DIM];

    // Local expansions
    double local_monopole;
    double local_dipole[DIM];
    double local_quadrupole[DIM][DIM];
} OctreeNode;

typedef struct{
    OctreeNode** arr;
    int size;
    int capacity;
} oList;

int N;            // Number of particles

// Basic operations of particle lists
pList* pListInit(){
    pList* newList = (pList*)malloc(sizeof(pList));
    newList->capacity = INIT_LEN;
    newList->size = 0;
    newList->arr = (Particle*)malloc(INIT_LEN * sizeof(Particle));
    return newList;
}

void pListInsert(pList* list, Particle particle){
    list->arr[list->size] = particle;
    list->size++;
    if(list->size == list->capacity){
        list->capacity *= 2;
        list->arr = (Particle*)realloc(list->arr, list->capacity * sizeof(Particle));
    }
}

// Basic operations of octree node lists
oList* oListInit(){
    oList* newList = (oList*)malloc(sizeof(oList));
    newList->capacity = INIT_LEN;
    newList->size = 0;
    newList->arr = (OctreeNode**)malloc(INIT_LEN * sizeof(OctreeNode*));
    return newList;
}

void oListInsert(oList* list, OctreeNode* node){
    list->arr[list->size] = node;
    list->size++;
    if(list->size == list->capacity){
        list->capacity *= 2;
        list->arr = (OctreeNode**)realloc(list->arr, list->capacity * sizeof(OctreeNode*));
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
    double* child_center = (double*)malloc(DIM * sizeof(double));
    child_center[0] = center[0] + offset * x;
    child_center[1] = center[1] + offset * y;
    child_center[2] = center[2] + offset * z;
    return child_center;
}

// Building the octree of given particles
OctreeNode* BuildOctree(pList* particle_list, double* center, double half_size, int depth, unsigned long long id){
    OctreeNode* node = (OctreeNode*)calloc(1, sizeof(OctreeNode));
    for(int i = 0; i < DIM; i++){
        node->center[i] = center[i];
    }
    node->half_size = half_size;
    node->particle_list = particle_list;
    node->children_list = NULL;
    node->id = id;

    if(particle_list->size < MAX_PARTICLES_PER_LEAF || depth > MAX_DEPTH){
        // Is a leaf node
        printf("Construct a leaf node id = %llu with %d particles\n", node->id, particle_list->size);
        for(int i = 0; i < particle_list->size; i++){
            //printf("%.10e\n", particle_list->arr[i].mass);
        }
        return node;
    }

    // Is an internal node, split children
    node->children_list = (OctreeNode**)malloc(8 * sizeof(OctreeNode*));
    pList** children_particles = (pList**)malloc(8 * sizeof(pList*));
    for(int i = 0; i < 8; i++){
        children_particles[i] = pListInit();
    }

    for(int i = 0; i < particle_list->size; i++){
        int id = ComputeOctreeIndex(particle_list->arr[i], center);
        pListInsert(children_particles[id], particle_list->arr[i]);
    }

    for(int i = 0; i < 8; i++){
        if(children_particles[i]->size > 0){
            double* child_center = ComputeChildCenter(center, half_size, i);
            double child_half_size = half_size / 2.0;
            node->children_list[i] = BuildOctree(children_particles[i], child_center, child_half_size, depth + 1, id * 8 + i);
        }
        else{
            node->children_list[i] = NULL;
        }
    }
    printf("Construct an internal node id = %llu with %d nodes\n", node->id, particle_list->size);
    return node;
}

// Compute the multipole expansions for each octree node
void ComputeMultipoles(OctreeNode* node){
    //printf("visit node with center x = %.10e\n", node->center[0]);

    if(node->children_list == NULL){
        for(int i = 0; i < node->particle_list->size; i++){
            Particle particle = node->particle_list->arr[i];
            double m = particle.mass;
            double dx[DIM];
            for(int j = 0; j < DIM; j++){
                dx[j] = particle.pos[j] - node->center[j];
            }

            // monopole
            node->monopole += m;

            // dipole
            for(int j = 0; j < DIM; j++){
                node->dipole[j] += m * dx[j];
            }

            // quadrupole
            for(int row = 0; row < DIM; row++){
                for(int col = 0; col < DIM; col++){
                    node->quadrupole[row][col] += m * dx[row] * dx[col];
                }
            }
        }
        printf("leaf id = %llu multipole finished, monopole = %.10e\n", node->id, node->monopole);
    }
    else{ // M2M translation
        for(int i = 0; i < 8; i++){
            OctreeNode* child = node->children_list[i];
            if(child == NULL) continue;

            ComputeMultipoles(child);
            double m = child->monopole;
            double dx[DIM];
            for(int j = 0; j < DIM; j++){
                dx[j] = child->center[j] - node->center[j];
            }

            // monopole
            node->monopole += m;

            // dipole
            for(int j = 0; j < DIM; j++){
                node->dipole[j] += child->dipole[j] + m * dx[j];
            }

            // quadrupole
            for(int row = 0; row < DIM; row++){
                for(int col = 0; col < DIM; col++){
                    node->quadrupole[row][col] += child->quadrupole[row][col]
                                                + dx[row] * child->dipole[col]
                                                + child->dipole[row] * dx[col];
                                                + m * dx[row] * dx[col];
                }
            }

        }
        printf("internal node id = %llu multipole finished, monopole = %.10e\n", node->id, node->monopole);
    }
}

void GetInteractionList(OctreeNode* target, OctreeNode* src, oList* list){
    if(src == NULL || src == target) return;

    double dx[3];
    double dist2 = 0;
    for(int i = 0; i < DIM; i++){
        dx[i] = src->center[i] - target->center[i];
        dist2 += dx[i] * dx[i];
    }
    double size2 = src->half_size * src->half_size * 4.0;
    if(src->children_list == NULL && size2 / dist2 > THETA * THETA){
        // Not well-seperated
        for(int i = 0; i < 8; i++){
            GetInteractionList(target, src->children_list[i], list);
        }
    }
    else{
        // Well-seperated
        oListInsert(list, src);
    }
}

void ComputeLocalExpansions(OctreeNode* node, OctreeNode* root){
    oList* interaction_list = oListInit();
    GetInteractionList(node, root, interaction_list);

    // M2L translation
    for(int i = 0; i < interaction_list->size; i++){
        OctreeNode* src = interaction_list->arr[i];

        double dx[DIM];
        double r2 = 0;
        for(int j = 0; j < DIM; j++){
            dx[j] = node->center[j] - src->center[j];
            r2 += dx[j] * dx[j];
        }
        double r = sqrt(r2);
        double rinv = 1.0 / r;
        double rinv3 = rinv * rinv * rinv;
        double rinv5 = rinv3 * rinv * rinv;

        // monopole
        node->local_monopole += src->monopole * rinv;

        // dipole
        for(int j = 0; j < DIM; j++){
            node->local_dipole[j] += -src->monopole * dx[j] * rinv3 + src->dipole[j] * rinv3;
        }

        // quadrupole
        for(int row = 0; row < DIM; row++){
            for(int col = 0; col < DIM; col++){
                node->local_quadrupole[row][col] += 3.0 * src->monopole * dx[row] * dx[col] * rinv5
                                                 - src->dipole[row] * dx[col] * rinv5
                                                 - dx[row] * src->dipole[col] * rinv5
                                                 + src->quadrupole[row][col] * rinv5;
            }
        }
    }
    printf("node id = %llu M2L finished, local monopole = %.10e\n", node->id, node->local_monopole);

    // internal nodes L2L to children
    if(node->children_list != NULL){
        for(int i = 0; i < 8; i++){
            OctreeNode* child = node->children_list[i];
            if(child == NULL) continue;

            double dx[DIM];
            for(int j = 0; j < DIM; j++){
                dx[j] = child->center[j] - node->center[j];
            }

            // monopole
            child->local_monopole = node->local_monopole;
            for(int row = 0; row < DIM; row++){
                child->local_monopole += node->local_dipole[row] * dx[row];
                for(int col = 0; col < DIM; col++){
                    child->local_monopole += 0.5 * node->local_quadrupole[row][col] * dx[row] * dx[col];
                }
            }

            // dipole
            for(int row = 0; row < DIM; row++){
                child->local_dipole[row] = node->local_dipole[row];
                for(int col = 0; col < DIM; col++){
                    child->local_dipole[row] += node->local_quadrupole[row][col] * dx[col];
                }
            }

            // quadrupole
            for(int row = 0; row < DIM; row++){
                for(int col = 0; col < DIM; col++){
                    child->local_quadrupole[row][col] = node->local_quadrupole[row][col];
                }
            }
            
            //printf("node id = %llu child node L2L finished\n", node->id);
            ComputeLocalExpansions(child, root);
        }
    }
}

void naive_compute_gravity(pList* particles, double forces[][DIM]){
    for(int i = 0; i < N; i++){
        forces[i][0] = 0.0;
        forces[i][1] = 0.0;
        forces[i][2] = 0.0;

        for(int j = 0; j < N; j++){
            if(i == j) continue;

            double dx = particles->arr[j].pos[0] - particles->arr[i].pos[0];
            double dy = particles->arr[j].pos[1] - particles->arr[i].pos[1];
            double dz = particles->arr[j].pos[2] - particles->arr[i].pos[2];

            double dist_sq = dx*dx + dy*dy + dz*dz + EPSILON;
            double dist = sqrt(dist_sq);
            double force_mag = G * particles->arr[i].mass * particles->arr[j].mass / dist_sq;

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

    double forces[N][DIM];

    double init_center[] = {0.0, 0.0, 0.0};
    printf("Start BuildOctree\n");
    OctreeNode* root = BuildOctree(particles, init_center, HALF_BOX_SIZE, 0, 1);
    printf("Finish BuildOctree\n\n");

    printf("Start Upward Pass\n");
    ComputeMultipoles(root);
    printf("Finish Upward Pass\n\n");

    printf("Start Downward Pass\n");
    ComputeLocalExpansions(root, root);
    printf("Finish Downward Pass\n\n");

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
    size_t written = fwrite(forces, sizeof(double), DIM * N, fbin);
    if(written != DIM * N){
        printf("Failed to write all data!\n");
    }
    fclose(fbin);
    */
    return 0;
}
