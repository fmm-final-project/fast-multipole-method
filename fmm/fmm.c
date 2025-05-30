#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define G 1.0                      // Gravitational constant
#define EPSILON 0.0                // Softening factor to prevent division by zero

#define INIT_LEN 16
#define MAX_PARTICLES_PER_LEAF 10
#define MAX_DEPTH 20
#define HALF_BOX_SIZE 100.0
#define DIM 3
#define THETA 0.3

typedef struct{
    double mass;
    double pos[DIM];
} Particle;

typedef struct{
    int* arr;
    int size;
    int capacity;
} pList;

typedef struct octreenode{
    unsigned long long id;
    double center[DIM];
    double half_size;
    pList* pid_list;
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

int N;                 // Number of particles
Particle* particles;   // All particles
double** forces;       // All forces
oList* leafNodes;      // All leaf octree nodes

// Basic operations of particle lists
pList* pListInit(){
    pList* newList = (pList*)malloc(sizeof(pList));
    newList->capacity = INIT_LEN;
    newList->size = 0;
    newList->arr = (int*)malloc(INIT_LEN * sizeof(int));
    return newList;
}

void pListInsert(pList* list, int pid){
    list->arr[list->size] = pid;
    list->size++;
    if(list->size == list->capacity){
        list->capacity *= 2;
        list->arr = (int*)realloc(list->arr, list->capacity * sizeof(int));
    }
}

void pListFree(pList* list){
    free(list->arr);
    free(list);
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

void oListFree(oList* list){
    free(list->arr);
    free(list);
}

void FreeOctreeNode(OctreeNode* node){
    pListFree(node->pid_list);
    free(node->children_list);
    free(node);
}

void FreeOctree(OctreeNode* node){
    if(node->children_list != NULL){
        for(int i = 0; i < 8; i++){
            if(node->children_list[i] != NULL) FreeOctree(node->children_list[i]);
        }
    }
    FreeOctreeNode(node);
}

// Compute which octree child node the particle should go
int ComputeOctreeIndex(int pid, double* center){
    int x = particles[pid].pos[0] > center[0] ? 1 : 0;
    int y = particles[pid].pos[1] > center[1] ? 1 : 0;
    int z = particles[pid].pos[2] > center[2] ? 1 : 0;
    return x + 2 * y + 4 * z;
}

// Compute the center of a child octree node from its parent's center
void ComputeChildCenter(double* center, double half_size, int index, double child_center[]){
    double offset = half_size / 2;
    int x = (index & 1) == 1 ? 1 : -1;
    int y = ((index / 2) & 1) == 1 ? 1 : -1;
    int z = ((index / 4) & 1) == 1 ? 1 : -1;
    child_center[0] = center[0] + offset * x;
    child_center[1] = center[1] + offset * y;
    child_center[2] = center[2] + offset * z;
}

// Building the octree of given particles
OctreeNode* BuildOctree(pList* pid_list, double center[], double half_size, int depth, unsigned long long id){
    OctreeNode* node = (OctreeNode*)calloc(1, sizeof(OctreeNode));
    for(int i = 0; i < DIM; i++){
        node->center[i] = center[i];
    }
    node->half_size = half_size;
    node->pid_list = pid_list;
    node->children_list = NULL;
    node->id = id;

    if(pid_list->size < MAX_PARTICLES_PER_LEAF || depth > MAX_DEPTH){
        // Is a leaf node
        printf("Construct a leaf node id = %llu with %d particles\n", node->id, pid_list->size);
        // Insert into LeafNodes
        oListInsert(leafNodes, node);
        return node;
    }

    // Is an internal node, split children
    node->children_list = (OctreeNode**)malloc(8 * sizeof(OctreeNode*));
    pList** children_pid_lists = (pList**)malloc(8 * sizeof(pList*));
    for(int i = 0; i < 8; i++){
        children_pid_lists[i] = pListInit();
    }

    for(int i = 0; i < pid_list->size; i++){
        int octreeIndex = ComputeOctreeIndex(pid_list->arr[i], center);
        pListInsert(children_pid_lists[octreeIndex], pid_list->arr[i]);
    }

    for(int i = 0; i < 8; i++){
        if(children_pid_lists[i]->size > 0){
            double child_center[DIM];
            ComputeChildCenter(center, half_size, i, child_center);
            double child_half_size = half_size / 2.0;
            node->children_list[i] = BuildOctree(children_pid_lists[i], child_center, child_half_size, depth + 1, id * 8 + i);
        }
        else{
            pListFree(children_pid_lists[i]);
            node->children_list[i] = NULL;
        }
    }
    free(children_pid_lists);
    printf("Construct an internal node id = %llu with %d nodes\n", node->id, pid_list->size);
    return node;
}

// Compute the multipole expansions for each octree node
void ComputeMultipoles(OctreeNode* node){
    //printf("visit node with center x = %.10e\n", node->center[0]);

    if(node->children_list == NULL){
        for(int i = 0; i < node->pid_list->size; i++){
            Particle particle = particles[node->pid_list->arr[i]];
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
                                                + child->dipole[row] * dx[col]
                                                + m * dx[row] * dx[col];
                }
            }

        }
        printf("internal node id = %llu multipole finished, monopole = %.10e\n", node->id, node->monopole);
    }
}

bool BoxesOverlap(OctreeNode* a, OctreeNode* b){
    double ax_min = a->center[0] - a->half_size;
    double ax_max = a->center[0] + a->half_size;
    double ay_min = a->center[1] - a->half_size;
    double ay_max = a->center[1] + a->half_size;
    double az_min = a->center[2] - a->half_size;
    double az_max = a->center[2] + a->half_size;

    double bx_min = b->center[0] - b->half_size;
    double bx_max = b->center[0] + b->half_size;
    double by_min = b->center[1] - b->half_size;
    double by_max = b->center[1] + b->half_size;
    double bz_min = b->center[2] - b->half_size;
    double bz_max = b->center[2] + b->half_size;

    return !(ax_max < bx_min || ax_min > bx_max ||
             ay_max < by_min || ay_min > by_max ||
             az_max < bz_min || az_min > bz_max);
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
    if(src->children_list != NULL && size2 / dist2 > THETA * THETA){
        // Not well-seperated
        for(int i = 0; i < 8; i++){
            GetInteractionList(target, src->children_list[i], list);
        }
    }
    else{
        if(!BoxesOverlap(target, src)){
            // Well-seperated
            oListInsert(list, src);
        }
    }
}

void ComputeLocalExpansions(OctreeNode* node, OctreeNode* root){
    oList* interaction_list = oListInit();
    GetInteractionList(node, root, interaction_list);

    // M2L translation
    for(int i = 0; i < interaction_list->size; i++){
        OctreeNode* src = interaction_list->arr[i];

        double dx[DIM];
        double r2 = EPSILON;
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
    oListFree(interaction_list);
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

void ComputeDirectForce(int p1, int p2){
    if(p1 == p2) return;
    double dx[DIM];
    double dist2 = EPSILON;
    for(int i = 0; i < DIM; i++){
        dx[i] = particles[p2].pos[i] - particles[p1].pos[i];
        dist2 += dx[i] * dx[i];
    }
    double dist = sqrt(dist2);
    double force_mag = G * particles[p1].mass * particles[p2].mass / dist2;

    // Directional force components
    for(int i = 0; i < DIM; i++){
        forces[p1][i] += force_mag * dx[i] / dist;
    }
}

void ComputeForceFromLocalExpansion(int pid, OctreeNode* node){
    double dx[DIM];
    for(int i = 0; i < DIM; i++){
        dx[i] = particles[pid].pos[i] - node->center[i];
    }
    double r2 = EPSILON;
    for(int j = 0; j < DIM; j++) {
        r2 += dx[j] * dx[j];
    }
    double r = sqrt(r2);
    double rinv3 = 1.0 / (r2 * r);

    for(int i = 0; i < DIM; i++){
        
        double grad_phi = node->local_dipole[i];
        for(int j = 0; j < DIM; j++){
            grad_phi += node->local_quadrupole[i][j] * dx[j];
        }
        grad_phi += node->local_monopole * dx[i] * rinv3;
        forces[pid][i] -= G * particles[pid].mass * grad_phi;
    }
}

void EvaluateForces(OctreeNode* node){
    for(int i = 0; i < node->pid_list->size; i++){
        int p1 = node->pid_list->arr[i];

        // Near-field
        for(int j = 0; j < leafNodes->size; j++){
            OctreeNode* neighbor = leafNodes->arr[j];
            if(!BoxesOverlap(node, neighbor)) continue;

            for(int k = 0; k < neighbor->pid_list->size; k++){
                int p2 = neighbor->pid_list->arr[k];
                ComputeDirectForce(p1, p2);
            }
        }

        // Far-field
        ComputeForceFromLocalExpansion(p1, node);
    }
}

int main() {
    // Input particle data
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

    particles = (Particle*)malloc(N * sizeof(Particle));
    pList* pid_list = (pList*)malloc(sizeof(pList));
    pid_list->arr = (int*)malloc(N * sizeof(int));
    pid_list->capacity = N;
    pid_list->size = N;
    for(int i = 0; i < N; i++){
        pid_list->arr[i] = i;
        particles[i].mass = buffer[i * 7];
        particles[i].pos[0] = buffer[i * 7 + 1];
        particles[i].pos[1] = buffer[i * 7 + 2];
        particles[i].pos[2] = buffer[i * 7 + 3];
    }
    forces = (double**)malloc(N * sizeof(double*));
    forces[0] = (double*)calloc(DIM * N, sizeof(double));
    for(int i = 1; i < N; i++){
        forces[i] = forces[0] + DIM * i;
    }
    fclose(fptr);
    free(buffer);

    // Compute gravitational forces

    // Build Octree
    double init_center[] = {0.0, 0.0, 0.0};
    leafNodes = oListInit();
    printf("Start BuildOctree\n");
    OctreeNode* root = BuildOctree(pid_list, init_center, HALF_BOX_SIZE, 0, 1);
    printf("Finish BuildOctree\n\n");

    // Upward Pass
    printf("Start Upward Pass\n");
    ComputeMultipoles(root);
    printf("Finish Upward Pass\n\n");

    // Downward Pass
    printf("Start Downward Pass\n");
    ComputeLocalExpansions(root, root);
    printf("Finish Downward Pass\n\n");

    // Evaluate forces for particles inside every leaf nodes
    printf("Start Evaluation\n");
    for(int i = 0; i < leafNodes->size; i++){
        EvaluateForces(leafNodes->arr[i]);
    }
    printf("Finish Evaluation\n");
    oListFree(leafNodes);
    FreeOctree(root);

    // Output the force vectors
    FILE* fcsv;
    fcsv = fopen("force_fmm.csv", "w");
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
    size_t written = fwrite(forces[0], sizeof(double), DIM * N, fbin);
    if(written != DIM * N){
        printf("Failed to write all data!\n");
    }
    fclose(fbin);
    
    free(particles);
    free(forces[0]);
    free(forces);
    return 0;
}
