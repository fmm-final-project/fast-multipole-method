#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PARTICLES_PER_LEAF 10
#define MAX_LEVELS 10
#define MAX_CELLS 100000
#define MAX_CHILDREN 8

// Structure for a particle
typedef struct {
    double x[3];
    double mass;
    double potential;
    double force[3];
} Particle;

// Octree cell structure
typedef struct Cell {
    double center[3];
    double half_width;
    int *particles;
    int n_particles;
    struct Cell *children[MAX_CHILDREN];
    struct Cell *parent;
    int child_index;      // index in parent's children array
    int level;
    int is_leaf;
    // Multipole expansion up to quadrupole
    double M0;
    double M1[3];
    double M2[3][3];
    // Local expansion up to quadrupole
    double L0;
    double L1[3];
    double L2[3][3];
} Cell;

// Global storage of cells by level
Cell *cells_by_level[MAX_LEVELS][MAX_CELLS];
int num_cells_by_level[MAX_LEVELS] = {0};
int max_level = 0;

// Create a new cell
Cell* create_cell(double cx, double cy, double cz, double hw, Cell *parent, int cidx, int level) {
    Cell *c = (Cell*)malloc(sizeof(Cell));
    c->center[0]=cx; c->center[1]=cy; c->center[2]=cz;
    c->half_width = hw;
    c->parent = parent;
    c->child_index = cidx;
    c->level = level;
    c->n_particles = 0;
    c->particles = NULL;
    c->is_leaf = 0;
    // init children
    for(int i=0;i<MAX_CHILDREN;i++) c->children[i]=NULL;
    // add to global list
    if (level < MAX_LEVELS) {
        cells_by_level[level][ num_cells_by_level[level]++ ] = c;
        if(level > max_level) max_level = level;
    }
    return c;
}

// Build the octree recursively
Cell* build_tree_rec(Particle *particles, int *inds, int n, double cx, double cy, double cz,
                     double hw, Cell *parent, int cidx, int level) {
    Cell *cell = create_cell(cx, cy, cz, hw, parent, cidx, level);
    cell->n_particles = n;
    cell->particles = (int*)malloc(n * sizeof(int));
    memcpy(cell->particles, inds, n*sizeof(int));
    if(n <= MAX_PARTICLES_PER_LEAF) {
        cell->is_leaf = 1;
    } else {
        cell->is_leaf = 0;
        // partition into 8 octants
        int *tmp_inds[MAX_CHILDREN]; int counts[MAX_CHILDREN] = {0};
        for(int i=0;i<MAX_CHILDREN;i++) tmp_inds[i] = (int*)malloc(n * sizeof(int));
        for(int i=0;i<n;i++){
            int idx = inds[i];
            double dx = particles[idx].x[0] - cx;
            double dy = particles[idx].x[1] - cy;
            double dz = particles[idx].x[2] - cz;
            int oct = (dx>0) | ((dy>0)<<1) | ((dz>0)<<2);
            tmp_inds[oct][ counts[oct]++ ] = idx;
        }
        for(int o=0;o<MAX_CHILDREN;o++){
            if(counts[o]>0) {
                double shift[3] = { ((o&1)?hw/2:-hw/2), ((o&2)?hw/2:-hw/2), ((o&4)?hw/2:-hw/2) };
                cell->children[o] = build_tree_rec(particles, tmp_inds[o], counts[o],
                        cx + shift[0], cy + shift[1], cz + shift[2], hw/2, cell, o, level+1);
            }
        }
        for(int i=0;i<MAX_CHILDREN;i++) free(tmp_inds[i]);
    }
    return cell;
}

// Particle-to-Multipole for leaves
void P2M(Cell *cell, Particle *particles) {
    // zero
    cell->M0 = 0;
    for(int i=0;i<3;i++) cell->M1[i]=0;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) cell->M2[i][j]=0;
    // accumulate
    for(int pi=0; pi<cell->n_particles; pi++){
        Particle *p = &particles[ cell->particles[pi] ];
        double dx = p->x[0] - cell->center[0];
        double dy = p->x[1] - cell->center[1];
        double dz = p->x[2] - cell->center[2];
        cell->M0 += p->mass;
        cell->M1[0] += p->mass * dx;
        cell->M1[1] += p->mass * dy;
        cell->M1[2] += p->mass * dz;
        cell->M2[0][0] += p->mass * dx*dx;
        cell->M2[0][1] += p->mass * dx*dy;
        cell->M2[0][2] += p->mass * dx*dz;
        cell->M2[1][0] += p->mass * dy*dx;
        cell->M2[1][1] += p->mass * dy*dy;
        cell->M2[1][2] += p->mass * dy*dz;
        cell->M2[2][0] += p->mass * dz*dx;
        cell->M2[2][1] += p->mass * dz*dy;
        cell->M2[2][2] += p->mass * dz*dz;
    }
}

// Upward pass: P2M + M2M
void upward(Cell *cell, Particle *particles) {
    if(cell->is_leaf) {
        P2M(cell, particles);
    } else {
        // zero parent multipoles
        cell->M0 = 0;
        for(int i=0;i<3;i++) cell->M1[i]=0;
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) cell->M2[i][j]=0;
        // collect from children
        for(int c=0;c<MAX_CHILDREN;c++){
            Cell *ch = cell->children[c];
            if(!ch) continue;
            upward(ch, particles);
            double d[3] = { ch->center[0]-cell->center[0], ch->center[1]-cell->center[1], ch->center[2]-cell->center[2] };
            // M2M
            cell->M0 += ch->M0;
            for(int i=0;i<3;i++) cell->M1[i] += ch->M1[i] + ch->M0 * d[i];
            for(int i=0;i<3;i++) for(int j=0;j<3;j++)
                cell->M2[i][j] += ch->M2[i][j]
                                + ch->M1[i]*d[j] + d[i]*ch->M1[j]
                                + ch->M0 * d[i]*d[j];
        }
    }
}

// Multipole-to-Local translation (M2L) up to quadrupole
void M2L(Cell *target, Cell *source) {
    // vector between centers
    double d[3] = { target->center[0]-source->center[0],
                    target->center[1]-source->center[1],
                    target->center[2]-source->center[2] };
    double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    double r = sqrt(r2);
    double inv_r = 1.0/r;
    double inv_r3 = inv_r/r2;
    double inv_r5 = inv_r3/r2;
    double inv_r7 = inv_r5/r2;
    double inv_r9 = inv_r7/r2;
    // precompute
    double tr2 = source->M2[0][0] + source->M2[1][1] + source->M2[2][2];
    double M1dotd = source->M1[0]*d[0] + source->M1[1]*d[1] + source->M1[2]*d[2];
    double dTM2d = 0;
    double M2d[3] = {0,0,0};
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            dTM2d += source->M2[i][j]*d[i]*d[j];
            M2d[i] += source->M2[i][j]*d[j];
        }
    }
    // zero increments
    double L0 = 0;
    double L1[3] = {0,0,0};
    double L2[3][3] = {{0}};
    // monopole
    L0 += source->M0 * inv_r;
    for(int i=0;i<3;i++)
        L1[i] += -source->M0 * d[i] * inv_r3;
    // dipole
    for(int i=0;i<3;i++) {
        L0 += -source->M1[i] * d[i] * inv_r3;
        L1[i] += -source->M1[i] * inv_r3 + 3 * M1dotd * d[i] * inv_r5;
    }
    // quadrupole onto L0, L1
    L0 += 0.5 * (-tr2 * inv_r3 + 3 * dTM2d * inv_r5);
    for(int k=0;k<3;k++) {
        L1[k] += 0.5 * ( 3*tr2 * d[k] * inv_r5 + 6*M2d[k] * inv_r5 - 15*dTM2d * d[k] * inv_r7 );
    }
    // quadrupole onto L2
    for(int k=0;k<3;k++) for(int l=0;l<3;l++){
        // M0 term
        L2[k][l] += source->M0 * ( - (k==l?inv_r3:0) + 3 * d[k]*d[l] * inv_r5 );
        // M1 term
        for(int i=0;i<3;i++){
            double term = 3*( (k==l?d[i]:0) + (k==i?d[l]:0) + (l==i?d[k]:0) ) * inv_r5
                          - 15 * d[i]*d[k]*d[l] * inv_r7;
            L2[k][l] += source->M1[i] * term;
        }
        // M2 term
        double sum = 0;
        for(int i=0;i<3;i++) for(int j=0;j<3;j++){
            double term = 3*(-((i==j&&(k==l)) + (i==k&&(j==l)) + (i==l&&(j==k))))*inv_r5
                          + 15*( (i==j?d[k]*d[l]:0)
                                + (i==k?d[j]*d[l]:0)
                                + (i==l?d[j]*d[k]:0)
                                + (j==k?d[i]*d[l]:0)
                                + (j==l?d[i]*d[k]:0)
                                + (k==l?d[i]*d[j]:0) ) * inv_r7
                          - 105 * d[i]*d[j]*d[k]*d[l] * inv_r9;
            sum += source->M2[i][j] * term;
        }
        L2[k][l] += 0.5 * sum;
    }
    // accumulate into target local
    target->L0 += L0;
    for(int i=0;i<3;i++) target->L1[i] += L1[i];
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) target->L2[i][j] += L2[i][j];
}

// Local-to-Local translation (L2L)
void L2L(Cell *child, Cell *parent) {
    double d[3] = { child->center[0]-parent->center[0],
                    child->center[1]-parent->center[1],
                    child->center[2]-parent->center[2] };
    // L0
    double newL0 = parent->L0;
    for(int i=0;i<3;i++) newL0 += parent->L1[i]*d[i];
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) newL0 += 0.5*parent->L2[i][j]*d[i]*d[j];
    child->L0 = newL0;
    // L1
    for(int j=0;j<3;j++){
        double sum = parent->L1[j];
        for(int i=0;i<3;i++) sum += parent->L2[j][i] * d[i];
        child->L1[j] = sum;
    }
    // L2
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) child->L2[i][j] = parent->L2[i][j];
}

// Evaluate Local expansion at particles in a leaf (L2P)
void L2P(Cell *cell, Particle *particles) {
    for(int pi=0; pi<cell->n_particles; pi++){
        Particle *p = &particles[ cell->particles[pi] ];
        double rvec[3] = { p->x[0]-cell->center[0], p->x[1]-cell->center[1], p->x[2]-cell->center[2] };
        // potential
        double pot = cell->L0
                   + cell->L1[0]*rvec[0] + cell->L1[1]*rvec[1] + cell->L1[2]*rvec[2]
                   + 0.5*( cell->L2[0][0]*rvec[0]*rvec[0]
                          + (cell->L2[0][1]+cell->L2[1][0])*rvec[0]*rvec[1]
                          + (cell->L2[0][2]+cell->L2[2][0])*rvec[0]*rvec[2]
                          + cell->L2[1][1]*rvec[1]*rvec[1]
                          + (cell->L2[1][2]+cell->L2[2][1])*rvec[1]*rvec[2]
                          + cell->L2[2][2]*rvec[2]*rvec[2] );
        p->potential += pot;
        // force = -m * grad
        for(int k=0;k<3;k++){
            double fval = cell->L1[k];
            for(int i=0;i<3;i++) fval += cell->L2[i][k] * rvec[i];
            p->force[k] -= p->mass * fval;
        }
    }
}

// Direct particle-to-particle for near-field (P2P) among siblings
void direct_interaction(Particle *p, Particle *q) {
    double dx = p->x[0] - q->x[0];
    double dy = p->x[1] - q->x[1];
    double dz = p->x[2] - q->x[2];
    double r2 = dx*dx + dy*dy + dz*dz + 1e-12;
    double inv_r3 = 1.0/(r2 * sqrt(r2));
    double f = p->mass * q->mass * inv_r3;
    p->force[0] -= f * dx;
    p->force[1] -= f * dy;
    p->force[2] -= f * dz;
    q->force[0] += f * dx;
    q->force[1] += f * dy;
    q->force[2] += f * dz;
}

void P2P(Cell *cell, Particle *particles) {
    Cell *parent = cell->parent;
    if(!parent) return;
    int self_idx = cell->child_index;
    for(int c=self_idx; c<MAX_CHILDREN; c++){
        Cell *nbr = parent->children[c];
        if(!nbr || !nbr->is_leaf) continue;
        if(nbr == cell) {
            // same cell: i < j
            for(int a=0;a<cell->n_particles;a++){
                Particle *pi = &particles[ cell->particles[a] ];
                for(int b=a+1;b<cell->n_particles;b++){
                    Particle *pj = &particles[ cell->particles[b] ];
                    direct_interaction(pi,pj);
                }
            }
        } else {
            // neighbor leaf
            for(int a=0;a<cell->n_particles;a++){
                Particle *pi = &particles[ cell->particles[a] ];
                for(int b=0;b<nbr->n_particles;b++){
                    Particle *pj = &particles[ nbr->particles[b] ];
                    direct_interaction(pi,pj);
                }
            }
        }
    }
}

// Downward pass: L2L + L2P + P2P
void downward(Cell *cell, Particle *particles) {
    for(int c=0;c<MAX_CHILDREN;c++){
        Cell *ch = cell->children[c];
        if(!ch) continue;
        // init child local = parent local
        ch->L0 = cell->L0;
        for(int i=0;i<3;i++) ch->L1[i] = cell->L1[i];
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) ch->L2[i][j] = cell->L2[i][j];
        // translate
        L2L(ch, cell);
        downward(ch, particles);
    }
    if(cell->is_leaf) {
        L2P(cell, particles);
        P2P(cell, particles);
    }
}

// Main FMM driver
void FMM(Cell *root, Particle *particles) {
    upward(root, particles);
    // init root local
    root->L0 = 0;
    for(int i=0;i<3;i++) root->L1[i]=0;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) root->L2[i][j]=0;
    // M2L for well-separated pairs at each level
    for(int lvl=1; lvl<=max_level; lvl++){
        int ncell = num_cells_by_level[lvl];
        for(int i=0;i<ncell;i++){
            Cell *ti = cells_by_level[lvl][i];
            for(int j=0;j<ncell;j++){
                if(i==j) continue;
                Cell *sj = cells_by_level[lvl][j];
                double dx = ti->center[0]-sj->center[0];
                double dy = ti->center[1]-sj->center[1];
                double dz = ti->center[2]-sj->center[2];
                double dist = sqrt(dx*dx+dy*dy+dz*dz);
                if(dist > 2*ti->half_width) {
                    M2L(ti, sj);
                }
            }
        }
    }
    // downward translates + evaluation
    downward(root, particles);
}

// Example usage: read particles from stdin: N, then lines of x y z mass
int main(){

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

    int N = num_of_doubles / 7;

    double* buffer = (double*)malloc(file_size);
    if(!buffer){
        printf("Memory allocate failed\n");
        fclose(fptr);
        exit(1);
    }

    if(fread(buffer, sizeof(double), num_of_doubles, fptr)){};

    Particle* particles = (Particle*)malloc(N * sizeof(Particle));
    for(int i = 0; i < N; i++){
        particles[i].mass = buffer[i * 7];
        particles[i].x[0] = buffer[i * 7 + 1];
        particles[i].x[1] = buffer[i * 7 + 2];
        particles[i].x[2] = buffer[i * 7 + 3];
        particles[i].potential = 0;
        for(int k=0;k<3;k++) particles[i].force[k]=0;
    }

    int *inds = (int*)malloc(N*sizeof(int));
    for(int i=0;i<N;i++) inds[i]=i;
    double cx=0.0, cy=0.0, cz=0.0;
    double hw=100.0; // domain half-width; adjust to cover all particles
    Cell *root = build_tree_rec(particles, inds, N, cx, cy, cz, hw, NULL, 0, 0);
    FMM(root, particles);
    // output forces
    for(int i=0;i<N;i++){
        printf("%.10e %.10e %.10e\n", particles[i].force[0], particles[i].force[1], particles[i].force[2]);
    }
    return 0;
}