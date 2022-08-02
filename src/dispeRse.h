#ifndef DISPERSE_H
#define DISPERSE_H

#define TURNOFF -1
#define BARRIER 1
#define CORRIDOR 2

/**
 * @brief Structure to represent xy coordinates.
 */
typedef struct Coord {
    int x;
    int y;
} Coord;

/**
 * @brief Structure to represent a spatial grid with various layers
 * (environment, terrain, etc.).
 */
typedef struct Grid {
    int nrow;
    int ncol;
    double* population;
    double* environment;
    int* terrain;
    int* arrival;
} Grid;

/**
 * @brief Structure to represent a dispersal model. It contains the demographic
 * and mobility parameters, also keeping track of active cells and time.
 */
typedef struct Model {
    double r;
    double phi;
    double t;
    int accel;
    double gamma;
    Coord* agents;
    int* active;
    int agent_count;
    int tick;
} Model;

/* For the sake of efficiency, the XY distances to the cells of a neighborhood
of up to 4 cells is hard-coded below. */

// Distance of 1 (Moore neighborhood)
static const Coord CELLS1[8] = {
    {-1,-1}, {-1, 0}, {-1, 1}, { 0,-1},
    { 0, 1}, { 1,-1}, { 1, 0}, { 1, 1}
};

// Distance of 2. Note that cells of the Moore neighborhood are excluded.
static const Coord CELLS2[12] = {
    {-2, 1}, {-2, 0}, {-2,-1}, {-1, 2},
    {-1,-2}, { 0, 2}, { 0,-2}, { 1, 2},
    { 1,-2}, { 2, 1}, { 2, 0}, { 2,-1}
};

// Distance of 3. Note that cells of the Moore neighborhood are excluded.
static const Coord CELLS3[28] = {
    {-3,-1}, {-3, 0}, {-3, 1}, {-2,-2},
    {-2,-1}, {-2, 0}, {-2, 1}, {-2, 2},
    {-1,-3}, {-1,-2}, {-1, 2}, {-1, 3},
    { 0,-3}, { 0,-2}, { 0, 2}, { 0, 3},
    { 1,-3}, { 1,-2}, { 1, 2}, { 1, 3},
    { 2,-2}, { 2,-1}, { 2, 0}, { 2, 1},
    { 2, 2}, { 3,-1}, { 3, 0}, { 3, 1}
};

// Distance of 4. Note that cells of the Moore neighborhood are excluded.
static const Coord CELLS4[60] = {
    {-4,-2}, {-4,-1}, {-4, 0}, {-4, 1},
    {-4, 2}, {-3,-3}, {-3,-2}, {-3,-1},
    {-3, 0}, {-3, 1}, {-3, 2}, {-3, 3},
    {-2,-4}, {-2,-3}, {-2,-2}, {-2,-1},
    {-2, 0}, {-2, 1}, {-2, 2}, {-2, 3},
    {-2, 4}, {-1,-4}, {-1,-3}, {-1,-2},
    {-1, 2}, {-1, 3}, {-1, 4}, { 0,-4},
    { 0,-3}, { 0,-2}, { 0, 2}, { 0, 3},
    { 0, 4}, { 1,-4}, { 1,-3}, { 1,-2},
    { 1, 2}, { 1, 3}, { 1, 4}, { 2,-4},
    { 2,-3}, { 2,-2}, { 2,-1}, { 2, 0},
    { 2, 1}, { 2, 2}, { 2, 3}, { 2, 4},
    { 3,-3}, { 3,-2}, { 3,-1}, { 3, 0},
    { 3, 1}, { 3, 2}, { 3, 3}, { 4,-2},
    { 4,-1}, { 4, 0}, { 4, 1}, { 4, 2}
};

// Array with max number of free cells within a distance of 2, 3 and 4.
static const int NCELL[3] = {20,36,68};

// Function declarations.
Coord* get_neighbors(Coord coord, Grid *grid);
double log_growth(double n, double r, double k, double t);
double fission_ta(double n, double cta, double k);
void grow(Model *model, Grid *grid);
void fission(Model *model, Grid *grid);
void run_model(int *nrow, int *ncol, double *environment, int *terrain,
               double *population, int *arrival, int *x, int *y, int *start,
               int *num_origins, int *num_iter, double *r, double *phi,
               double *t, int *accel, double *gamma, int *updates);

#endif
