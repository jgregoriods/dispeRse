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
    double dist;
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
of up to 4 cells is hard-coded below. The inverse of the square distance is included.*/

// Distance of 1 (Moore neighborhood)
static const Coord CELLS1[8] = {
    {-1,-1, 0.5}, {-1, 0, 1.0}, {-1, 1, 0.5}, { 0,-1, 1.0},
    { 0, 1, 1.0}, { 1,-1, 0.5}, { 1, 0, 1.0}, { 1, 1, 0.5}
};

// Distance of 2. Note that cells of the Moore neighborhood are excluded.
static const Coord CELLS2[12] = {
    {-2, 1, 0.2}, {-2, 0, 0.25}, {-2,-1, 0.2 }, {-1, 2, 0.2},
    {-1,-2, 0.2}, { 0, 2, 0.25}, { 0,-2, 0.25}, { 1, 2, 0.2},
    { 1,-2, 0.2}, { 2, 1, 0.2 }, { 2, 0, 0.25}, { 2,-1, 0.2}
};

// Distance of 3. Note that cells of the Moore neighborhood are excluded.
static const Coord CELLS3[28] = {
    {-3,-1, 0.1 }, {-3, 0, 0.11}, {-3, 1, 0.1 }, {-2,-2, 0.12},
    {-2,-1, 0.2 }, {-2, 0, 0.25}, {-2, 1, 0.2 }, {-2, 2, 0.12},
    {-1,-3, 0.1 }, {-1,-2, 0.2 }, {-1, 2, 0.2 }, {-1, 3, 0.1 },
    { 0,-3, 0.11}, { 0,-2, 0.25}, { 0, 2, 0.25}, { 0, 3, 0.11},
    { 1,-3, 0.1 }, { 1,-2, 0.2 }, { 1, 2, 0.2 }, { 1, 3, 0.1 },
    { 2,-2, 0.12}, { 2,-1, 0.2 }, { 2, 0, 0.25}, { 2, 1, 0.2 },
    { 2, 2, 0.12}, { 3,-1, 0.1 }, { 3, 0, 0.11}, { 3, 1, 0.1 }
};

// Distance of 4. Note that cells of the Moore neighborhood are excluded.
static const Coord CELLS4[60] = {
    {-4,-2, 0.05}, {-4,-1, 0.06}, {-4, 0, 0.06}, {-4, 1, 0.06},
    {-4, 2, 0.05}, {-3,-3, 0.05}, {-3,-2, 0.08}, {-3,-1, 0.1 },
    {-3, 0, 0.11}, {-3, 1, 0.1 }, {-3, 2, 0.08}, {-3, 3, 0.05},
    {-2,-4, 0.05}, {-2,-3, 0.08}, {-2,-2, 0.12}, {-2,-1, 0.2 },
    {-2, 0, 0.25}, {-2, 1, 0.2 }, {-2, 2, 0.12}, {-2, 3, 0.08},
    {-2, 4, 0.05}, {-1,-4, 0.06}, {-1,-3, 0.1 }, {-1,-2, 0.2 },
    {-1, 2, 0.2 }, {-1, 3, 0.1 }, {-1, 4, 0.06}, { 0,-4, 0.06},
    { 0,-3, 0.11}, { 0,-2, 0.25}, { 0, 2, 0.25}, { 0, 3, 0.11},
    { 0, 4, 0.06}, { 1,-4, 0.06}, { 1,-3, 0.1 }, { 1,-2, 0.2 },
    { 1, 2, 0.2 }, { 1, 3, 0.1 }, { 1, 4, 0.06}, { 2,-4, 0.05},
    { 2,-3, 0.08}, { 2,-2, 0.12}, { 2,-1, 0.2 }, { 2, 0, 0.25},
    { 2, 1, 0.2 }, { 2, 2, 0.12}, { 2, 3, 0.08}, { 2, 4, 0.05},
    { 3,-3, 0.05}, { 3,-2, 0.08}, { 3,-1, 0.1 }, { 3, 0, 0.11},
    { 3, 1, 0.1 }, { 3, 2, 0.08}, { 3, 3, 0.05}, { 4,-2, 0.05},
    { 4,-1, 0.06}, { 4, 0, 0.06}, { 4, 1, 0.06}, { 4, 2, 0.05}
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
