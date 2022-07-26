#ifndef MAIN_H
#define MAIN_H

#define TURNOFF -1
#define BARRIER 1
#define CORRIDOR 2

typedef struct Coord {
    int x;
    int y;
} Coord;

typedef struct Grid {
    double* population;
    double* env;
    int* terr;
    int* arrival;
    int nrow;
    int ncol;
} Grid;

typedef struct Model {
    double r;
    double phi;
    double t;
    Coord* agents;
    int* active;
    int agent_count;
    int tick;
    double gamma;
} Model;

static const Coord CELLS1[8] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};
static const Coord CELLS2[20] = {
    {-2, 1}, {-2, 0}, {-2, -1},
    {-1, 2}, {-1, 1}, {-1, 0}, {-1, -1}, {-1, -2},
    {0, 2}, {0, 1}, {0, -1}, {0, -2},
    {1, 2}, {1, 1}, {1, 0}, {1, -1}, {1, -2},
    {2, 1}, {2, 0}, {2, -1}
};
static const Coord CELLS3[28] = {
    {-3, -1}, {-3, 0}, {-3, 1}, {-2, -2}, {-2, -1}, {-2, 0},
    {-2, 1}, {-2, 2}, {-1, -3}, {-1, -2},
    {-1, 2}, {-1, 3}, {0, -3}, {0, -2},
    {0, 2}, {0, 3}, {1, -3}, {1, -2},
    {1, 2}, {1, 3}, {2, -2}, {2, -1}, {2, 0}, {2, 1},
    {2, 2}, {3, -1}, {3, 0}, {3, 1}
};

Coord* get_neighbors(Coord coord, Grid* grid);
double log_growth(double n, double r, double k, double t);
double fission_ta(double n, double cta, double k);
void grow(Model* model, Grid* grid);
void fission(Model* model, Grid* grid);
void run_model(int* nrow, int* ncol, double* population, double* env, int* arrival,
               double* r, double* phi, int* start,
               int* x, int* y, int* iter, int* num_origins, double* t, int* terrain);

#endif