#include <stdlib.h>

#include "main.h"

Coord* get_neighbors_far(Coord coord, Grid* grid, int accel) {
    int i;
    
    int num_cells = NCELL[accel-2];

    const Coord* DIST_CELLS;
    
    if (accel == 2) DIST_CELLS = CELLS2;
    else if (accel == 3) DIST_CELLS = CELLS3;
    else if (accel == 4) DIST_CELLS = CELLS4;

    Coord* neighbors = malloc(sizeof(Coord) * num_cells);
    int k = 0;
    int new_x, new_y;
    for (i = 0; i < 8; ++i) {
        new_x = coord.x + CELLS1[i].x;
        new_y = coord.y + CELLS1[i].y;
        if (new_x >= 0 && new_x < grid->ncol && new_y >= 0 && new_y < grid->nrow) {
            Coord neighbor = {new_x, new_y};
            neighbors[k++] = neighbor;
        }
    }

    for (i = 0; i < (num_cells - 8); ++i) {
        new_x = coord.x + DIST_CELLS[i].x;
        new_y = coord.y + DIST_CELLS[i].y;
        if (new_x >= 0 && new_x < grid->ncol && new_y >= 0 && new_y < grid->nrow && grid->terr[new_y * grid->ncol + new_x] == CORRIDOR) {
            Coord neighbor = {new_x, new_y};
            neighbors[k++] = neighbor;
        }
    }
    if (k < num_cells) neighbors[k].x = TURNOFF;
    return neighbors;
}

Coord* get_neighbors(Coord coord, Grid* grid) {
    int i;
    int num_cells = 8;

    Coord* neighbors = malloc(sizeof(Coord) * num_cells);

    int k = 0;
    int new_x, new_y;

    for (i = 0; i < num_cells; ++i) {

        new_x = coord.x + CELLS1[i].x;
        new_y = coord.y + CELLS1[i].y;

        if (new_x >= 0 && new_x < grid->ncol && new_y >= 0 && new_y < grid->nrow) {
            Coord neighbor = {new_x, new_y};
            neighbors[k++] = neighbor;
        }
    }

    if (k < num_cells) neighbors[k].x = TURNOFF;
    return neighbors;
}

double fission_ta(double n, double cta, double k) {
    if (cta >= n/k) return 0;
    else return n * (1.0 - (cta / (n/k)));
    //else return n - cta;
}

void fission(Model* model, Grid* grid) {
    int i, j;
    int n = model->agent_count;
    int num_cells = 8;
    int far_cells = NCELL[model->accel - 2];

    for (i = 0; i < n; ++i) {
        if (model->active[i] == 0) continue;

        Coord coord = model->agents[i];

        if (grid->terr[coord.y * grid->ncol + coord.x] == BARRIER) continue;

        double n = grid->population[coord.y * grid->ncol + coord.x];
        double local_k = grid->env[coord.y * grid->ncol + coord.x];
        double migrants = fission_ta(n, model->phi, local_k);

        if (migrants > 0) {
            Coord* nbr;
            Coord* free_nbr;
            int ncell = num_cells;
            if (grid->terr[coord.y * grid->ncol + coord.x] == CORRIDOR) {
                nbr = get_neighbors_far(coord, grid, model->accel);
                free_nbr = malloc(sizeof(Coord) * far_cells);
                ncell = far_cells;
            } else {
                nbr = get_neighbors(coord, grid);
                free_nbr = malloc(sizeof(Coord) * num_cells);
            }

            int len = 0;
            for (j = 0; j < ncell; ++j) {
                if (nbr[j].x == TURNOFF) {
                    break;
                }
                double nbr_k = grid->env[nbr[j].y * grid->ncol + nbr[j].x];
                if (grid->population[nbr[j].y * grid->ncol + nbr[j].x] < (nbr_k * model->phi) && grid->env[nbr[j].y * grid->ncol + nbr[j].x] > 0) {
                    free_nbr[len++] = nbr[j];
                }                    
            }

            if (len > 0) {
                grid->population[coord.y * grid->ncol + coord.x] -= migrants;

                int best_cell = 0;
                double best_env = 0;
                for (j = 0; j < len; ++j) {
                    double score = grid->env[free_nbr[j].y * grid->ncol + free_nbr[j].x];
                    if (score > best_env) {
                        best_env = score;
                        best_cell = j;
                    }
                }

                if (grid->arrival[free_nbr[best_cell].y * grid->ncol + free_nbr[best_cell].x] == 0) {
                    grid->arrival[free_nbr[best_cell].y * grid->ncol + free_nbr[best_cell].x] = model->tick;
                    model->active[model->agent_count] = 1;
                    model->agents[model->agent_count++] = free_nbr[best_cell];
                    model->agents[model->agent_count].x = TURNOFF;
                }
                grid->population[free_nbr[best_cell].y * grid->ncol + free_nbr[best_cell].x] += migrants;

            } else {
                model->active[i] = 0;
            }

            free(nbr);
            free(free_nbr);

        }
    }
}