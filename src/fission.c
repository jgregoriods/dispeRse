#include <stdlib.h>
#include <math.h>

#include "dispeRse.h"

/**
 * @brief Returns cells with terrain values of corridor wihin a distance (in
 * cells) defined by the parameter accel.
 * 
 * @param coord The coordinates of the central cell.
 * @param grid The model grid structure.
 * @param accel An integer in {2,3,4} that determines the max distance of
 * the search.
 * @return An array of coord structures.
 */
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
    double dist;

    // consider all cells in the Moore neighborhood
    for (i = 0; i < 8; i++) {
        new_x = coord.x + CELLS1[i].x;
        new_y = coord.y + CELLS1[i].y;
        dist = CELLS1[i].dist;

        if (new_x >= 0 && new_x < grid->ncol && new_y >= 0 && new_y < grid->nrow) {
            Coord neighbor = {new_x, new_y, dist};
            neighbors[k++] = neighbor;
        }
    }

    // add distant cells if they are in a corridor
    for (i = 0; i < (num_cells - 8); i++) {
        new_x = coord.x + DIST_CELLS[i].x;
        new_y = coord.y + DIST_CELLS[i].y;
        dist = DIST_CELLS[i].dist;

        if (new_x >= 0 && new_x < grid->ncol &&
                new_y >= 0 && new_y < grid->nrow && grid->terrain[new_y * grid->ncol + new_x] == CORRIDOR) {
            Coord neighbor = {new_x, new_y, dist};
            neighbors[k++] = neighbor;
        }
    }
    if (k < num_cells) neighbors[k].x = TURNOFF;    // where to stop iterations
    return neighbors;
}

/**
 * @brief Returns cells in the Moore neighborhood of a given cell.
 * 
 * @param coord The coordinates of the central cell.
 * @param grid The model grid structure.
 * @return An array of coord structures.
 */
Coord* get_neighbors(Coord coord, Grid* grid) {
    int i;
    int num_cells = 8;

    Coord* neighbors = malloc(sizeof(Coord) * num_cells);

    int k = 0;
    int new_x, new_y;
    double dist;

    for (i = 0; i < num_cells; i++) {

        new_x = coord.x + CELLS1[i].x;
        new_y = coord.y + CELLS1[i].y;
        dist = CELLS1[i].dist;

        if (new_x >= 0 && new_x < grid->ncol && new_y >= 0 && new_y < grid->nrow) {
            Coord neighbor = {new_x, new_y, dist};
            neighbors[k++] = neighbor;
        }
    }

    if (k < num_cells) neighbors[k].x = TURNOFF;    // where to stop iterations
    return neighbors;
}

/**
 * @brief Asymptotic threshold model of emigration.
 * 
 * @param n The population.
 * @param phi The threshold as a fraction of the local carrying capacity.
 * @param k The local carrying capacity.
 * @return The number of migrants. 
 */
double fission_ta(double n, double phi, double k) {
    if (phi >= n/k) return 0;
    else return n - (phi * k);
}

/**
 * @brief Applies migration to each populated cell in the model's grid.
 * 
 * @param model The model structure with the fission parameter.
 * @param grid The grid structure on which migration will be applied.
 */
void fission(Model* model, Grid* grid) {
    int i, j;

    int n = model->agent_count; // copy value as it will grow during iteration

    int num_cells = 8;
    int far_cells = NCELL[model->accel - 2];

    for (i = 0; i < n; i++) {
        if (model->active[i] == 0) continue;    // skip inactive cells

        Coord coord = model->agents[i];

        // prevent fission of population in barriers
        if (grid->terrain[coord.y * grid->ncol + coord.x] == BARRIER) continue;

        double n = grid->population[coord.y * grid->ncol + coord.x];
        double local_k = grid->environment[coord.y * grid->ncol + coord.x];
        double migrants = fission_ta(n, model->phi, local_k);

        if (migrants > 0) {
            Coord* nbr;         // stores all neighbors
            Coord* free_nbr;    // stores available neighbors
            int ncell = num_cells;

            if (grid->terrain[coord.y * grid->ncol + coord.x] == CORRIDOR) {
                nbr = get_neighbors_far(coord, grid, model->accel);
                free_nbr = malloc(sizeof(Coord) * far_cells);
                ncell = far_cells;
            } else {
                nbr = get_neighbors(coord, grid);
                free_nbr = malloc(sizeof(Coord) * num_cells);
            }

            int len = 0;        // keep track of how many available cells
            double tot = 0;     // to store the sum of inv sq distances

            for (j = 0; j < ncell; j++) {
                if (nbr[j].x == TURNOFF) {
                    break;
                }

                double nbr_k = grid->environment[nbr[j].y * grid->ncol + nbr[j].x];

                // only add to free neighbors if the environment is positive and
                // population is below the emigration threshold (given the local K)
                if (grid->population[nbr[j].y * grid->ncol + nbr[j].x] < (nbr_k * model->phi) &&
                        grid->environment[nbr[j].y * grid->ncol + nbr[j].x] > 0) {
                    free_nbr[len++] = nbr[j];
                    tot += nbr[j].dist;
                }                    
            }

            if (len > 0) {
                grid->population[coord.y * grid->ncol + coord.x] -= migrants;

                for (j = 0; j < len; j++) {
                    double percent = free_nbr[j].dist / tot;
                    if (grid->arrival[free_nbr[j].y * grid->ncol + free_nbr[j].x] == 0) {
                        grid->arrival[free_nbr[j].y * grid->ncol + free_nbr[j].x] = model->tick;
                        model->active[model->agent_count] = 1;
                        model->agents[model->agent_count++] = free_nbr[j];
                        model->agents[model->agent_count].x = TURNOFF;
                    }
                    // move migrants proportionally to the inv sq distance
                    grid->population[free_nbr[j].y * grid->ncol + free_nbr[j].x] += (migrants * percent);
                }
            } else {
                model->active[i] = 0;
            }

            free(nbr);
            free(free_nbr);
        }
    }
}