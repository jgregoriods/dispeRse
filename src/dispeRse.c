#include <stdlib.h>
#include <math.h>

#include "dispeRse.h"

/**
 * @brief Simulates population dispersal from one or more geographical origins
 * using density-dependent growth and migration. Carrying capacity, growth rate
 * and mobility are allowed to vary with the environment and terrain.
 *
 * @param nrow The number of rows in the environment and terrain grids.
 * @param ncol The number of colunms in the environment and terrain grids.
 * @param environment An array. Result of flattening a 2d grid. Environmental
 * values that affect carrying capacity (K) and growth rate. Typically given as
 * a fraction [0,1] of the max K. Notice that the array length can be a multiple
 * of nrow * ncol if the environment is to be updated during the experiment.
 * @param terrain An array. Result of flattening a 2d grid. Cells with value 1
 * are barriers (prevent fission). Cells with value 2 are corridors - migrants
 * from those cells can move longer distances than the Moore neighborhood along
 * the corridor. Notice that the array length can be a multiple of nrow * ncol
 * if the terrain is to be updated during the experiment.
 * @param population An array that will store the population of each cell as a
 * fraction of max K.
 * @param arrival An array that will store the arrival time at each cell.
 * @param x An array with the X coordinates of each origin.
 * @param y An array with the Y coordinates of each origin.
 * @param start An array with the start dates of each origin.
 * @param num_origins Number of origins.
 * @param num_iter Number of steps to run the model for.
 * @param r The annual growth rate as a decimal.
 * @param phi The fission threshold as a fraction of (local) carrying capacity.
 * Migration starts above this threshold.
 * @param t The duration, in years, of a generation (model time step).
 * @param accel The factor by which the usual miration distance is increased
 * along corridors. Must range from 2 to 4.
 * @param gamma A power that controls the shape of the dependency between r and
 * the environment.
 * @param updates An array specifying the model steps at which the environment
 * and terrain grids will be updated.
 */
void run_model(int *nrow, int *ncol, double *environment, int *terrain,
               double *population, int *arrival, int *x, int *y, int *start,
               int *num_origins, int *num_iter, double *r, double *phi,
               double *t, int *accel, double *gamma, int *updates) {
    int i, j;
    int ncell = *nrow * *ncol;
    int max_agents = *nrow * *ncol * 2;
    double local_k;
    int index = 0;  // to keep track of environment updates

    Grid grid = {*nrow, *ncol, population, environment, terrain, arrival};
    Model model = {*r, *phi, *t, *accel, *gamma};

    // init all cells with start time = max start time
    int max_age = start[0];
    for (i = 0; i < *num_origins; i++) {
        if (start[i] > max_age) {
            max_age = start[i];
        }
    }

    model.agents = malloc(sizeof(Coord) * max_agents);

    // keep track of inactive cells (those that cannot fission)
    model.active = malloc(sizeof(int) * max_agents);
    for (i = 0; i < max_agents; i++) model.active[i] = 0;

    model.tick = max_age;

    int active_agents = 0;
    for (i = 0; i < *num_origins; i++) {
        if (start[i] == max_age) {
            Coord init_cell = {x[i], y[i]};
            model.active[active_agents] = 1;
            model.agents[active_agents++] = init_cell;
            grid.arrival[y[i] * *ncol + x[i]] = max_age;
            local_k = pow(grid.environment[y[i] * *ncol + x[i]], model.gamma);
            grid.population[y[i] * *ncol + x[i]] = local_k;
        }
    }
    model.agents[active_agents].x = TURNOFF;    // where to stop iterations
    model.agent_count = active_agents;          // keep track of how many settled cells

    // main loop
    for (i = 0; i < *num_iter; i++) {

        grow(&model, &grid);
        fission(&model, &grid);
        model.tick -= model.t;

        // check if other origin cells should start their expansion
        // at the current time step
        for (j = 1; j < *num_origins; j++) {
            // only if the cell has not already been reached
            if (start[j] >= model.tick && grid.population[y[j] * *ncol + x[j]] == 0) {
                Coord init_cell = {x[j], y[j]};
                model.active[model.agent_count] = 1;
                model.agents[model.agent_count++] = init_cell;
                grid.arrival[y[j] * *ncol + x[j]] = model.tick;
                local_k = pow(grid.environment[y[j] * *ncol + x[j]], model.gamma);
                grid.population[y[j] * *ncol + x[j]] = local_k;
                model.agents[model.agent_count].x = TURNOFF;
            }
        }

        // check if environment must be updated
        if (i == updates[index]) {
            index++;
            grid.environment = malloc(sizeof(double) * ncell);
            grid.terrain = malloc(sizeof(int) * ncell);
            for (j = 0; j < ncell; j++) {
                grid.environment[j] = environment[(ncell * index)+j];
                grid.terrain[j] = terrain[(ncell * index)+j];
            }
            // reactivate all cells, as they may be able to fission now
            for (j = 0; j < model.agent_count; j++) {
                model.active[j] = 1;
            }
        }
    }

    free(model.agents);
}
