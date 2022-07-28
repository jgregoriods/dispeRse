#include <stdlib.h>
#include <math.h>

#include "main.h"

void run_model(int* nrow, int* ncol, double* population, double* env, int* arrival,
               double* r, double* phi, int* start,
               int* x, int* y, int* iter, int* num_origins, double* t, int* terrain, int* accel, double* gamma) {
    int i, j;
    int max_agents = *nrow * *ncol * 2;
    double local_k;
    Grid grid;

    grid.population = population;
    grid.arrival = arrival;

    grid.env = env;
    grid.terr = terrain;
    grid.nrow = *nrow;
    grid.ncol = *ncol;

    Model model = {*r, *phi, *t};

    int max_age = start[0];
    for (i = 0; i < *num_origins; ++i) {
        if (start[i] > max_age) {
            max_age = start[i];
        }
    }

    model.agents = malloc(sizeof(Coord) * max_agents);
    model.active = malloc(sizeof(int) * max_agents);
    for (i = 0; i < max_agents; ++i) model.active[i] = 0;
    model.tick = max_age;

    model.gamma = *gamma;
    model.accel = *accel;

    int active_agents = 0;
    for (i = 0; i < *num_origins; ++i) {
        if (start[i] == max_age) {
            Coord init_cell = {x[i], y[i]};
            model.active[active_agents] = 1;
            model.agents[active_agents++] = init_cell;
            grid.arrival[y[i] * *ncol + x[i]] = max_age;
            local_k = grid.env[y[i] * *ncol + x[i]];
            grid.population[y[i] * *ncol + x[i]] = local_k;
        }
    }
    model.agents[active_agents].x = TURNOFF;
    model.agent_count = active_agents;

    for (i = 0; i < *iter; ++i) {
        grow(&model, &grid);
        fission(&model, &grid);
        model.tick -= model.t;

        for (j = 1; j < *num_origins; ++j) {
            if (start[j] >= model.tick && grid.population[y[j] * *ncol + x[j]] == 0) {
                Coord init_cell = {x[j], y[j]};
                model.active[model.agent_count] = 1;
                model.agents[model.agent_count++] = init_cell;
                grid.arrival[y[j] * *ncol + x[j]] = model.tick;
                local_k = grid.env[y[j] * *ncol + x[j]];
                grid.population[y[j] * *ncol + x[j]] = local_k;
                model.agents[model.agent_count].x = TURNOFF;
            }
        }
    }

    free(model.agents);
}
