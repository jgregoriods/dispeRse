#include <math.h>

#include "dispeRse.h"

/**
 * @brief Logistic growth model.
 * 
 * @param n Initial population.
 * @param r Annual growth rate.
 * @param k Carrying capacity.
 * @param t Number of years of growth.
 * @return The new population after t years.
 */
double log_growth(double n, double r, double k, double t) {
    return (k*n) / ((k-n) * exp(-r*t) + n);
}

/**
 * @brief Applies growth to each populated cell in the model's grid. A local
 * carrying capacity and growth rate are obtained from the grid's environment.
 * 
 * @param model The model structure with the growth parameters.
 * @param grid The grid structure to which growth will be applied.
 */
void grow(Model* model, Grid* grid) {
    for (int i = 0; i < model->agent_count; i++) {
        double n = grid->population[model->agents[i].y * grid->ncol + model->agents[i].x];
        double local_k = grid->environment[model->agents[i].y * grid->ncol + model->agents[i].x];
        double local_r = model->r * pow(grid->environment[model->agents[i].y * grid->ncol + model->agents[i].x], model->gamma);
        grid->population[model->agents[i].y * grid->ncol + model->agents[i].x] = log_growth(n, local_r, local_k, model->t);
    }
}
