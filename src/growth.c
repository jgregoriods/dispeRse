#include <math.h>

#include "main.h"

double log_growth(double n, double r, double k, double t) {
    return (k*n) / ((k-n) * exp(-r*t) + n);
}

void grow(Model* model, Grid* grid) {
    for (int i = 0; i < model->agent_count; ++i) {
        double n = grid->population[model->agents[i].y * grid->ncol + model->agents[i].x];
        double local_k = grid->env[model->agents[i].y * grid->ncol + model->agents[i].x];
        double local_r = model->r * pow(grid->env[model->agents[i].y * grid->ncol + model->agents[i].x], model->gamma);
        grid->population[model->agents[i].y * grid->ncol + model->agents[i].x] = log_growth(n, local_r, local_k, model->t);
    }
}