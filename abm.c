#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TURNOFF -1

typedef struct Coord {
    int x;
    int y;
} Coord;

typedef struct Agent {
    Coord coord;
    double population;
    int active;
} Agent;

typedef struct Model {
    int agent_count;
    double K;
    double r;
    double cta;
    int step;
    int generation;
    int dist;
} Model;

typedef struct Grid {
    int height;
    int width;
    int *base;
    double *env;
    int *agent;
} Grid;

int NCELLS[20];
Coord *CELLS[20];

Coord *get_neighborhood(Coord coord, int dist) {
    Coord *neighborhood = malloc(sizeof(Coord) * NCELLS[dist - 1]);
    int index = 0;

    for (int i = 0; i < NCELLS[dist - 1]; i++) {
        Coord neighbor = {coord.x + CELLS[dist - 1][i].x, coord.y + CELLS[dist - 1][i].y};
        neighborhood[index++] = neighbor;
    }

    return neighborhood;
}


// **********************************************************
// Function for logistic growth
//
/*double log_growth_t(double n, double r, double k, int t) {
    return (k * n) / (n + (k - n) * exp(-r * t));
}*/
double log_growth_t(double n, double r, int t) {
    return n / (n + (1 - n) * exp(-r * t));
}

// **********************************************************
// Function for exponential growth
double exp_growth_t(double n, double r, int t) {
    return n * exp(r * t);
}

void grow_population(Agent *agent_list[], Model *model) {
    for (int i = 0; i < model->agent_count; i++) {
        if (!agent_list[i]->active) break;
        //agent_list[i]->population = log_growth_t(agent_list[i]->population, model->r, model->generation);
        agent_list[i]->population = exp_growth_t(agent_list[i]->population, model->r, model->generation);
        if (agent_list[i]->population >= 1) agent_list[i]->population = 1;
    }
}
//
// **********************************************************


Coord choose_cell(Coord *cells, Model *model, Grid *grid) {
    int num = NCELLS[model->dist - 1];
    Coord empty_cells[num];
    for (int i = 0; i < num; i++) {
        empty_cells[i].x = TURNOFF;
        empty_cells[i].y = TURNOFF;
    }
    int j = 0;
    for (int i = 0; i < num; i++) {
        int index = (grid->width * cells[i].y) + cells[i].x;
        if (cells[i].x >= 0 && cells[i].x < grid->width
            && cells[i].y >= 0 && cells[i].y < grid->height
            && !*(grid->agent + index) && *(grid->env + index)) {
            empty_cells[j++] = cells[i];
        }
    }
    if (!j) return empty_cells[0];

    double best = 0;
    int r = 0;
    for (int i = 0; i < j; i++) {
        int index = (grid->width * empty_cells[i].y) + empty_cells[i].x;
        if (*(grid->env + index) > best) {
            best = *(grid->env + index);
            r = i;
        }
    }

    //int r = rand() % j;
    return empty_cells[r];
}


// **********************************************************
// Functions for density-dependent dispersal
//
/*double ta_dispersal(double n, double cta, double k) {
    if (n / k > cta) return (1 - (cta / (n / k)));
    return 0;
}*/
double ta_dispersal(double n, double cta) {
    if (n > cta) return n - cta;
    return 0;
}
//
//
//
void disperse_population(Agent *agent_list[], Model *model, Grid *grid) {
    for (int i = 0; i < model->agent_count; i++) {
        if (!agent_list[i]->active) break;
        double migrants = agent_list[i]->population * ta_dispersal(agent_list[i]->population, model->cta);
        if (migrants) {
            Coord *neighborhood = get_neighborhood(agent_list[i]->coord, model->dist);
            Coord cell = choose_cell(neighborhood, model, grid);
            if (cell.x != TURNOFF && cell.y != TURNOFF) {
                agent_list[i]->population -= migrants;
                
                agent_list[model->agent_count]->coord = cell;
                agent_list[model->agent_count]->population = migrants;
                agent_list[(model->agent_count)++]->active = 1;
                int index = (grid->width * cell.y) + cell.x;
                *(grid->agent + index) = 1;
                *(grid->base + index) = model->step;
            }
            free(neighborhood);
        }
    }
}
//
// **********************************************************


Coord* dist_to_cells(int radius, int* len) {
    Coord* cells = (Coord *)malloc(pow(radius * 2 + 1, 2) * sizeof(Coord));
    int idx = 0;
    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
            double dist = round(hypot(i, j));
            if (dist <= radius && (i || j)) {
                cells[idx].x = i;
                cells[idx].y = j;
                ++idx;
            }
        }
    }
    free(cells);
    Coord* res = (Coord *)malloc(idx * sizeof(Coord));
    for (int i = 0; i < idx; ++i) {
        res[i] = cells[i];
    }
    *len = idx;
    return res;
}

void run_model(int *height, int *width,
               int *base_grid, double *env_grid,
               int *start_x, int *start_y,
               double *k, double *r, double *cta,
               int *dist, int *num_iter) {

    srand(time(NULL));

    int i, j;
    int max_agents = *height * *width;

    int *agent_grid = (int *)malloc(sizeof(int) * max_agents);

    for (i = 0; i < 20; ++i) {
        int len;
        CELLS[i] = dist_to_cells(i + 1, &len);
        NCELLS[i] = len;
    }

    Model model = {0, *k, *r, *cta, 0, 30, *dist};
    Grid grid = {*height, *width, base_grid, env_grid, agent_grid};
    Agent *agent_list[max_agents];

    // Initialize agents
    for (i = 0; i < max_agents; i++) {
        agent_list[i] = malloc(sizeof(Agent));
        agent_list[i]->active = 0;
        grid.agent[i] = 0;
    }

    agent_list[0]->coord.x = *start_x;
    agent_list[0]->coord.y = *start_y;
    agent_list[0]->population = model.cta;
    agent_list[0]->active = 1;
    model.agent_count++;
    int index = (*start_y * (*width)) + *start_x;
    *(agent_grid + index) = 1;
    *(base_grid + index) = 0;

    for (i = 0; i < *num_iter; i++) {
        grow_population(agent_list, &model);
        disperse_population(agent_list, &model, &grid);
        model.step++;
    }

    free(agent_grid);
    for (i = 0; i < max_agents; i++) {
        free(agent_list[i]);

    }
}

int main() {
    int height = 500;
    int width = 500;
    int *grid = (int *)malloc(sizeof(int) * 500 * 500);
    double *env_grid = (double *)malloc(sizeof(double) * 500 * 500);
    for (int i = 0; i < 500 * 500; i++) {
        grid[i] = -1;
        env_grid[i] = 1.0;
    }
    int start_x = 250;
    int start_y = 250;
    int num_iter = 200;
    double k = 1;
    double r = 0.025;
    double cta = 0.5;
    int dist = 3;
    run_model(&height, &width, grid, env_grid, &start_x, &start_y, &k, &r, &cta, &dist, &num_iter);
    free(grid);
    free(env_grid);
    return 0;
}
