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
    int population;
    int active;
} Agent;

typedef struct Model {
    int agent_count;
    int K;
    double r;
    double epsilon;
    int step;
} Model;

typedef struct Grid {
    int height;
    int width;
    int *base;
    double *env;
    int *agent;
} Grid;

const Coord KNIGHT[16] = {
    {-1, 2}, {1, 2},
    {-2, 1}, {-1, 1}, {0, 1}, {1, 1}, {2, 1},
    {-1, 0}, {1, 0},
    {-2, -1}, {-1, -1}, {0, -1}, {1, -1}, {2, -1},
    {-1, -2}, {1, -2}
};

Coord *get_neighborhood(Coord coord) {
    Coord *neighborhood = malloc(sizeof(Coord) * 16);
    int index = 0;
    //for (int i = -1; i < 2; i++) {
    //    for (int j = -1; j < 2; j++) {
    //        if (i || j) {
    for (int i = 0; i < 16; i++) {
        Coord neighbor = {coord.x + KNIGHT[i].x, coord.y + KNIGHT[i].y};
        neighborhood[index++] = neighbor;
    }
    //    }
    //}
    return neighborhood;
}

// Logistic growth
void grow_population(Agent *agent_list[], Model *model) {
    for (int i = 0; i < model->agent_count; i++) {
        if (!agent_list[i]->active) break;
        int N = agent_list[i]->population;
        agent_list[i]->population += round(model->r * (double)N * (1 - ((double)N / (double)model->K)));
    }
}

Coord choose_cell(Coord cells[16], Model *model, Grid *grid) {
    Coord empty_cells[16];
    for (int i = 0; i < 16; i++) {
        empty_cells[i].x = TURNOFF;
        empty_cells[i].y = TURNOFF;
    }
    int j = 0;
    for (int i = 0; i < 16; i++) {
        int index = (grid->width * cells[i].y) + cells[i].x;
        if (cells[i].x >= 0 && cells[i].x < grid->width
            && cells[i].y >= 0 && cells[i].y < grid->height
            && !*(grid->agent + index) && *(grid->env + index)) {
            empty_cells[j++] = cells[i];
        }
    }
    if (!j) return empty_cells[0];
    int r = rand() % j;
    return empty_cells[r];
}

// Density-dependent emigration
void disperse_population(Agent *agent_list[], Model *model, Grid *grid) {
    for (int i = 0; i < model->agent_count; i++) {
        if (!agent_list[i]->active) break;
        int N = agent_list[i]->population;
        int migrants = round((double)N * model->epsilon * ((double)N / (double)model->K));
        if (migrants > 150) {
            Coord *neighborhood = get_neighborhood(agent_list[i]->coord);
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

void run_model(int *height, int *width,
               int *base_grid, double *env_grid,
               int *start_x, int *start_y,
               int *num_iter) {
    srand(time(NULL));

    int i, j;
    int max_agents = *height * *width;

    int *agent_grid = (int *)malloc(sizeof(int) * max_agents);

    Model model = {0, 500, 0.05, 0.6, 0};
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
    agent_list[0]->population = 500;
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

/*
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
    int num_iter = 6000;
    run_model(&height, &width, grid, env_grid, &start_x, &start_y, &num_iter);
    free(grid);
    free(env_grid);
    return 0;
}
*/