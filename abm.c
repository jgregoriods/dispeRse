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
    int height;
    int width;
    int *grid;
} Model;

Coord *get_neighborhood(Coord coord) {
    Coord *neighborhood = malloc(sizeof(Coord) * 8);
    int index = 0;
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            if (i || j) {
                Coord neighbor = {coord.x + i, coord.y + j};
                neighborhood[index++] = neighbor;
            }
        }
    }
    return neighborhood;
}

// Logistic growth
void grow_population(Agent *agent_list[], float r, int K, int *agent_count) {
    for (int i = 0; i < *agent_count; i++) {
        if (!agent_list[i]->active) break;
        int N = agent_list[i]->population;
        agent_list[i]->population += round(r * (float)N * (1 - ((float)N / (float)K)));
    }
}

Coord choose_cell(Coord cells[8], Model *model) {
    Coord empty_cells[8];
    for (int i = 0; i < 8; i++) {
        empty_cells[i].x = TURNOFF;
        empty_cells[i].y = TURNOFF;
    }
    int j = 0;
    for (int i = 0; i < 8; i++) {
        int index = (model->width * cells[i].y) + cells[i].x;
        if (cells[i].x >= 0 && cells[i].x < model->width && cells[i].y >= 0 && cells[i].y < model->height && !*(model->grid + index)) {
            empty_cells[j++] = cells[i];
        }
    }
    if (!j) return empty_cells[0];
    int r = rand() % j;
    return empty_cells[r];
}

// Density-dependent growth
void disperse_population(Agent *agent_list[], int fiss, int *agent_count, Model *model) {
    for (int i = 0; i < *agent_count; i++) {
        if (!agent_list[i]->active) break;
        if (agent_list[i]->population > fiss) {
            Coord *neighborhood = get_neighborhood(agent_list[i]->coord);
            Coord cell = choose_cell(neighborhood, model);
            if (cell.x != TURNOFF && cell.y != TURNOFF) {
                agent_list[i]->population /= 2;
                
                agent_list[*agent_count]->coord = cell;
                agent_list[*agent_count]->population = agent_list[i]->population;
                agent_list[(*agent_count)++]->active = 1;
                int index = (model->width * cell.y) + cell.x;
                *(model->grid + index) = 1;
            }
        }
    }
}

void run_model(int *height, int *width, int *grid, int *num_iter) {
    srand(time(NULL));

    int i, j;
    int max_agents = *height * *width;

    int agent_count = 0;

    Model model = {*height, *width, grid};

    Agent *agent_list[max_agents];

    // Initialize agents
    for (i = 0; i < max_agents; i++) {
        agent_list[i] = malloc(sizeof(Agent));
        agent_list[i]->active = 0;
    }

    Coord xy = {12, 25};
    agent_list[0]->coord = xy;
    agent_list[0]->population = 100;
    agent_list[0]->active = 1;
    agent_count++;
    int index = (25 * (*height)) + 12;
    *(grid + index) = 1;

    for (i = 0; i < *num_iter; i++) {
        grow_population(agent_list, 0.02, 500, &agent_count);
        disperse_population(agent_list, 120, &agent_count, &model);
    }
}