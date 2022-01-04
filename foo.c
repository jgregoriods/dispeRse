#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define TURNOFF -1

typedef struct Agent {
    Coord coord;
    double population;
} Agent;

typedef struct Coord {
    int x;
    int y;
} Coord;

typedef struct Grid {
    int width;
    int height;
    int* agent;
} Grid;

double exp_growth(double n, double r) {
    return n + r * n;
}

double exp_growth_t(double n, double r, int t) {
    return n * exp(r * t);
}

double log_growth(double n, double r, double k) {
    return (n + r * n * (1 - n / k));
}

double log_growth_t(double n, double r, double k, int t) {
    return (k * n) / (n + (k - n) * exp(-r * t));
}

double ek_dispersal(double n, double ek, double k, double gamma) {
    return ek * pow((n / k), gamma);
}

double ta_dispersal(double n, double cta, double k) {
    if (n / k > cta) return (1 - (cta / (n / k)));
    return 0;
}

void print_grid(int* grid, int nrow, int ncol) {
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            printf("%d ", grid[i * ncol + j]);
        }
        printf("\n");
    }
}

Coord* dist_to_cells(int radius, int* len) {
    Coord* cells = (Coord *)malloc(pow(radius * 2 + 1, 2) * sizeof(Coord));
    int idx = 0;
    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
            double dist = hypot(i, j);
            if (dist <= radius && (i || j)) {
                cells[idx].x = i;
                cells[idx].y = j;
                ++idx;
            }
        }
    }
    Coord* res = (Coord *)malloc(idx * sizeof(Coord));
    for (int i = 0; i < idx; ++i) {
        res[i] = cells[i];
    }
    *len = idx;
    return res;
}

Coord* cells_in_radius(Coord coord, Grid* grid, int radius) {
    int ncell = pow(2 * radius + 1, 2);
    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
            int newx = x + i;
            int newy = y + j;
            if ((i || j) && newx )
        }
    }
}

int main(void) {

    Coord* cells[3];
    int ncells[3];

    for (int i = 0; i < 3; ++i) {
        int len;
        cells[i] = dist_to_cells(i+1, &len);
        ncells[i] = len;
    }

    Grid grid = {10, 10};
    grid.agent = (int *)malloc(100 * sizeof(int));

    int a;
    while (1) {
        scanf("%d", &a);
        for (int i = 0; i < ncells[a - 1]; ++i) {
            printf("%d %d\n", cells[a-1][i].x, cells[a-1][i].y);
        }
    }

    Agent *agents = (Agent *)malloc(10 * sizeof(Agent));

    for (int i = 0; i < 10; ++i) {
        agents[i].population = TURNOFF;
    }

    agents[0].coord.x = 5;
    agents[0].coord.y = 5;
    agents[0].population = 0.1;

    grid[5 * ncol + 5] = 1;

    return 0;
}
