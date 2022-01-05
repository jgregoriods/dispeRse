#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define TURNOFF -1


typedef struct Coord {
    int x;
    int y;
} Coord;

typedef struct Agent {
    Coord coord;
    double population;
} Agent;

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
            double dist = round(hypot(i, j));
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

int main(void) {
    int ncells[3];
    Coord *cells[3];

    double n, cta, k;

    while (1) {
        scanf("%lf %lf %lf", &n, &cta, &k);
        printf("%.2f\n\n", ta_dispersal(n, cta, k));
    }

    return 0;
}
