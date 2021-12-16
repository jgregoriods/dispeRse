#include <math.h>
#include <stdio.h>

double exp_growth_t(double n, double r, int t) {
    return n * exp(r * t);
}

double exp_growth(double n, double r) {
    return n + r * n;
}

double log_growth_t(double n, double r, double k, int t) {
    return (k*n) / (n + (k - n) * exp(-r * t));
}

double log_growth(double n, double r, double k) {
    return (n + r * n * (1 - n / k));
}

int main(void) {
    double N = 50;
    printf("%.2f\n", log_growth_t(N, 0.025, 200, 500));
    return 0;
}

