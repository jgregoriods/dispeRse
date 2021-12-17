#include <math.h>
#include <stdio.h>

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

int main(void) {
    double N = 50;
    double K = 200;
    double R = 0.025;
    double EK = 0.25;

    for (int i = 0; i < 20; i++) {
        printf("%.2f\n", N);
        N = log_growth_t(N, R, K, 25);
        N -= N * ek_dispersal(N, EK, K, 1);
    }

    return 0;
}

