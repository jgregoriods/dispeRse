#include <math.h>

void calibc(int* calcurve_c14, int* calcurve_error, int* mu, int* sigma, double *res, int* len) {
    for (int i = 0; i < *len; ++i) {
        double aux = pow(*sigma, 2) + pow(calcurve_error[i], 2);
        res[i] = exp( - ( pow(*mu - calcurve_c14[i], 2) / (2*aux) ) ) / sqrt(2*M_PI*aux);
    }
}