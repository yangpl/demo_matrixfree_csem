#include <stdio.h>
#include <stdlib.h>

/*
 * Conservative averaging from source grid (xe, f) to target grid (Xe)
 * Both grids must cover the same total domain.
 *
 * Inputs:
 *   Ns  - number of source cells
 *   Nt  - number of target cells
 *   xe  - array of length Ns+1 (source cell edges)
 *   f   - array of length Ns   (source cell averages)
 *   Xe  - array of length Nt+1 (target cell edges)
 *
 * Output:
 *   F   - array of length Nt   (target cell averages)
 *
 * Conservation property:
 *   sum_i f[i]*(xe[i+1]-xe[i]) == sum_j F[j]*(Xe[j+1]-Xe[j])
 */

void conservative_average_same_length(
    int Ns, int Nt,
    const double *xe, const double *f,
    const double *Xe, double *F)
{
    int i = 0;
    for (int j = 0; j < Nt; j++) {
        double X1 = Xe[j];
        double X2 = Xe[j + 1];
        double dx_target = X2 - X1;
        double s = 0.0;

        /* advance i to skip source cells entirely left of this target cell */
        while (i < Ns && xe[i + 1] <= X1)
            i++;

        int k = i;
        /* accumulate overlaps with source cells intersecting [X1, X2] */
        while (k < Ns && xe[k] < X2) {
            double overlap = 0.0;
            double left = (xe[k] > X1) ? xe[k] : X1;
            double right = (xe[k + 1] < X2) ? xe[k + 1] : X2;
            overlap = right - left;
            if (overlap > 0.0) s += f[k] * overlap;
            k++;
        }

        F[j] = s / dx_target;
    }
}

/* ================= Example Usage ================= */
int main(void)
{
    int Ns = 3;
    int Nt = 4;
    double xe[] = {0.0, 1.0, 2.5, 4.0};//length=Ns+1
    double f[]  = {2.0, 3.0, 1.0};//length=Ns
    double Xe[] = {0.0, 1.5, 2.0, 3.0, 4.0};//length=Nt+1
    double F[3];

    conservative_average_same_length(Ns, Nt, xe, f, Xe, F);

    /* check conservation */
    double source_total = 0.0, target_total = 0.0;
    for (int i = 0; i < Ns; i++) source_total += f[i] * (xe[i + 1] - xe[i]);
    for (int j = 0; j < Nt; j++) target_total += F[j] * (Xe[j + 1] - Xe[j]);

    printf("Target averages:\n");
    for (int j = 0; j < Nt; j++) printf("  F[%d] = %g\n", j, F[j]);

    printf("Source integral = %g\n", source_total);
    printf("Target integral = %g\n", target_total);

    return 0;
}
