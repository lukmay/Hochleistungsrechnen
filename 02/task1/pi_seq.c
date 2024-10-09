#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double calculate_pi(int num_points) {
    int points_in_circle = 0;

    for (int count = 0; count < num_points; count++) {
        double random_x = (double)rand() / RAND_MAX;
        double random_y = (double)rand() / RAND_MAX;

        if (random_x * random_x + random_y * random_y <= 1.0) {
            points_in_circle++;
        }
    }

    return 4.0 * points_in_circle / num_points;
}

int main(int argc, char *argv[]) {
    int total_points = atoi(argv[1]);

    double pi_approx = calculate_pi(total_points);

    printf("Approximation of Pi: %f\n", pi_approx);

    return 0;
}
