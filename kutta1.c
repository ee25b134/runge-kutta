#include <stdio.h>
#include <math.h>


double EulerMethod(double start_theta, double start_speed,
                   double linear_acceleration, double circle_radius,
                   double time_step, double final_time) {

    // Initialize state
    double x = circle_radius * cos(start_theta);
    double y = circle_radius * sin(start_theta);
    double vx = -start_speed * sin(start_theta);
    double vy =  start_speed * cos(start_theta);

    int steps = (int)(final_time / time_step);
    for (int i = 0; i < steps; i++) {
        // compute tangential acceleration
        double r = sqrt(x*x + y*y);
        double tx = -y / r;
        double ty =  x / r;
        double ax = linear_acceleration * tx;
        double ay = linear_acceleration * ty;

        // x and v are updated using the euler method.
        x  += vx * time_step;
        y  += vy * time_step;
        vx += ax * time_step;
        vy += ay * time_step;

        // projecting back to circle by adjusting the distance.
        double rn = sqrt(x*x + y*y);
        x = x * circle_radius / rn;
        y = y * circle_radius / rn;
    }

    // to keep theta strictly above 0, we are using a if clause to add 2pi.
    double theta = atan2(y, x);
    if (theta < 0) theta += 2*M_PI;
    return theta;
}


int main() {
    double theta_final = EulerMethod(
        0.0,    
        0.5,  
        0.2,   
        1.0, 
        0.01, 
        5.0   
    );

    printf("Final angular position θ = %.6f rad\n", theta_final);
    return 0;
}
