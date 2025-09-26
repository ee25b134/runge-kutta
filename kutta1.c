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


/* ---------- Heun method for circular motion ---------- */
double HeunMethod(double start_theta, double start_speed,
                  double linear_acceleration, double circle_radius,
                  double time_step, double final_time) {

    // Initialize state
    double x = circle_radius * cos(start_theta);
    double y = circle_radius * sin(start_theta);
    double vx = -start_speed * sin(start_theta);
    double vy =  start_speed * cos(start_theta);

    int steps = (int)(final_time / time_step);

    for (int i = 0; i < steps; i++) {
        // --- Predictor (Euler step) ---
        double r = sqrt(x*x + y*y);
        double tx = -y / r;
        double ty =  x / r;
        double ax = linear_acceleration * tx;
        double ay = linear_acceleration * ty;

        // provisional step
        double x_pred  = x  + vx * time_step;
        double y_pred  = y  + vy * time_step;
        double vx_pred = vx + ax * time_step;
        double vy_pred = vy + ay * time_step;

        // --- Corrector (compute slope at provisional step) ---
        double r_pred = sqrt(x_pred*x_pred + y_pred*y_pred);
        double tx_pred = -y_pred / r_pred;
        double ty_pred =  x_pred / r_pred;
        double ax_pred = linear_acceleration * tx_pred;
        double ay_pred = linear_acceleration * ty_pred;

        // --- Update using average slope ---
        x  += 0.5 * (vx + vx_pred) * time_step;
        y  += 0.5 * (vy + vy_pred) * time_step;
        vx += 0.5 * (ax + ax_pred) * time_step;
        vy += 0.5 * (ay + ay_pred) * time_step;

        // Project back to circle
        double rn = sqrt(x*x + y*y);
        x = x * circle_radius / rn;
        y = y * circle_radius / rn;
    }

    // final angle
    double theta = atan2(y, x);
    if (theta < 0) theta += 2*M_PI;
    return theta;
}
double MidpointMethod(double start_theta, double start_speed,
                      double linear_acceleration, double circle_radius,
                      double time_step, double final_time) {

    // Initialize state
    double x = circle_radius * cos(start_theta);
    double y = circle_radius * sin(start_theta);
    double vx = -start_speed * sin(start_theta);
    double vy =  start_speed * cos(start_theta);

    int steps = (int)(final_time / time_step);

    for (int i = 0; i < steps; i++) {
        // --- Step 1: slope at current state ---
        double r = sqrt(x*x + y*y);
        double tx = -y / r;
        double ty =  x / r;
        double ax = linear_acceleration * tx;
        double ay = linear_acceleration * ty;

        // --- Step 2: midpoint prediction (half Euler step) ---
        double x_mid  = x  + 0.5 * vx * time_step;
        double y_mid  = y  + 0.5 * vy * time_step;
        double vx_mid = vx + 0.5 * ax * time_step;
        double vy_mid = vy + 0.5 * ay * time_step;

        // recompute acceleration at midpoint
        double r_mid = sqrt(x_mid*x_mid + y_mid*y_mid);
        double tx_mid = -y_mid / r_mid;
        double ty_mid =  x_mid / r_mid;
        double ax_mid = linear_acceleration * tx_mid;
        double ay_mid = linear_acceleration * ty_mid;

        // --- Step 3: update using midpoint slope ---
        x  += vx_mid * time_step;
        y  += vy_mid * time_step;
        vx += ax_mid * time_step;
        vy += ay_mid * time_step;

        // Project back to circle
        double rn = sqrt(x*x + y*y);
        x = x * circle_radius / rn;
        y = y * circle_radius / rn;
    }

    // final angle
    double theta = atan2(y, x);
    if (theta < 0) theta += 2*M_PI;
    return theta;
}

#include <stdio.h>
#include <math.h>

// RK4 Method to solve the ODE
double RK4Method(double start_theta, double start_speed, double linear_acceleration, 
                 double circle_radius, double time_step, double final_time) {
    // Initialize state
    double x = circle_radius * cos(start_theta);
    double y = circle_radius * sin(start_theta);
    double vx = -start_speed * sin(start_theta);
    double vy = start_speed * cos(start_theta);

    int steps = (int)(final_time / time_step);

    // RK4 integration loop
    for (int i = 0; i < steps; i++) {
        // Step 1: Calculate acceleration at current position
        double r = sqrt(x * x + y * y); // radius
        double tx = -y / r;  // Tangential unit vector in x direction
        double ty = x / r;   // Tangential unit vector in y direction
        double ax = linear_acceleration * tx;  // Tangential acceleration in x direction
        double ay = linear_acceleration * ty;  // Tangential acceleration in y direction

        // Calculate k1 (current slope)
        double k1x = vx * time_step;
        double k1y = vy * time_step;
        double k1vx = ax * time_step;
        double k1vy = ay * time_step;

        // Step 2: Calculate k2 (midpoint slope, using k1)
        double x_temp = x + 0.5 * k1x;
        double y_temp = y + 0.5 * k1y;
        double vx_temp = vx + 0.5 * k1vx;
        double vy_temp = vy + 0.5 * k1vy;

        r = sqrt(x_temp * x_temp + y_temp * y_temp);
        tx = -y_temp / r;
        ty = x_temp / r;
        ax = linear_acceleration * tx;
        ay = linear_acceleration * ty;

        double k2x = vx_temp * time_step;
        double k2y = vy_temp * time_step;
        double k2vx = ax * time_step;
        double k2vy = ay * time_step;

        // Step 3: Calculate k3 (midpoint slope, using k2)
        x_temp = x + 0.5 * k2x;
        y_temp = y + 0.5 * k2y;
        vx_temp = vx + 0.5 * k2vx;
        vy_temp = vy + 0.5 * k2vy;

        r = sqrt(x_temp * x_temp + y_temp * y_temp);
        tx = -y_temp / r;
        ty = x_temp / r;
        ax = linear_acceleration * tx;
        ay = linear_acceleration * ty;

        double k3x = vx_temp * time_step;
        double k3y = vy_temp * time_step;
        double k3vx = ax * time_step;
        double k3vy = ay * time_step;

        // Step 4: Calculate k4 (next slope, using k3)
        x_temp = x + k3x;
        y_temp = y + k3y;
        vx_temp = vx + k3vx;
        vy_temp = vy + k3vy;

        r = sqrt(x_temp * x_temp + y_temp * y_temp);
        tx = -y_temp / r;
        ty = x_temp / r;
        ax = linear_acceleration * tx;
        ay = linear_acceleration * ty;

        double k4x = vx_temp * time_step;
        double k4y = vy_temp * time_step;
        double k4vx = ax * time_step;
        double k4vy = ay * time_step;

        // Step 5: Update the position and velocity using the average of the slopes
        x += (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        vx += (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6;
        vy += (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6;

        // Project back to the circle to maintain constant radius
        double rn = sqrt(x * x + y * y);
        x = x * circle_radius / rn;
        y = y * circle_radius / rn;
    }

    // Final angular position θ = atan2(y, x)
    double theta = atan2(y, x);
    if (theta < 0) theta += 2 * M_PI;
    return theta;
}



/* ---------- main driver for testing ---------- */
int main() {
    double theta_heun = HeunMethod(
        0.0,    // start_theta
        0.5,    // start_speed
        0.2,    // linear_acceleration
        1.0,    // circle_radius
        0.01,   // time_step
        5.0     // final_time
    );
        double theta_euler = EulerMethod(
        0.0,    
        0.5,  
        0.2,   
        1.0, 
        0.01, 
        5.0   
    );
  double theta_midpoint = MidpointMethod(
        0.0,    
        0.5,  
        0.2,   
        1.0, 
        0.01, 
        5.0   
    );
 
    double theta_RK4Method = RK4Method(
        0.0,    
        0.5,  
        0.2,   
        1.0, 
        0.01, 
        5.0   
    );


    printf("Final angular position θ (Heun) = %.6f rad\n", theta_heun);

    printf("Final angular position θ (euler) = %.6f rad\n", theta_euler);
    
    printf("Final angular position θ (midpoint) = %.6f rad\n", theta_midpoint);
    
    printf("Final angular position θ (RK4method) = %.6f rad\n", theta_RK4Method);
    return 0;
}
