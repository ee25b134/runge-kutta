#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Tangential projection helper
void project_to_circle(double *x, double *y, double *vx, double *vy, double R) {
    double r = sqrt((*x)*(*x) + (*y)*(*y));
    if (r == 0.0) return;
    *x = (*x) * R / r;
    *y = (*y) * R / r;

    // Project velocity back onto tangent
    double tx = -(*y) / R;
    double ty =  (*x) / R;
    double vt = (*vx)*tx + (*vy)*ty;
    *vx = vt * tx;
    *vy = vt * ty;
}

// ===== Solver Implementations =====

void EulerMethod(double *x, double *y, double *vx, double *vy,
                 double a_t, double R, double dt) {
    double r = sqrt((*x)*(*x) + (*y)*(*y));
    double tx = -(*y) / r;
    double ty =  (*x) / r;
    double ax = a_t * tx;
    double ay = a_t * ty;

    *x  += (*vx) * dt;
    *y  += (*vy) * dt;
    *vx += ax * dt;
    *vy += ay * dt;

    project_to_circle(x, y, vx, vy, R);
}

void HeunMethod(double *x, double *y, double *vx, double *vy,
                double a_t, double R, double dt) {
    double r = sqrt((*x)*(*x) + (*y)*(*y));
    double tx = -(*y) / r;
    double ty =  (*x) / r;
    double ax = a_t * tx;
    double ay = a_t * ty;

    double x_pred = *x + (*vx) * dt;
    double y_pred = *y + (*vy) * dt;
    double vx_pred = *vx + ax * dt;
    double vy_pred = *vy + ay * dt;

    double r_pred = sqrt(x_pred*x_pred + y_pred*y_pred);
    double tx_pred = -y_pred / r_pred;
    double ty_pred =  x_pred / r_pred;
    double ax_pred = a_t * tx_pred;
    double ay_pred = a_t * ty_pred;

    *x  += 0.5 * ((*vx) + vx_pred) * dt;
    *y  += 0.5 * ((*vy) + vy_pred) * dt;
    *vx += 0.5 * (ax + ax_pred) * dt;
    *vy += 0.5 * (ay + ay_pred) * dt;

    project_to_circle(x, y, vx, vy, R);
}

void MidpointMethod(double *x, double *y, double *vx, double *vy,
                    double a_t, double R, double dt) {
    double r = sqrt((*x)*(*x) + (*y)*(*y));
    double tx = -(*y) / r;
    double ty =  (*x) / r;
    double ax = a_t * tx;
    double ay = a_t * ty;

    double x_mid  = *x + 0.5*(*vx)*dt;
    double y_mid  = *y + 0.5*(*vy)*dt;
    double vx_mid = *vx + 0.5*ax*dt;
    double vy_mid = *vy + 0.5*ay*dt;

    double r_mid = sqrt(x_mid*x_mid + y_mid*y_mid);
    double tx_mid = -y_mid / r_mid;
    double ty_mid =  x_mid / r_mid;
    double ax_mid = a_t * tx_mid;
    double ay_mid = a_t * ty_mid;

    *x  += vx_mid * dt;
    *y  += vy_mid * dt;
    *vx += ax_mid * dt;
    *vy += ay_mid * dt;

    project_to_circle(x, y, vx, vy, R);
}

void RK4Method(double *x, double *y, double *vx, double *vy,
               double a_t, double R, double dt) {
    double r, tx, ty, ax, ay;

    // k1
    r = sqrt((*x)*(*x) + (*y)*(*y));
    tx = -(*y) / r;  ty = (*x) / r;
    ax = a_t * tx;   ay = a_t * ty;

    double k1x = (*vx) * dt;
    double k1y = (*vy) * dt;
    double k1vx = ax * dt;
    double k1vy = ay * dt;

    // k2
    double x2 = *x + 0.5*k1x;
    double y2 = *y + 0.5*k1y;
    double vx2 = *vx + 0.5*k1vx;
    double vy2 = *vy + 0.5*k1vy;
    r = sqrt(x2*x2 + y2*y2);
    tx = -y2 / r; ty = x2 / r;
    ax = a_t * tx; ay = a_t * ty;

    double k2x = vx2 * dt;
    double k2y = vy2 * dt;
    double k2vx = ax * dt;
    double k2vy = ay * dt;

    // k3
    double x3 = *x + 0.5*k2x;
    double y3 = *y + 0.5*k2y;
    double vx3 = *vx + 0.5*k2vx;
    double vy3 = *vy + 0.5*k2vy;
    r = sqrt(x3*x3 + y3*y3);
    tx = -y3 / r; ty = x3 / r;
    ax = a_t * tx; ay = a_t * ty;

    double k3x = vx3 * dt;
    double k3y = vy3 * dt;
    double k3vx = ax * dt;
    double k3vy = ay * dt;

    // k4
    double x4 = *x + k3x;
    double y4 = *y + k3y;
    double vx4 = *vx + k3vx;
    double vy4 = *vy + k3vy;
    r = sqrt(x4*x4 + y4*y4);
    tx = -y4 / r; ty = x4 / r;
    ax = a_t * tx; ay = a_t * ty;

    double k4x = vx4 * dt;
    double k4y = vy4 * dt;
    double k4vx = ax * dt;
    double k4vy = ay * dt;

    *x  += (k1x + 2*k2x + 2*k3x + k4x)/6.0;
    *y  += (k1y + 2*k2y + 2*k3y + k4y)/6.0;
    *vx += (k1vx + 2*k2vx + 2*k3vx + k4vx)/6.0;
    *vy += (k1vy + 2*k2vy + 2*k3vy + k4vy)/6.0;

    project_to_circle(x, y, vx, vy, R);
}

// ===== Simulation wrapper =====
void simulate(const char *solver_name,
              double start_theta, double start_speed,
              double a_t, double R,
              double dt, double final_time,
              double *theta_traj, int steps) {
    // Initial conditions
    double x = R * cos(start_theta);
    double y = R * sin(start_theta);
    double vx = -start_speed * sin(start_theta);
    double vy =  start_speed * cos(start_theta);

    for (int i=0; i<steps; i++) {
        if (strcmp(solver_name, "Euler") == 0) {
            EulerMethod(&x, &y, &vx, &vy, a_t, R, dt);
        } else if (strcmp(solver_name, "Heun") == 0) {
            HeunMethod(&x, &y, &vx, &vy, a_t, R, dt);
        } else if (strcmp(solver_name, "Midpoint") == 0) {
            MidpointMethod(&x, &y, &vx, &vy, a_t, R, dt);
        } else if (strcmp(solver_name, "RK4") == 0) {
            RK4Method(&x, &y, &vx, &vy, a_t, R, dt);
        } else {
            fprintf(stderr, "Error: Invalid solver name '%s'.\n", solver_name);
            exit(1);
        }

        double theta = atan2(y,x);
        if (theta < 0) theta += 2*M_PI;
        theta_traj[i] = theta;
    }
}

// ===== Main =====
int main(int argc, char *argv[]) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s start_theta start_speed linear_acceleration circle_radius time_step final_time solver_name\n", argv[0]);
        return 1;
    }

    double start_theta = atof(argv[1]);
    double start_speed = atof(argv[2]);
    double a_t         = atof(argv[3]);
    double R           = atof(argv[4]);
    double dt          = atof(argv[5]);
    double final_time  = atof(argv[6]);
    char *solver_name  = argv[7];

    int steps = (int)round(final_time / dt);

    double *traj_solver = malloc(steps * sizeof(double));
    double *traj_rk4    = malloc(steps * sizeof(double));

    simulate(solver_name, start_theta, start_speed, a_t, R, dt, final_time, traj_solver, steps);
    simulate("RK4",      start_theta, start_speed, a_t, R, dt, final_time, traj_rk4, steps);

    // Compute RMS error
    double sumsq = 0.0;
    for (int i=0; i<steps; i++) {
        double diff = traj_solver[i] - traj_rk4[i];
        sumsq += diff*diff;
    }
    double rms = sqrt(sumsq/steps);

    printf("Final angle (%s) = %.12f radians\n", solver_name, traj_solver[steps-1]);
    printf("Final angle (RK4) = %.12f radians\n", traj_rk4[steps-1]);
    printf("RMS error vs RK4 = %.12e\n", rms);

    free(traj_solver);
    free(traj_rk4);
    return 0;
}
