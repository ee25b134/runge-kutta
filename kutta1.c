#include <stdio.h>
#include <math.h>

typedef struct {
    double x;
    double y;
} Vector2;

// Function to compute vector magnitude
double norm(Vector2 v) {
    return sqrt(v.x * v.x + v.y * v.y);
}

// Function to multiply vector by scalar
Vector2 scalar_mul(Vector2 v, double s) {
    Vector2 res = {v.x * s, v.y * s};
    return res;
}

// Function to add two vectors
Vector2 vec_add(Vector2 a, Vector2 b) {
    Vector2 res = {a.x + b.x, a.y + b.y};
    return res;
}

// Function to normalize vector to length R
Vector2 project_to_circle(Vector2 pos, double R) {
    double n = norm(pos);
    Vector2 res = {pos.x * R / n, pos.y * R / n};
    return res;
}

// Euler step function
void euler_step(Vector2 *pos, Vector2 *vel, Vector2 acc, double h, double R) {
    // Step 1: move along velocity
    Vector2 pos_new = vec_add(*pos, scalar_mul(*vel, h));
    
    // Step 2: update velocity
    Vector2 vel_new = vec_add(*vel, scalar_mul(acc, h));
    
    // Step 3: project back to circle
    pos_new = project_to_circle(pos_new, R);
    
    // Update original position and velocity
    *pos = pos_new;
    *vel = vel_new;
}
