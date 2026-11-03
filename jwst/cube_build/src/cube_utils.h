#ifndef CUBE_UTILS_H
#define CUBE_UTILS_H

// Memory allocation function declarations
int
alloc_flux_arrays(int nelem, double **fluxv, double **weightv, double **varv, double **ifluxv);

int
alloc_flux_dq_arrays(
    int nelem, double **fluxv, double **weightv, double **varv, double **ifluxv, int **dqv);

// Geometry and overlap function declarations
double
find_area_quad(double MinX, double MinY, double Xcorner[], double Ycorner[]);

double
find_area_poly(int nVertices, double xPixel[], double yPixel[]);

// Geometry / Overlap declarations
double
sh_find_overlap(
    double xcenter, double ycenter, double xlength, double ylength, double xPixelCorner[],
    double yPixelCorner[]);

#endif /* CUBE_UTILS_H */
