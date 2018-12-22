//
// Created by Verose on 12/14/2018.
//

#ifndef FREEFROMHOLDERS_MAIN_H
#define FREEFROMHOLDERS_MAIN_H

#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/random_points_on_mesh.h>
#include <igl/exact_geodesic.h>
#include <igl/colormap.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/PI.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <iostream>


void hello_mesh();

void draw_mesh();

void sample_random_point(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, igl::opengl::glfw::Viewer &viewer);

#endif //FREEFROMHOLDERS_MAIN_H
