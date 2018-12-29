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


void calc_grip(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &d,
               std::vector<Eigen::Vector2i> &cuts);

bool is_other_edge(const int v1, const int v2, const double t, const Eigen::VectorXd &d, int &smaller_v);

int get_closer_v_id_from_source(const Eigen::Vector2i e, const Eigen::VectorXd &d);

void save_grip_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &d, const double t);

void sample_random_point(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, igl::opengl::glfw::Viewer &viewer);

#endif //FREEFROMHOLDERS_MAIN_H
