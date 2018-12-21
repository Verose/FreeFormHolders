//
// Created by Verose on 12/14/2018.
//

#include "../include/main.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/random_points_on_mesh.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <iostream>

#include "model_path.h"


using namespace std;


int main(int argc, char *argv[]) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // Load a mesh
    igl::readOBJ(MODEL_PATH "/bunny.obj", V, F);

    // Print the vertices and faces matrices
    std::cout << "Vertices: " << std::endl << V << std::endl;
    std::cout << "Faces:    " << std::endl << F << std::endl;

    // sample 1 random point on mesh
    // given by barycentric coordinates of the sampled point in face FI
    Eigen::VectorXi FI;
    Eigen::SparseMatrix<double> B;
    igl::random_points_on_mesh(1, V, F, B, FI);

    // convert barycenter coordinates to original sampled point
    Eigen::MatrixXd V_sample = B * V;

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);  // copies the mesh into the viewer

    viewer.data().add_points(V_sample, Eigen::RowVector3d(1, 0, 0));  // print sampled point

    viewer.launch();  // creates a window, an OpenGL context and it starts the draw loop

    // Save the mesh
//    igl::writeOBJ(MODEL_PATH "/bunny_new.obj", V, F);


    return 0;
}

void draw_mesh() {// Inline mesh of a cube
    const Eigen::MatrixXd V = (Eigen::MatrixXd(8, 3) <<
            0.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
            0.0, 1.0, 0.0,
            0.0, 1.0, 1.0,
            1.0, 0.0, 0.0,
            1.0, 0.0, 1.0,
            1.0, 1.0, 0.0,
            1.0, 1.0, 1.0).finished();
    const Eigen::MatrixXi F = (Eigen::MatrixXi(12, 3) <<
            1, 7, 5,
            1, 3, 7,
            1, 4, 3,
            1, 2, 4,
            3, 8, 7,
            3, 4, 8,
            5, 7, 8,
            5, 8, 6,
            1, 5, 6,
            1, 6, 2,
            2, 6, 8,
            2, 8, 4).finished().array() - 1;

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();
}

void hello_mesh() {
    Eigen::MatrixXd E(4, 2);
    E << 0, 0,
         1, 0,
         1, 1,
         0, 1;
    Eigen::MatrixXi G(2, 3);
    G << 0, 1, 2,
         0, 2, 3;
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(E, G, L);
    cout << "Hello, mesh: " << endl << L * E << endl;
}
