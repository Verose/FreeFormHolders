//
// Created by Verose on 12/14/2018.
//

#include "main.h"
#include "model_path.h"


int main(int argc, char *argv[]) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // Load a mesh
    igl::readOBJ(MODEL_PATH "/bunny.obj", V, F);

    igl::opengl::glfw::Viewer viewer;

    const auto update_distance = [&](const int vid) {
        Eigen::VectorXi VS, FS, VT, FT;
        Eigen::VectorXd d;
        // The selected vertex is the source
        VS.resize(1);
        VS << vid;
        // All vertices are the targets
        VT.setLinSpaced(V.rows(), 0, V.rows() - 1);
        std::cout << "Computing geodesic distance to vertex " << vid << "..." << std::endl;
        igl::exact_geodesic(V, F, VS, FS, VT, FT, d);
        const double strip_size = 0.05;
        // The function should be 1 on each integer coordinate
        d = (d / strip_size * igl::PI).array().sin().abs().eval();
        // Compute per-vertex colors
        Eigen::MatrixXd C;
        igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, d, false, C);
        // Plot the mesh
        viewer.data().set_mesh(V, F);
        viewer.data().set_colors(C);
    };

    // Plot a distance when a vertex is picked
    viewer.callback_mouse_down =
            [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool {
                int fid;
                Eigen::Vector3f bc;
                // Cast a ray in the view direction starting from the mouse position
                double x = viewer.current_mouse_x;
                double y = viewer.core.viewport(3) - viewer.current_mouse_y;

                if (igl::unproject_onto_mesh(
                        Eigen::Vector2f(x, y),
                        viewer.core.view,
                        viewer.core.proj,
                        viewer.core.viewport,
                        V,
                        F,
                        fid,
                        bc)) {
                    int max;
                    bc.maxCoeff(&max);
                    int vid = F(fid, max);
                    update_distance(vid);
                    return true;
                }
                return false;
            };
    viewer.data().set_mesh(V, F);

    std::cout << "Click on mesh to define new source.\n" << std::endl;
    update_distance(0);
    viewer.launch();

    //sample_random_point(V, F, viewer);


    return 0;
}

void sample_random_point(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, igl::opengl::glfw::Viewer &viewer) {
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
    viewer.data().set_mesh(V, F);  // copies the mesh into the viewer

    viewer.data().add_points(V_sample, Eigen::RowVector3d(1, 0, 0));  // print sampled point

    viewer.launch();  // creates a window, an OpenGL context and it starts the draw loop

    // Save the mesh
//    writeOBJ(MODEL_PATH "/bunny_new.obj", V, F);

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
    std::cout << "Hello, mesh: " << std::endl << L * E << std::endl;
}
