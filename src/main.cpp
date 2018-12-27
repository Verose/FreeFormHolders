//
// Created by Verose on 12/14/2018.
//
#include <igl/triangle_triangle_adjacency.h>
#include <igl/adjacency_matrix.h>
#include <igl/HalfEdgeIterator.h>
#include <igl/cut_mesh.h>

#include <array>
#include <set>

#include "main.h"
#include "model_path.h"

//VT = v1 v2 v3
//FT = f1
//d  = 1  3   4  4


int main(int argc, char *argv[]) {
    bool init = false;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXd d;

    // Load a mesh
    igl::readOBJ(MODEL_PATH "/cube2.obj", V, F);

    igl::opengl::glfw::Viewer viewer;

    const auto update_distance = [&](const int vid) {
        Eigen::VectorXi VS, FS, VT, FT;
        // The selected vertex is the source
        VS.resize(1);
        VS << vid;
        // All vertices are the targets
        VT.setLinSpaced(V.rows(), 0, V.rows() - 1);
        std::cout << "Computing geodesic distance to vertex " << vid << "..." << std::endl;
        igl::exact_geodesic(V, F, VS, FS, VT, FT, d);
        // const double strip_size = 0.05;
        // The function should be 1 on each integer coordinate
        d = (d * 1000).array().eval();
        calc_grip(V, F, d);

        // Compute per-vertex colors
        Eigen::MatrixXd C;
        igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, d, true, C);
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
//                    if (init) { // Only once
//                        return true;
//                    }
                    init = true;
                    int max;
                    bc.maxCoeff(&max);
                    int vid = F(fid, max);
                    update_distance(vid);
                    return true;
                }
                return false;
            };
    viewer.data().set_mesh(V, F);

    std::cout << "Click on mesh to define the source.\n" << std::endl;
    update_distance(0);
    viewer.launch();

    // sample_random_point(V, F, viewer);



    return 0;
}

void sample_random_point(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, igl::opengl::glfw::Viewer &viewer) {
    // Print the vertices and faces matrices
    std::cout << "Vertices: " << std::endl << V << std::endl;
    std::cout << "Faces:    " << std::endl << F << std::endl;

    // Sample 1 random point on mesh
    // Given by barycentric coordinates of the sampled point in face FI
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
    // writeOBJ(MODEL_PATH "/bunny_new.obj", V, F);

}


void calc_grip(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &d) {
    std::cout << "Calc Grip" << std::endl;

    double max_distance = d.maxCoeff();
    double t = max_distance * 0.2;

    // #F by #3 adjacent matrix, the element i,j is the id of the triangle
    // adjacent to the j edge of triangle i
    Eigen::MatrixXi FF;
    // #F by #3 adjacent matrix, the element i,j is the id of edge of the
    // triangle FF(i,j) that is adjacent with triangle i
    Eigen::MatrixXi FFi;
    igl::triangle_triangle_adjacency(F, FF, FFi);

    // max(F) by max(F) cotangent matrix, each row i corresponding to V(i,:)
    Eigen::SparseMatrix<double> Adj;
    igl::adjacency_matrix(F, Adj);

    int num_faces = F.rows();
    Eigen::MatrixXi cut_mask(num_faces, 3);
    std::vector<std::array<int, 2>> cuts;

    std::set<int> visited_faces_ids;
    int current_face_id = 0;
    int e1_id = 0;
    int e2_id = 0;


    // Find first face:
    for (int i = 0; i < num_faces; i++) {
        current_face_id = i;
        std::array<int, 3> fv{F(i, 0), F(i, 1), F(i, 2)};
        std::array<double, 3> vd{d[fv[0]], d[fv[1]], d[fv[2]]};
        std::array<int, 2> e0{F(i, 0), F(i, 1)};
        std::array<int, 2> e1{F(i, 1), F(i, 2)};
        std::array<int, 2> e2{F(i, 2), F(i, 0)};

        bool found_cut = true;
        // We are interested in faces where 2 vertices' distances are closer then t and one is bigger.
        if ((vd[0] < t && vd[1] < t && vd[2] < t) || (vd[0] > t && vd[1] > t && vd[2] > t)) {
            // Not interesting
            found_cut = false;
        }
        if (vd[0] < t && vd[1] < t && vd[2] > t) {
            // We want to add e(v0,v1) = e0
            // e1,e2
            e1_id = 1;
            e2_id = 2;
            cut_mask(i, 0) = 1;
            cuts.push_back(e0);
        } else if (vd[0] > t && vd[1] < t && vd[2] < t) {
            // We want e(v1,v2) = e1
            // e2,20
            e1_id = 2;
            e2_id = 0;
            cuts.push_back(e1);
            cut_mask(i, 1) = 1;
        } else if (vd[0] < t && vd[1] > t && vd[2] < t) {
            // We want e(v0,v2) = e2
            // e0,e1
            e1_id = 0;
            e2_id = 1;
            cuts.push_back(e2);
            cut_mask(i, 2) = 1;
        } else {
            // only vertex
            found_cut = false;
        }
        if (found_cut == true) {
            visited_faces_ids.insert(i);

            std::cout << "=================================" << std::endl;
            std::cout << "face id: " << current_face_id << std::endl;
            std::cout << "v0 id: " << fv[0] << " v0 dist: " << d[fv[0]] << std::endl;
            std::cout << "v1 id: " << fv[1] << " v1 dist: " << d[fv[1]] << std::endl;
            std::cout << "v2 id: " << fv[2] << " v2 dist: " << d[fv[2]] << std::endl;
            std::cout << "e0: (" << fv[0] << "," << fv[1] << ")" << std::endl;
            std::cout << "e1: (" << fv[1] << "," << fv[2] << ")" << std::endl;
            std::cout << "e2: (" << fv[2] << "," << fv[0] << ")" << std::endl;
            std::cout << "=================================" << std::endl;
            break;
        }
    }

    int first_face = current_face_id;

    std::cout << "Found first face: " << current_face_id << std::endl;

    int next_face_id = FF(current_face_id, e1_id);
    int source_edge_id = FFi(current_face_id, e1_id);

    std::cout << "Mooving from face [" << current_face_id << "], edge [" << e1_id << "] to face ["
              << next_face_id << "] edge [" << source_edge_id << "]" << std::endl;

    while (next_face_id != first_face) {
        current_face_id = next_face_id;
        if (visited_faces_ids.find(current_face_id) != visited_faces_ids.end()) {
            // Sanity check
            std::cout << "Been in this face already - is this a bug?" << std::endl;
        }
        visited_faces_ids.insert(current_face_id);


        // find the other edge s.t. d(v1)<t, d(v2) > t
        std::array<int, 3> fv{F(current_face_id, 0), F(current_face_id, 1), F(current_face_id, 2)};
        std::array<double, 3> vd{d[fv[0]], d[fv[1]], d[fv[2]]};

        std::cout << "=================================" << std::endl;
        std::cout << "face id: " << current_face_id << std::endl;
        std::cout << "v0 id: " << fv[0] << " v0 dist: " << d[fv[0]] << std::endl;
        std::cout << "v1 id: " << fv[1] << " v1 dist: " << d[fv[1]] << std::endl;
        std::cout << "v2 id: " << fv[2] << " v2 dist: " << d[fv[2]] << std::endl;
        std::cout << "e0: (" << fv[0] << "," << fv[1] << ")" << std::endl;
        std::cout << "e1: (" << fv[1] << "," << fv[2] << ")" << std::endl;
        std::cout << "e2: (" << fv[2] << "," << fv[0] << ")" << std::endl;
        std::cout << "=================================" << std::endl;

        std::array<int, 2> e0{F(current_face_id, 0), F(current_face_id, 1)};
        std::array<int, 2> e1{F(current_face_id, 1), F(current_face_id, 2)};
        std::array<int, 2> e2{F(current_face_id, 2), F(current_face_id, 0)};
        std::array<std::array<int, 2>, 3> E{e0, e1, e2};

        int closer_v = 0;
        int closer_v_source = get_closer_v_id_from_source(E[source_edge_id], d);
        int other_edge_id = 0;
        int candidate_cut_edge_id = 0;
        bool found_other_edge = false;
        for (int i = 0; i < E.size(); i++) {
            if (source_edge_id == i) {
                continue;
            }
            if (is_other_edge(E[i][0], E[i][1], t, d, closer_v)) {
                other_edge_id = i;
                found_other_edge = true;
                // Get third edge as candidate for cut (will be added if 2 vertices are in wanted dist)
                for (int cand = 0; cand < E.size(); cand++) {
                    if (cand != source_edge_id && cand != other_edge_id) {
                        candidate_cut_edge_id = cand;
                    }
                }
                break;
            }
        }
        std::cout << "Found other edge: " << found_other_edge << std::endl;

        // If add vertex - just go to next face
        // else add edge
        e1_id = source_edge_id;
        e2_id = other_edge_id;
        if (closer_v != closer_v_source) {
            // Add third edge (not source, and not other) -  (closer_v , closer_v_source) to cuts
            cut_mask(current_face_id, candidate_cut_edge_id) = 1;
            switch (candidate_cut_edge_id) {
                case (0):
                    cuts.push_back(e0);
                    break;
                case (1):
                    cuts.push_back(e1);
                    break;
                case (2):
                    cuts.push_back(e2);
                    break;
            }
        } else {
            // "Add vertex" =  just go to next face
        }

        // Go to next face
        next_face_id = FF(current_face_id, other_edge_id);
        source_edge_id = FFi(current_face_id, other_edge_id);
        std::cout << "Mooving from face [" << current_face_id << "], edge [" << other_edge_id << "] to face ["
                  << next_face_id << "] edge [" << source_edge_id << "]" << std::endl;
    }
    for (auto cut: cuts) {
        std::cout << "(" << cut[0] << "," << cut[1] << ") -> ";
    }
}

bool is_other_edge(const int v1, const int v2, const double t, const Eigen::VectorXd &d, int &closer_v) {
    std::cout << "v1 dist: " << d[v1] << " v2 dist: " << d[v2];
    if (d[v1] > t && d[v2] < t) {
        closer_v = v2;
        std::cout << " Found other edge" << std::endl;
        return true;
    } else if (d[v2] > t && d[v1] < t) {
        closer_v = v1;
        std::cout << " Found other edge" << std::endl;
        return true;
    } else {
        std::cout << " Not other edge" << std::endl;
        return false;
    }
}

int get_closer_v_id_from_source(const std::array<int, 2> e, const Eigen::VectorXd &d) {
    double v1 = e[0];
    double v2 = e[1];

    double dist1 = d[e[0]];
    double dist2 = d[e[1]];

    if (d[e[0]] < d[e[1]]) {
        return e[0];
    } else {
        return e[1];
    }
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