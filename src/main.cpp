//
// Created by Verose on 12/14/2018.
//
#include <igl/triangle_triangle_adjacency.h>
#include <igl/adjacency_matrix.h>
#include <igl/HalfEdgeIterator.h>
#include <igl/combine.h>

#include <iostream>
#include <fstream>
#include <array>
#include <unordered_map>
#include <set>

#include "main.h"
#include "model_path.h"

//VT = v1 v2 v3
//FT = f1
//d  = 1  3   4  4

#ifndef DEBUG
#define DEBUG false
#endif


int main(int argc, char *argv[]) {
//    bool init = false;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXd d;

    // Load a mesh
    igl::readOBJ(MODEL_PATH "/Cow.obj", V, F);

    igl::opengl::glfw::Viewer viewer;

    const auto show_holder_with_distances = [&](const int vid) {
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

        // Generating the mesh without using the calculated cut.
        double max_distance = d.maxCoeff();
        double t = max_distance * 0.5;
        save_grip_mesh(V, F, d, t);

        std::vector<Eigen::Vector2i> cuts;
        calc_grip(V, F, d, t, cuts);

        display_cut(V, F, d, viewer, cuts);
    };

    // Plot a distance when a vertex is picked
    viewer.callback_mouse_down =
            [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool {
                int fid;
                Eigen::Vector3f bc;
                // Cast a ray in the view direction starting from the mouse position
                double x = viewer.current_mouse_x;
                double y = viewer.core.viewport(3) - viewer.current_mouse_y;

                if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view, viewer.core.proj,
                                             viewer.core.viewport, V, F, fid, bc)) {
                    int max;
                    bc.maxCoeff(&max);
                    int vid = F(fid, max);
                    viewer.data().clear();
                    show_holder_with_distances(vid);

                    std::cout << "Press key '1' to finalize the mesh or pick a new point" << std::endl;
                }
                return false;
            };
    viewer.callback_key_down =
            [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int) -> bool {
                if (key == '1') {
                    viewer.data().clear();
                    move_gripper_in_normal_direction();
                    invert_gripper_normal_direction();
                    combine_meshes();

                    Eigen::MatrixXd V_holder;
                    Eigen::MatrixXi F_holder;
                    igl::readOBJ(MODEL_PATH "/holder.obj", V_holder, F_holder);

                    viewer.data().set_mesh(V_holder, F_holder);

                    std::cout << "Your new cut is ready!" << std::endl;
                    std::cout << "To reset press key '2'" << std::endl;
                }
                else if (key == '2') {
                    viewer.data().clear();
                    viewer.data().set_mesh(V, F);
                }
                return false;
            };

    viewer.data().set_mesh(V, F);
    std::cout << "Click on mesh to define source for holder\n" << std::endl;
    viewer.launch();

    return 0;
}

void display_cut(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd d,
        igl::opengl::glfw::Viewer &viewer, const std::vector<Eigen::Vector2i> &cuts) {
    for (const auto &cut : cuts) {
        double point_index_start = cut[0];
        double point_index_end = cut[1];
        Eigen::RowVector3d point_start(V(point_index_start, 0), V(point_index_start, 1), V(point_index_start, 2));
        Eigen::RowVector3d point_end(V(point_index_end, 0), V(point_index_end, 1), V(point_index_end, 2));
        viewer.data().add_points(point_start, Eigen::RowVector3d(1, 0, 0));  // show the first point of each edge
        viewer.data().add_edges(point_start, point_end, Eigen::RowVector3d(1, 0, 0));  // print edge (it's printed but not visible)
    }

    // Compute per-vertex colors
    Eigen::MatrixXd C;
    igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, d, true, C);

    // Plot the mesh
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
}

void save_grip_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &d, const double t) {
    std::unordered_map<int, int> old_to_new_vertex_id_map;
    std::set<int> in_vertices_ids;

    std::ofstream grip_mesh;
    grip_mesh.open(MODEL_PATH "/grip_mesh.obj");
    grip_mesh << "# grip_mesh \n";

    int num_vertices = 0;
    for (int i = 0; i < V.rows(); i++) {
        if (d[i] < t) {
            // v v0 v1 v2
            grip_mesh << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
            grip_mesh.flush();
            num_vertices++;
            old_to_new_vertex_id_map[i] = num_vertices;
            in_vertices_ids.insert(i);
        }
    }
    for (int i = 0; i < F.rows(); i++) {
        Eigen::Vector3i fv{F(i, 0), F(i, 1), F(i, 2)};
        bool save_face =
                in_vertices_ids.find(fv[0]) != in_vertices_ids.end() &&
                in_vertices_ids.find(fv[1]) != in_vertices_ids.end() &&
                in_vertices_ids.find(fv[2]) != in_vertices_ids.end();

        if (save_face) {
            int v0 = old_to_new_vertex_id_map[fv[0]];
            int v1 = old_to_new_vertex_id_map[fv[1]];
            int v2 = old_to_new_vertex_id_map[fv[2]];

            // f 45//45 193//193 281//281
            grip_mesh << "f "
                      << v0 << "//" << v0 << " "
                      << v1 << "//" << v1 << " "
                      << v2 << "//" << v2 << "\n";
            grip_mesh.flush();
        }
    }
    grip_mesh.close();
}

void move_gripper_in_normal_direction() {
    Eigen::MatrixXd V_grip;
    Eigen::MatrixXi F_grip;
    igl::readOBJ(MODEL_PATH "/grip_mesh.obj", V_grip, F_grip);

    Eigen::MatrixXd N;

    // Load a mesh
    igl::per_vertex_normals(
            V_grip,
            F_grip,
            igl::PerVertexNormalsWeightingType::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA,
            N);
    double stride = 0.02;
    Eigen::MatrixXd result(V_grip.rows(), 3);

    for (int i = 0; i < V_grip.rows(); i++) {
        Eigen::RowVector3d point(V_grip(i, 0), V_grip(i, 1), V_grip(i, 2));
        Eigen::RowVector3d normal(N(i, 0), N(i, 1), N(i, 2));
        Eigen::RowVector3d sub = normal - point;
        sub.normalize();
        Eigen::RowVector3d res = point + stride*sub;
        result(i, 0) = res(0);
        result(i, 1) = res(1);
        result(i, 2) = res(2);
    }

    igl::writeOBJ(MODEL_PATH "/grip_mesh_out.obj", result, F_grip);
}

void invert_gripper_normal_direction() {
    Eigen::MatrixXd V_grip;
    Eigen::MatrixXi F_grip;
    igl::readOBJ(MODEL_PATH "/grip_mesh.obj", V_grip, F_grip);
    Eigen::MatrixXi F_grip_inv(F_grip.rows(), 3);

    Eigen::MatrixXi N;

    for (int i = 0; i < F_grip.rows(); i++) {
        F_grip_inv(i, 0) = F_grip(i, 2);
        F_grip_inv(i, 2) = F_grip(i, 0);
        F_grip_inv(i, 1) = F_grip(i, 1);
    }

    igl::writeOBJ(MODEL_PATH "/grip_mesh_in.obj", V_grip, F_grip_inv);
}

void combine_meshes() {
    Eigen::MatrixXd V_grip_in, V_grip_out, V;
    Eigen::MatrixXi F_grip_in, F_grip_out, F;
    igl::readOBJ(MODEL_PATH "/grip_mesh_in.obj", V_grip_in, F_grip_in);
    igl::readOBJ(MODEL_PATH "/grip_mesh_out.obj", V_grip_out, F_grip_out);
    std::vector<Eigen::MatrixXd> V_list{V_grip_in, V_grip_out};
    std::vector<Eigen::MatrixXi> F_list{F_grip_in, F_grip_out};

    igl::combine(V_list, F_list, V, F);

    igl::writeOBJ(MODEL_PATH "/holder.obj", V, F);
}

void calc_grip(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &d, const double t,
               std::vector<Eigen::Vector2i> &cuts) {
    std::cout << "Calc Grip" << std::endl;

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

    std::set<int> visited_faces_ids;
    int current_face_id = 0;
    int e1_id = 0;
    int e2_id = 0;

    // Find first face:
    for (int i = 0; i < num_faces; i++) {
        current_face_id = i;
        Eigen::Vector3i fv{F(i, 0), F(i, 1), F(i, 2)};
        Eigen::Vector3d vd{d[fv[0]], d[fv[1]], d[fv[2]]};
        Eigen::Vector2i e0{F(i, 0), F(i, 1)};
        Eigen::Vector2i e1{F(i, 1), F(i, 2)};
        Eigen::Vector2i e2{F(i, 2), F(i, 0)};

        bool found_cut = true;
        // We are interested in faces where 2 vertices' distances are closer than t and one is bigger.
        if ((vd[0] < t && vd[1] < t && vd[2] < t) ||
            (vd[0] > t && vd[1] > t && vd[2] > t)) {
            // Not interesting
            found_cut = false;
        }
        if (vd[0] < t &&
            vd[1] < t &&
            vd[2] > t) {
            // We want to add e(v0,v1) = e0
            // e1,e2
            e1_id = 1;
            e2_id = 2;
            cut_mask(i, 0) = 1;
            cuts.push_back(e0);
        } else if (vd[0] > t &&
                   vd[1] < t &&
                   vd[2] < t) {
            // We want e(v1,v2) = e1
            // e2,e0
            e1_id = 2;
            e2_id = 0;
            cuts.push_back(e1);
            cut_mask(i, 1) = 1;
        } else if (vd[0] < t &&
                   vd[1] > t &&
                   vd[2] < t) {
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
        if (found_cut) {
            visited_faces_ids.insert(i);

#if DEBUG
            std::cout << "=================================" << std::endl;
            std::cout << "face id: " << current_face_id << std::endl;
            std::cout << "v0 id: " << fv[0] << " v0 dist: " << d[fv[0]] << std::endl;
            std::cout << "v1 id: " << fv[1] << " v1 dist: " << d[fv[1]] << std::endl;
            std::cout << "v2 id: " << fv[2] << " v2 dist: " << d[fv[2]] << std::endl;
            std::cout << "e0: (" << fv[0] << "," << fv[1] << ")" << std::endl;
            std::cout << "e1: (" << fv[1] << "," << fv[2] << ")" << std::endl;
            std::cout << "e2: (" << fv[2] << "," << fv[0] << ")" << std::endl;
            std::cout << "=================================" << std::endl;
#endif
            break;
        }
    }

    int first_face = current_face_id;
    int next_face_id = FF(current_face_id, e1_id);
    int source_edge_id = FFi(current_face_id, e1_id);

#if DEBUG
    std::cout << "Found first face: " << current_face_id << std::endl;
    std::cout << "Moving from face [" << current_face_id << "], edge [" << e1_id << "] to face ["
              << next_face_id << "] edge [" << source_edge_id << "]" << std::endl;
#endif

    while (next_face_id != first_face) {
        current_face_id = next_face_id;
        if (visited_faces_ids.find(current_face_id) != visited_faces_ids.end()) {
            // Sanity check
            std::cout << "Been in this face already - is this a bug?" << std::endl;
        }
        visited_faces_ids.insert(current_face_id);


        // find the other edge s.t. d(v1)<t, d(v2) > t
        Eigen::Vector3i fv{F(current_face_id, 0), F(current_face_id, 1), F(current_face_id, 2)};
        Eigen::Vector3d vd{d[fv[0]], d[fv[1]], d[fv[2]]};
#if DEBUG
        std::cout << "=================================" << std::endl;
        std::cout << "face id: " << current_face_id << std::endl;
        std::cout << "v0 id: " << fv[0] << " v0 dist: " << d[fv[0]] << std::endl;
        std::cout << "v1 id: " << fv[1] << " v1 dist: " << d[fv[1]] << std::endl;
        std::cout << "v2 id: " << fv[2] << " v2 dist: " << d[fv[2]] << std::endl;
        std::cout << "e0: (" << fv[0] << "," << fv[1] << ")" << std::endl;
        std::cout << "e1: (" << fv[1] << "," << fv[2] << ")" << std::endl;
        std::cout << "e2: (" << fv[2] << "," << fv[0] << ")" << std::endl;
        std::cout << "=================================" << std::endl;
#endif
        Eigen::Vector2i e0{F(current_face_id, 0), F(current_face_id, 1)};
        Eigen::Vector2i e1{F(current_face_id, 1), F(current_face_id, 2)};
        Eigen::Vector2i e2{F(current_face_id, 2), F(current_face_id, 0)};
        std::array<Eigen::Vector2i, 3> E{e0, e1, e2};

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
#if DEBUG
        std::cout << "Found other edge: " << found_other_edge << std::endl;
#endif
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
#if DEBUG
        std::cout << "Moving from face [" << current_face_id << "], edge [" << other_edge_id << "] to face ["
                  << next_face_id << "] edge [" << source_edge_id << "]" << std::endl;
#endif
    }
#if DEBUG
    for (auto cut: cuts) {
        std::cout << "(" << cut[0] << "," << cut[1] << ") -> ";
    }
    std::cout << std::endl;
#endif
}

bool is_other_edge(const int v1, const int v2, const double t, const Eigen::VectorXd &d, int &closer_v) {
#if DEBUG
    std::cout << "v1 dist: " << d[v1] << " v2 dist: " << d[v2];
#endif
    if (d[v1] > t && d[v2] < t) {
        closer_v = v2;
#if DEBUG
        std::cout << " Found other edge" << std::endl;
#endif
        return true;
    } else if (d[v2] > t && d[v1] < t) {
        closer_v = v1;
#if DEBUG
        std::cout << " Found other edge" << std::endl;
#endif
        return true;
    } else {
#if DEBUG
        std::cout << " Not other edge" << std::endl;
#endif
        return false;
    }
}

int get_closer_v_id_from_source(const Eigen::Vector2i e, const Eigen::VectorXd &d) {
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
