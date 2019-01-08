//
// Created by Verose on 12/14/2018.
//
#include <igl/triangle_triangle_adjacency.h>
#include <igl/adjacency_matrix.h>
#include <igl/HalfEdgeIterator.h>
#include <igl/combine.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <iostream>
#include <fstream>
#include <array>
#include <unordered_map>
#include <set>

#include "main.h"

//VT = v1 v2 v3
//FT = f1
//d  = 1  3   4  4

#ifndef DEBUG
#define DEBUG false
#endif

#ifndef MODEL_PATH
#define MODEL_PATH "../models"
#endif

#ifndef DIST_MULTIPLAYER
#define DIST_MULTIPLAYER 1000
#endif


int main(int argc, char *argv[]) {
//    bool init = false;
    Eigen::MatrixXd V, V_in, V_out, V_holder;
    Eigen::MatrixXi F, F_in, F_out, F_holder;
    Eigen::VectorXd d;
    std::string save_path = MODEL_PATH "/holder.obj";
    std::string load_path = MODEL_PATH "/Cow.obj";

    // Load a mesh
    igl::readOBJ(load_path, V, F);

    igl::opengl::glfw::Viewer viewer;

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    // Customize the menu
    double distance_mod = 0.5;
    double holder_width = 0.02;
    bool cut_selected = false;

    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();
    };

    // Draw additional windows
    menu.callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(250, 195), ImGuiSetCond_FirstUseEver);
        ImGui::Begin("Holder Parameters", nullptr, ImGuiWindowFlags_NoSavedSettings);

        // Expose the variables directly
        ImGui::PushItemWidth(-100);
        ImGui::InputDouble("Max Distance %", &distance_mod, 0.05, 0.1, "%.2f");
        ImGui::InputDouble("Holder Width", &holder_width, 0.005, 0.01, "%.3f");
        ImGui::InputText("Load Holder Path", load_path);
        ImGui::InputText("Save Holder Path", save_path);
        ImGui::PopItemWidth();

        ImGui::Text("(1) Cut");
        ImGui::Text("(2) Clear");
        ImGui::Text("(3) Load");
        ImGui::Text("(4) Save");

        ImGui::End();
    };

    const auto show_holder_with_distances = [&](const int vid) {
        Eigen::VectorXi VS, FS, VT, FT;
        // The selected vertex is the source
        VS.resize(1);
        VS << vid;
        // All vertices are the targets
        VT.setLinSpaced(V.rows(), 0, V.rows() - 1);
        igl::exact_geodesic(V, F, VS, FS, VT, FT, d);
        d = (d * DIST_MULTIPLAYER).array().eval();

        // Generating the mesh without using the calculated cut.
        double max_distance = d.maxCoeff();
        double t = max_distance * distance_mod;
        save_grip_mesh(V, F, d, t);

        std::vector<Eigen::Vector2i> cuts;
        calc_grip(V, F, d, t, cuts);

        display_cut(V, F, d, viewer, cuts);
    };

    // Plot a distance when a vertex is picked
    viewer.callback_mouse_down =
            [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool {
                int fid;

                // Don't allow to pick new source if cut is selected
                if (cut_selected) {
                    return false;
                }

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

                    Eigen::MatrixXi F_grip_out;
                    Eigen::MatrixXd V_grip;
                    igl::readOBJ(MODEL_PATH "/grip_mesh.obj", V_grip, F_grip_out);
                    Eigen::MatrixXd V_grip_out(V_grip.rows(), 3);
                    move_gripper_in_normal_direction(holder_width, V_grip, F_grip_out, V_grip_out);

                    Eigen::MatrixXi F_grip_in(F_grip_out.rows(), 3);
                    invert_gripper_normal_direction(V_grip, F_grip_out, F_grip_in);
                    combine_meshes(V_grip, F_grip_in, V_grip_out, F_grip_out, V_holder, F_holder);
                    viewer.data().set_mesh(V_holder, F_holder);

                    std::cout << "Your new cut is ready!" << std::endl;
                    std::cout << "To reset press '2', to load new model press '3', to save press '4'" << std::endl;
                    cut_selected = true;
                }
                else if (key == '2') {
                    viewer.data().clear();
                    viewer.data().set_mesh(V, F);
                    std::cout << "Resetting to original mesh" << std::endl;
                    std::cout << "To choose new source for holder press '1', to load new model press '3', to save press '4'" << std::endl;
                    cut_selected = false;
                }
                else if (key == '3') {
                    igl::readOBJ(load_path, V, F);
                    viewer.data().clear();
                    viewer.data().set_mesh(V, F);
                    std::cout << "Loaded " << load_path.c_str() << std::endl;
                    std::cout << "To choose new source for holder press '1', to reset press '2', to save press '4'" << std::endl;
                    cut_selected = false;
                }
                else if (key == '4') {
                    igl::writeOBJ(save_path, V_holder, F_holder);
                    std::cout << "Saved " << save_path.c_str() << std::endl;
                    std::cout << "To choose new source for holder press '1', to reset press '2', to load new model press '3'" << std::endl;
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
        long point_index_start = cut[0];
        long point_index_end = cut[1];
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

void move_gripper_in_normal_direction(double holder_width, Eigen::MatrixXd &V_mesh, Eigen::MatrixXi &F_grip_out,
        Eigen::MatrixXd &V_grip_out) {
    Eigen::MatrixXd N;

    // Load a mesh
    igl::per_vertex_normals(
            V_mesh,
            F_grip_out,
            igl::PerVertexNormalsWeightingType::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA,
            N);

    for (int i = 0; i < V_mesh.rows(); i++) {
        Eigen::RowVector3d point(V_mesh(i, 0), V_mesh(i, 1), V_mesh(i, 2));
        Eigen::RowVector3d normal(N(i, 0), N(i, 1), N(i, 2));
        Eigen::RowVector3d sub = normal - point;
        sub.normalize();
        Eigen::RowVector3d res = point + holder_width*sub;
        V_grip_out(i, 0) = res(0);
        V_grip_out(i, 1) = res(1);
        V_grip_out(i, 2) = res(2);
    }
}

void invert_gripper_normal_direction(Eigen::MatrixXd &V_grip, Eigen::MatrixXi &F_grip, Eigen::MatrixXi &F_grip_inv) {
    for (int i = 0; i < F_grip.rows(); i++) {
        F_grip_inv(i, 0) = F_grip(i, 2);
        F_grip_inv(i, 2) = F_grip(i, 0);
        F_grip_inv(i, 1) = F_grip(i, 1);
    }
}

void combine_meshes(Eigen::MatrixXd &V_grip_in, Eigen::MatrixXi &F_grip_in,
        Eigen::MatrixXd &V_grip_out, Eigen::MatrixXi &F_grip_out,
        Eigen::MatrixXd &V_holder, Eigen::MatrixXi &F_holder) {
    std::vector<Eigen::MatrixXd> V_list{V_grip_in, V_grip_out};
    std::vector<Eigen::MatrixXi> F_list{F_grip_in, F_grip_out};

    igl::combine(V_list, F_list, V_holder, F_holder);
}

void calc_grip(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &d, const double t,
               std::vector<Eigen::Vector2i> &cuts) {
    // #F by #3 adjacent matrix, the element i,j is the id of the triangle
    // adjacent to the j edge of triangle i
    Eigen::MatrixXi FF;
    // #F by #3 adjacent matrix, the element i,j is the id of edge of the
    // triangle FF(i,j) that is adjacent with triangle i
    Eigen::MatrixXi FFi;
    igl::triangle_triangle_adjacency(F, FF, FFi);

    long num_faces = F.rows();

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
        if (vd[0] < t && vd[1] < t && vd[2] > t) {
            // We want to add e(v0,v1) = e0
            // e1,e2
            e1_id = 1;
            e2_id = 2;
            cuts.push_back(e0);
        } else if (vd[0] > t && vd[1] < t && vd[2] < t) {
            // We want e(v1,v2) = e1
            // e2,e0
            e1_id = 2;
            e2_id = 0;
            cuts.push_back(e1);
        } else if (vd[0] < t && vd[1] > t && vd[2] < t) {
            // We want e(v0,v2) = e2
            // e0,e1
            e1_id = 0;
            e2_id = 1;
            cuts.push_back(e2);
        } else {
            // only vertex
            found_cut = false;
        }
        if (found_cut) {
            visited_faces_ids.insert(i);
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
            std::cout << "Sanity check - Been in this face already." << std::endl;
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
        int closer_v_source = get_edge_closer_v_id(E[source_edge_id], d);
        int other_edge_id = 0;
        int candidate_cut_edge_id = 0;

        for (int i = 0; i < E.size(); i++) {
            if (source_edge_id == i) {
                continue;
            }
            if (is_other_edge(E[i][0], E[i][1], t, d)) {
                closer_v = get_edge_closer_v_id(E[i], d);
                other_edge_id = i;
                // Get third edge as candidate for cut (will be added if 2 vertices are in wanted dist)
                for (int cand = 0; cand < E.size(); cand++) {
                    if (cand != source_edge_id && cand != other_edge_id) {
                        candidate_cut_edge_id = cand;
                    }
                }
                break;
            }
        }

        if (closer_v != closer_v_source) {
            // Add third edge (not source, and not other) -  i.e. (closer_v , closer_v_source) to cuts
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
            // "Add vertex",  just go to next face
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

bool is_other_edge(const int v1, const int v2, const double t, const Eigen::VectorXd &d) {
    if ((d[v1] > t && d[v2] < t) || (d[v2] > t && d[v1] < t)) {
        return true;
    } else {
        return false;
    }
}

int get_edge_closer_v_id(const Eigen::Vector2i e, const Eigen::VectorXd &d) {
    if (d[e[0]] < d[e[1]]) {
        return e[0];
    } else {
        return e[1];
    }
}

int get_edge_further_v_id(const Eigen::Vector2i e, const Eigen::VectorXd &d) {
    return (get_edge_closer_v_id(e,d) == e[0] ? e[1] : e[0]);
}


void find_new_p_on_edge(const int closer_v_id, const double step_size, const int further_v_id, const Eigen::MatrixXd &V,
                        Eigen::RowVector3d &new_point) {
    Eigen::RowVector3d closer(V(closer_v_id, 0), V(closer_v_id, 1), V(closer_v_id, 2));
    Eigen::RowVector3d further(V(further_v_id, 0), V(further_v_id, 1), V(further_v_id, 2));
    Eigen::RowVector3d direction = further - closer;
    direction.normalize();
    Eigen::RowVector3d res = closer + 0.00001 * direction;
    new_point(0) = res(0);
    new_point(1) = res(1);
    new_point(2) = res(2);
}

void calc_smooth_gripper(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &d, const double t,
                         Eigen::MatrixXd &V_for_smooth, Eigen::MatrixXi &F_for_smooth) {
    std::unordered_map<int, Eigen::Vector3d> old_v_id_to_new_v;
    std::set<int> altered_vertices_ids;
    // #F by #3 adjacent matrix, the element i,j is the id of the triangle
    // adjacent to the j edge of triangle i
    Eigen::MatrixXi FF;
    // #F by #3 adjacent matrix, the element i,j is the id of edge of the
    // triangle FF(i,j) that is adjacent with triangle i
    Eigen::MatrixXi FFi;
    igl::triangle_triangle_adjacency(F, FF, FFi);


    // Find first face:
    int num_faces = F.rows();
    int current_face_id = 0;
    int e1_id = 0;
    int e2_id = 0;
    std::set<int> visited_faces_ids;
    for (int i = 0; i < num_faces; i++) {
        current_face_id = i;
        Eigen::Vector3i fv{F(i, 0), F(i, 1), F(i, 2)};
        Eigen::Vector3d vd{d[fv[0]], d[fv[1]], d[fv[2]]};
        Eigen::Vector2i e0{F(i, 0), F(i, 1)};
        Eigen::Vector2i e1{F(i, 1), F(i, 2)};
        Eigen::Vector2i e2{F(i, 2), F(i, 0)};

        bool found_cut = true;
        // We are interested in faces where 2 vertices' distances are closer than t and one is bigger to start with.
        if (vd[0] < t && vd[1] < t && vd[2] > t) {
            // We want to transform v2 to be within t
            e1_id = 1;
            e2_id = 2;
        } else if (vd[0] > t && vd[1] < t && vd[2] < t) {
            // We want to transform v0 to be within t
            e1_id = 2;
            e2_id = 0;
        } else if (vd[0] < t && vd[1] > t && vd[2] < t) {
            // We want to transform v1 to be within t
            e1_id = 0;
            e2_id = 1;
        } else {
            found_cut = false;
        }
        if (found_cut) {
            break;
        }
    }

    int first_face = current_face_id;
    visited_faces_ids.insert(current_face_id);
    int next_face_id = FF(current_face_id, e1_id); // Go to the neighbor face using an interesting edge.
    int source_edge_id = FFi(current_face_id, e1_id); // Remember which edge is the one we used for moving.

#if DEBUG
    std::cout << "Moving from face [" << current_face_id << "], edge [" << other_edge_id << "] to face ["
                  << next_face_id << "] edge [" << source_edge_id << "]" << std::endl;
#endif

    while (next_face_id != first_face) {
        current_face_id = next_face_id;
        if (visited_faces_ids.find(current_face_id) != visited_faces_ids.end()) {
            // Sanity check.
            std::cout << "Been in this face already." << std::endl;
        }
        visited_faces_ids.insert(current_face_id);

        // find the other edge s.t. d(v1)<t, d(v2)>t
        Eigen::Vector3i fv{F(current_face_id, 0), F(current_face_id, 1), F(current_face_id, 2)};
        Eigen::Vector3d vd{d[fv[0]], d[fv[1]], d[fv[2]]};
        Eigen::Vector2i e0{F(current_face_id, 0), F(current_face_id, 1)};
        Eigen::Vector2i e1{F(current_face_id, 1), F(current_face_id, 2)};
        Eigen::Vector2i e2{F(current_face_id, 2), F(current_face_id, 0)};
        std::array<Eigen::Vector2i, 3> E{e0, e1, e2};

        // Find the other edge s.t. one vertex is within t and one is out of t.
        // O - Vertex distance is smaller than t
        // X - Vertex distance is bigger than t

        //     O            X
        //    / \          / \
        //   /   \        /   \
        //  /     \      /     \
        // X-------X    O-------O
        //     1            2

        // ID of the vertex that was within t on the source edge.
        int closer_v_id_source = get_edge_closer_v_id(E[source_edge_id], d);
        int further_v_id_source = get_edge_further_v_id(E[source_edge_id], d);

        // Will hold the ID of the vertex that is within t on the other edge.
        // For face 1 in the example, it will be the same as closer_v_source, for face 2 it will be a different vertex.
        int other_edge_id = 0;
        for (int i = 0; i < E.size(); i++) {
            if (source_edge_id == i) {
                // This is the edge we started from, we already have the coordinate on it.
                continue;
            }
            if (is_other_edge(E[i][0], E[i][1], t, d)) {
                // This is the other edge that has one vertex in and one out.
                other_edge_id = i;
                break;
            }
        }
        int closer_v_id_other_edge = get_edge_closer_v_id(E[other_edge_id], d);
        int further_v_id_other_edge = get_edge_further_v_id(E[other_edge_id], d);

        if (closer_v_id_other_edge == closer_v_id_source) {
            //
            //     O              O
            //    / \     ->     / \
            //   /   \          /   \
            //  /     \        O-----O
            // X-------X
            //
            double step_size = (t - d[closer_v_id_source]) / DIST_MULTIPLAYER;

            Eigen::RowVector3d p1;
            find_new_p_on_edge(closer_v_id_source, step_size, further_v_id_source, V, p1);
            old_v_id_to_new_v[further_v_id_source] = p1;
            altered_vertices_ids.insert(further_v_id_source);

            Eigen::RowVector3d p2;
            find_new_p_on_edge(closer_v_id_other_edge, step_size, further_v_id_other_edge, V, p2);
            old_v_id_to_new_v[further_v_id_other_edge] = p2;
            altered_vertices_ids.insert(further_v_id_other_edge);


        } else {
            //
            //     X
            //    / \     ->      O
            //   /   \           / \
            //  /     \         /   \
            // O-------O       O-----O
            //
            // Will be taken care of using the first case.
        }

        // Go to next face
        next_face_id = FF(current_face_id, other_edge_id);
        source_edge_id = FFi(current_face_id, other_edge_id);
        std::cout << "Moving from face [" << current_face_id << "], edge [" << other_edge_id << "] to face ["
                  << next_face_id << "] edge [" << source_edge_id << "]" << std::endl;
    }
    std::cout << std::endl;
    gen_smooth_gripper_mesh(V, F,old_v_id_to_new_v, altered_vertices_ids, V_for_smooth, F_for_smooth);
}

void gen_smooth_gripper_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                std::unordered_map<int, Eigen::Vector3d> &old_v_id_to_new_v,
                                const std::set<int> &altered_vertices_ids, Eigen::MatrixXd &V_for_smooth,
                                Eigen::MatrixXi &F_for_smooth) {
    V_for_smooth = V;
    F_for_smooth = F;
    for (int i = 0; i < V_for_smooth.rows(); i++) {
        if (altered_vertices_ids.find(i) != altered_vertices_ids.end()) {
            V_for_smooth(i, 0) = old_v_id_to_new_v[i](0);
            V_for_smooth(i, 2) = old_v_id_to_new_v[i](1);
            V_for_smooth(i, 1) = old_v_id_to_new_v[i](2);
        }
    }
}