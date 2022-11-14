#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/edge_topology.h>
#include <igl/triangle_triangle_adjacency.h>
#include <math.h>
/*** insert any libigl headers here ***/

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

void subdivide_sqrt3(const Eigen::MatrixXd &V,
					 const Eigen::MatrixXi &F,
					 Eigen::MatrixXd &Vout,
					 Eigen::MatrixXi &Fout){
    Eigen::MatrixXd BC;
	igl::barycenter(V, F, BC);
    Vout.resize(V.rows() + BC.rows(), 3);
    Vout << V, BC;
    // cout << V.rows() << " " << BC.rows() << " " << Vout.rows() << endl;
    Eigen::MatrixXi EV, EF, FE;
    igl::edge_topology(V, F, EV, FE, EF);
    // step 2
    igl::adjacency_list(F, VV);
    for (int i = 0; i < V.rows(); ++i)
    {
        int n = VV[i].size();
        double an = (4 - 2 * std::cos(2 * M_PI / n)) / 9;
        Eigen::RowVector3d sum;
        sum << 0.0, 0.0, 0.0;
        for (int t = 0; t < n; ++t)
        {
            sum += V.row(VV[i][t]);
        }
        Vout.row(i) << (1 - an) * V.row(i) + (an / n) * sum;
    }
    // step 3
    int row_id = 0;
    int V_row = V.rows();
    int Fout_row = F.rows() * 3;
    Fout.setZero(Fout_row, 3);
    for (int i = 0; i < EV.rows(); ++i)
    {
        // if (EF(i, 0) == -1 || EF(i, 1) == -1) // not EV
        //     continue;
        if (EF(i, 0) == -1)
        {
            Fout.row(row_id) << EF(i, 1) + V_row, EV(i, 1), EV(i, 0);
            row_id++;
            continue;
        }
        if (EF(i, 1) == -1)
        {
            Fout.row(row_id) << EV(i, 0), EV(i, 1), EF(i, 0) + V_row;
            row_id++;
            continue;
        }
        Fout.row(row_id) << EV(i, 0), EF(i, 1) + V_row, EF(i, 0) + V_row; // order matters
        row_id++;
        Fout.row(row_id) << EF(i, 1) + V_row, EV(i, 1), EF(i, 0) + V_row;
        row_id++;
    }
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
        cout << "Vertex   Connect to Face" << endl;
        for (int i = 0; i < VF.size(); ++i)
        {
            cout << i << "        ";
            for (int j = 0; j < VF[i].size(); ++j)
            {
                cout << VF[i][j] << " ";
            }
            cout << endl;
        }
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
        igl::adjacency_list(F, VV);
        cout << "Vertex   Connect to Vertex" << endl;
        for (int i = 0; i < VV.size(); ++i)
        {
            cout << i << "        ";
            for (int j = 0; j < VV[i].size(); ++j)
            {
                cout << VV[i][j] << " ";
            }
            cout << endl;
        }
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
        igl::per_face_normals(V, F, FN);
        // Set the viewer normals.
        viewer.data().set_normals(FN);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-vertex normals here: store in VN.
        igl::per_vertex_normals(V, F, VN);
        // Set the viewer normals.
        viewer.data().set_normals(VN);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
        float thres = 90;
        igl::per_corner_normals(V, F, thres, CN);
        //Set the viewer normals
        viewer.data().set_normals(CN);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
        igl::facet_components(F, cid);
        int component_num = cid.maxCoeff() + 1;
        std::vector<int> face_num(component_num, 0);
        for (int i = 0; i < cid.rows(); ++i)
        {
            face_num[cid(i)]++;
        }
        cout << "component_num:  " << component_num << endl;
        cout << "face num:  " << endl;
        for (int i = 0; i < face_num.size(); ++i)
        {
            cout << face_num[i] << "  ";
        }

        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
        igl::jet(cid, 1, component_colors_per_face);

        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7') {
		Eigen::MatrixXd Vout;
		Eigen::MatrixXi Fout;
        // Fill the subdivide_sqrt3() function with your code for sqrt(3) subdivision.
		subdivide_sqrt3(V,F,Vout,Fout);
        cout << "out" << endl;
        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOFF(filename,V,F);
  viewer.data().clear();
  viewer.data().set_mesh(V,F);
  viewer.data().compute_normals();
  viewer.core().align_camera_center(V, F);
  return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]); // Mesh provided as command line argument
    }
    else {
        filename = std::string("../data/bunny.off"); // Default mesh
    }
	
    load_mesh(viewer,filename,V,F);

    // callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    viewer.launch();
}
