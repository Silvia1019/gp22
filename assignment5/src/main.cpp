#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>

#include "Lasso.h"
#include "Colors.h"

#include <igl/slice.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_face_normals.h>
#include <igl/grad.h>

//activate this for alternate UI (easier to debug)
// #define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

//vertex array, #V x3
Eigen::MatrixXd V(0,3), V_cp(0, 3);
//face array, #F x3
Eigen::MatrixXi F(0,3);

//mouse interaction
enum MouseMode { SELECT, TRANSLATE, ROTATE, NONE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;

//for selecting vertices
std::unique_ptr<Lasso> lasso;
//list of currently selected vertices
Eigen::VectorXi selected_v(0,1);

//for saving constrained vertices
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0,1);
//list of all vertices belonging to handles, #HV x1
Eigen::VectorXi handle_vertices(0,1);
//centroids of handle regions, #H x1
Eigen::MatrixXd handle_centroids(0,3);
//updated positions of handle vertices, #HV x3
Eigen::MatrixXd handle_vertex_positions(0,3);
//index of handle being moved
int moving_handle = -1;
//rotation and translation for the handle being moved
Eigen::Vector3f translation(0,0,0);
Eigen::Vector4f rotation(0,0,0,1.);
typedef Eigen::Triplet<double> T;
//per vertex color array, #V x3
Eigen::MatrixXd vertex_colors;

Eigen::SparseMatrix<double> L, M, Minv, A, A_ff, A_fc;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver;
Eigen::MatrixXd V_c, V_f, R;
Eigen::MatrixXd B, B_p, S_p;
Eigen::VectorXi free_vertices, edge_index;
Eigen::MatrixXd d;
int mesh_shape = 3; // S, B, B_p, S_p

//When false, use standard displacement vectors for details, when true use Deformation Transfer from part 2
bool use_deformation_transfer = true;

//function declarations (see below for implementation)
bool solve(Viewer& viewer);
void get_new_handle_locations();
Eigen::Vector3f computeTranslation (Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
Eigen::Vector4f computeRotation(Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
void compute_handle_centroids();
Eigen::MatrixXd readMatrix(const char *filename);

bool callback_mouse_down(Viewer& viewer, int button, int modifier);
bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y);
bool callback_mouse_up(Viewer& viewer, int button, int modifier);
bool callback_pre_draw(Viewer& viewer);
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);
void onNewHandleID();
void applySelection();

void get_smooth()
{
  // cout << "get_smooth" << endl;
  // 1.2
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  igl::invert_diag(M,Minv);
  A = L * Minv * L;
	igl::slice(A, free_vertices, free_vertices, A_ff);
	igl::slice(A, free_vertices, handle_vertices, A_fc);
  igl::slice(V, handle_vertices, 1, V_c);
  // cout << A_ff << endl;
  // cout << V_c.rows() << "   " << V_c.cols() << endl;
  // cout << V_c << endl;
  solver.compute(A_ff);
  R = -A_fc * V_c;
  V_f = solver.solve(R);
	
  B.setZero(V.rows(), 3);
	B = V;
	igl::slice_into(V_f, free_vertices, 1, B);
}

void encode_displacement()
{
  // cout << "encode_displacement" << endl;
  edge_index.setZero(V.rows());
  d.setZero(V.rows(), 3); // di_x, di_y, di_z
  Eigen::MatrixXd D, N;
  D.setZero(V.rows(), 3);
  D = V - B;
  igl::per_vertex_normals(B, F, N);
  std::vector<std::vector<int> > VV;
  igl::adjacency_list(F, VV);
  Eigen::RowVector3d n, x, vi, vj, y; 
  for (int i = 0; i < V.rows(); ++i)
  {
    n = N.row(i);
    vi = B.row(i);
    double max_project = 0;
    for (int j = 0; j < VV[i].size(); ++j)
    {
      vj = B.row(VV[i][j]) - vi;
      vj -= n.dot(vj) * n;
      if (vj.squaredNorm() > max_project)
      {
        edge_index(i) = VV[i][j];
        max_project = vj.squaredNorm();
        x = vj;
      }
    }
    x.normalize();
    y = n.cross(x);
    d(i, 0) = D.row(i).dot(x);
    d(i, 1) = D.row(i).dot(y);
    d(i, 2) = D.row(i).dot(n);
  }
}

bool solve(Viewer& viewer)
{
  /**** Add your code for computing the deformation from handle_vertex_positions and handle_vertices here (replace following line) ****/
  // igl::slice_into(handle_vertex_positions, handle_vertices, 1, V);
  // cout << "solve called" << endl;
  // 1.3
  R = -A_fc * handle_vertex_positions;
	V_f = solver.solve(R);
  B_p.setZero(V.rows(), 3);
	igl::slice_into(V_f, free_vertices, 1, B_p);
	igl::slice_into(handle_vertex_positions, handle_vertices, 1, B_p);

  // part 2
  if (use_deformation_transfer)
  {
    // compute Si
    Eigen::MatrixXd SS, NN, N_p, Si;
    S_p.setZero(V.rows(), 3);
    SS.setZero(F.rows() * 3, 3);
    Eigen::Vector3d q1, q2, q3, q1_p, q2_p, q3_p, n, n_p;
    igl::per_face_normals(B, F, NN);
    igl::per_face_normals(B_p, F, N_p);
    Si.setZero(3, 3);
    for (int i = 0; i < F.rows(); ++i)
    {
      q1 = B.row(F(i, 0));
      q2 = B.row(F(i, 1));
      q3 = B.row(F(i, 2));
      q1_p = B_p.row(F(i, 0));
      q2_p = B_p.row(F(i, 1));
      q3_p = B_p.row(F(i, 2));
      n = NN.row(i);
      n_p = N_p.row(i);

      Si << q1_p - q3_p, q2_p - q3_p, n_p;
      Eigen::Matrix3d t;
      t << q1 - q3, q2 - q3, n;
      Si = Si * t.inverse();
      SS.block(i * 3, 0, 3, 3) = Si.transpose();
    }
    // compute G
    Eigen::SparseMatrix<double> G;
    Eigen::SparseMatrix<double, RowMajor> GG;
    // Eigen::MatrixXd GG;
    igl::grad(V, F, G);
    // shuffle G
    GG.resize(G.rows(), G.cols());
    for (int i = 0; i < F.rows(); ++i)
    {
      GG.row(3 * i) = G.row(i);
      GG.row(3 * i + 1) = G.row(i + F.rows());
      GG.row(3 * i + 2) = G.row(i + F.rows() * 2);
    }

    // compute D
    VectorXd area_vec;
		igl::doublearea(V, F, area_vec);
    // cout << "area size: " << endl;
		// cout << area_vec.size() << endl;
    VectorXd area_vec_3(area_vec.size() * 3);
    for (int i = 0; i < area_vec.size(); ++i)
    {
      area_vec(i) /= 2;
      area_vec_3(3 * i) = area_vec(i);
      area_vec_3(3 * i + 1) = area_vec(i);
      area_vec_3(3 * i + 2) = area_vec(i);
    }
    // cout << area_vec_3.size() << endl;
    SparseMatrix<double> D_area;
    igl::diag(area_vec_3, D_area);
    // cout << "D_area size:  " << endl;
    // cout << D_area.rows() << "  " << D_area.cols() << endl;

    // compute S'
    Eigen::MatrixXd Rhs, R_f, V_f2;
    Eigen::SparseMatrix<double> Lhs, G_f, G_c;
    Rhs = GG.transpose() * D_area * SS;
    Lhs = GG.transpose() * D_area * GG;
    igl::slice(Lhs, free_vertices, handle_vertices, G_c);
    igl::slice(Lhs, free_vertices, free_vertices, G_f);
    // cout << "G_c size:  " << endl;
    // cout << G_c.rows() << "  " << G_c.cols() << endl;
    // cout << "G_f size:  " << endl;
    // cout << G_f.rows() << "  " << G_f.cols() << endl;
    igl::slice(Rhs, free_vertices, 1, R_f);
    R_f = R_f - G_c * handle_vertex_positions;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> deform_solver(G_f);
    V_f2 = deform_solver.solve(R_f);
    // cout << "V_F2 size: " << endl;
    // cout << V_f2.rows() << "  " << V_f2.cols() << endl;
    igl::slice_into(V_f2, free_vertices, 1, S_p);
    igl::slice_into(handle_vertex_positions, handle_vertices, 1, S_p);
  }     

  // cout << "exit part2 " << endl;
  // 1.4
  if (use_deformation_transfer == 0) {
    Eigen::MatrixXd D, N;
    S_p.setZero(V.rows(), 3);
    igl::per_vertex_normals(B_p, F, N);
    Eigen::RowVector3d n, x, vi, vj, y;
    for (int i = 0; i < V.rows(); ++i)
    {
      n = N.row(i);
      vi = B_p.row(i);
      vj = B_p.row(edge_index[i]);
      x = vj - vi;
      x -= n.dot(x) * n;
      x.normalize();
      y = n.cross(x);
      S_p.row(i) = B_p.row(i) + d(i, 0) * x + d(i, 1) * y + d(i, 2) * n;
    }
  }

  switch (mesh_shape)
  {
  case 1:
    // cout << "B" << endl;
    V = B;
    break;
  case 2:
    // cout << "B_p" << endl;
    V = B_p;
    break;
  case 3:
    // cout << "S_p" << endl;
    V = S_p;
    break;
  default:
    break;
  }

  Eigen::MatrixXd normals;
	igl::per_vertex_normals(V, F, normals);
	viewer.data().set_normals(normals);

  return true;
};

void get_new_handle_locations()
{
  int count = 0;
  for (long vi = 0; vi < V.rows(); ++vi)
    if (handle_id[vi] >= 0)
    {
      Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();
      if (handle_id[vi] == moving_handle) {
        if (mouse_mode == TRANSLATE)
          goalPosition += translation;
        else if (mouse_mode == ROTATE) {
            Eigen::RowVector3f  goalPositionCopy = goalPosition;
            goalPosition -= handle_centroids.row(moving_handle).cast<float>();
            igl::rotate_by_quat(goalPosition.data(), rotation.data(), goalPositionCopy.data());
            goalPosition = goalPositionCopy;
            goalPosition += handle_centroids.row(moving_handle).cast<float>();
        }
      }
      handle_vertex_positions.row(count++) = goalPosition.cast<double>();
    }
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  viewer.data().clear();
  viewer.data().set_mesh(V, F);

  viewer.core().align_camera_center(V);
  V_cp = V;
  handle_id.setConstant(V.rows(), 1, -1);
  // Initialize selector
  lasso = std::unique_ptr<Lasso>(new Lasso(V, F, viewer));

  selected_v.resize(0,1);

  return true;
}

int main(int argc, char *argv[])
{
  if(argc != 2) {
    cout << "Usage assignment5 mesh.off>" << endl;
    load_mesh("../data/woody-lo.off");
  }
  else
  {
    load_mesh(argv[1]);
  }

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();
    ImGui::InputInt("mesh_shape", &mesh_shape, 0, 0);
    // Add new group
    if (ImGui::CollapsingHeader("Deformation Controls", ImGuiTreeNodeFlags_DefaultOpen))
    {
          int mouse_mode_type = static_cast<int>(mouse_mode);

          if (ImGui::Combo("Mouse Mode", &mouse_mode_type, "SELECT\0TRANSLATE\0ROTATE\0NONE\0"))
          {
            mouse_mode = static_cast<MouseMode>(mouse_mode_type);
          }

          if (ImGui::Button("Clear Selection", ImVec2(-1,0)))
          {
            selected_v.resize(0,1);
          }

          if (ImGui::Button("Apply Selection", ImVec2(-1,0)))
          {
            applySelection();
          }

          if (ImGui::Button("Clear Constraints", ImVec2(-1,0)))
          {
            handle_id.setConstant(V.rows(),1,-1);
          }
          if(ImGui::Checkbox("Deformation Transfer", &use_deformation_transfer)){}
    }
  };

  viewer.callback_key_down = callback_key_down;
  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;
  viewer.callback_pre_draw = callback_pre_draw;

  viewer.data().point_size = 10;
  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
  viewer.launch();
}


bool callback_mouse_down(Viewer& viewer, int button, int modifier)
{
  if (button == (int) Viewer::MouseButton::Right)
    return false;

  down_mouse_x = viewer.current_mouse_x;
  down_mouse_y = viewer.current_mouse_y;

  if (mouse_mode == SELECT)
  {
    if (lasso->strokeAdd(viewer.current_mouse_x, viewer.current_mouse_y) >=0)
      doit = true;
    else
      lasso->strokeReset();
  }
  else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
    int vi = lasso->pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
    if(vi>=0 && handle_id[vi]>=0)  //if a region was found, mark it for translation/rotation
    {
      moving_handle = handle_id[vi];
      get_new_handle_locations();
      doit = true;
    }
  }
  return doit;
}

bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y)
{
  // cout << "mouse_move" << endl;
  if (!doit)
    return false;
  if (mouse_mode == SELECT)
  {
    lasso->strokeAdd(mouse_x, mouse_y);
    return true;
  }
  if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
    if (mouse_mode == TRANSLATE) {
      translation = computeTranslation(viewer,
                                       mouse_x,
                                       down_mouse_x,
                                       mouse_y,
                                       down_mouse_y,
                                       handle_centroids.row(moving_handle));
    }
    else {
      rotation = computeRotation(viewer,
                                 mouse_x,
                                 down_mouse_x,
                                 mouse_y,
                                 down_mouse_y,
                                 handle_centroids.row(moving_handle));
    }
    get_new_handle_locations();
#ifndef UPDATE_ONLY_ON_UP
    // cout << "mouse move and up" << endl;
    solve(viewer);
    down_mouse_x = mouse_x;
    down_mouse_y = mouse_y;
#endif
    return true;

  }
  return false;
}

bool callback_mouse_up(Viewer& viewer, int button, int modifier)
{
  // cout << "mouse_up" << endl;
  if (!doit)
    return false;
  doit = false;
  if (mouse_mode == SELECT)
  {
    selected_v.resize(0,1);
    lasso->strokeFinish(selected_v);
    return true;
  }

  if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
#ifdef UPDATE_ONLY_ON_UP
    if(moving_handle>=0)
      solve(viewer);
#endif
    translation.setZero();
    rotation.setZero(); rotation[3] = 1.;
    moving_handle = -1;

    compute_handle_centroids();

    return true;
  }

  return false;
};


bool callback_pre_draw(Viewer& viewer)
{
  // initialize vertex colors
  vertex_colors = Eigen::MatrixXd::Constant(V.rows(),3,.9);

  // first, color constraints
  int num = handle_id.maxCoeff();
  if (num == 0)
    num = 1;
  for (int i = 0; i<V.rows(); ++i)
    if (handle_id[i]!=-1)
    {
      int r = handle_id[i] % MAXNUMREGIONS;
      vertex_colors.row(i) << regionColors[r][0], regionColors[r][1], regionColors[r][2];
    }
  // then, color selection
  for (int i = 0; i<selected_v.size(); ++i)
    vertex_colors.row(selected_v[i]) << 131./255, 131./255, 131./255.;

  viewer.data().set_colors(vertex_colors);
  viewer.data().V_material_specular.fill(0);
  viewer.data().V_material_specular.col(3).fill(1);
  viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;


  //clear points and lines
  viewer.data().set_points(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
  viewer.data().set_edges(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXi::Zero(0,3), Eigen::MatrixXd::Zero(0,3));

  //draw the stroke of the selection
  for (unsigned int i = 0; i<lasso->strokePoints.size(); ++i)
  {
    viewer.data().add_points(lasso->strokePoints[i],Eigen::RowVector3d(0.4,0.4,0.4));
    if(i>1)
      viewer.data().add_edges(lasso->strokePoints[i-1], lasso->strokePoints[i], Eigen::RowVector3d(0.7,0.7,0.7));
  }

  // update the vertex position all the time
  viewer.data().V.resize(V.rows(),3);
  viewer.data().V << V;

  viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_POSITION;

#ifdef UPDATE_ONLY_ON_UP
  //draw only the moving parts with a white line
  if (moving_handle>=0)
  {
    Eigen::MatrixXd edges(3*F.rows(),6);
    int num_edges = 0;
    for (int fi = 0; fi<F.rows(); ++fi)
    {
      int firstPickedVertex = -1;
      for(int vi = 0; vi<3 ; ++vi)
        if (handle_id[F(fi,vi)] == moving_handle)
        {
          firstPickedVertex = vi;
          break;
        }
      if(firstPickedVertex==-1)
        continue;


      Eigen::Matrix3d points;
      for(int vi = 0; vi<3; ++vi)
      {
        int vertex_id = F(fi,vi);
        if (handle_id[vertex_id] == moving_handle)
        {
          int index = -1;
          // if face is already constrained, find index in the constraints
          (handle_vertices.array()-vertex_id).cwiseAbs().minCoeff(&index);
          points.row(vi) = handle_vertex_positions.row(index);
        }
        else
          points.row(vi) =  V.row(vertex_id);

      }
      edges.row(num_edges++) << points.row(0), points.row(1);
      edges.row(num_edges++) << points.row(1), points.row(2);
      edges.row(num_edges++) << points.row(2), points.row(0);
    }
    edges.conservativeResize(num_edges, Eigen::NoChange);
    viewer.data().add_edges(edges.leftCols(3), edges.rightCols(3), Eigen::RowVector3d(0.9,0.9,0.9));

  }
#endif
  return false;

}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers)
{
  bool handled = false;
  if (key == 'S')
  {
    mouse_mode = SELECT;
    handled = true;
  }

  if ((key == 'T') && (modifiers == IGL_MOD_ALT))
  {
    mouse_mode = TRANSLATE;
    handled = true;
  }

  if ((key == 'R') && (modifiers == IGL_MOD_ALT))
  {
    mouse_mode = ROTATE;
    handled = true;
  }
  if (key == 'A')
  {
    applySelection();
    callback_key_down(viewer, '1', 0);
    handled = true;
  }
  // if (key == '2')
  // viewer.core().align_camera_position(B);
  //  }

  //viewer.ngui->refresh();
  return handled;
}


void onNewHandleID()
{
  //store handle vertices too
  int numFree = (handle_id.array() == -1).cast<int>().sum();
  int num_handle_vertices = V.rows() - numFree;
  handle_vertices.setZero(num_handle_vertices);
  handle_vertex_positions.setZero(num_handle_vertices,3);
  free_vertices.setZero(numFree);
  V_f.setZero(numFree, 3);

  int count = 0;
  int count_free = 0;
  for (long vi = 0; vi<V.rows(); ++vi)
    if(handle_id[vi] >=0)
      handle_vertices[count++] = vi;
    else
      free_vertices[count_free++] = vi;

  // cout << "count_free.  " << count_free << "  num_free.  " << numFree << endl;
  // cout << "count_handle.  " << count << "  num_handle.  " << num_handle_vertices << endl;
  compute_handle_centroids();
}

void applySelection()
{
  int index = handle_id.maxCoeff()+1;
  for (int i =0; i<selected_v.rows(); ++i)
  {
    const int selected_vertex = selected_v[i];
    if (handle_id[selected_vertex] == -1)
      handle_id[selected_vertex] = index;
  }
  selected_v.resize(0,1);

  onNewHandleID();
  get_smooth(); // 1.2
  encode_displacement(); // 1.4
}

void compute_handle_centroids()
{
  //compute centroids of handles
  int num_handles = handle_id.maxCoeff()+1;
  handle_centroids.setZero(num_handles,3);

  Eigen::VectorXi num; num.setZero(num_handles,1);
  for (long vi = 0; vi<V.rows(); ++vi)
  {
    int r = handle_id[vi];
    if ( r!= -1)
    {
      handle_centroids.row(r) += V.row(vi);
      num[r]++;
    }
  }

  for (long i = 0; i<num_handles; ++i)
    handle_centroids.row(i) = handle_centroids.row(i).array()/num[i];

}

//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector3f computeTranslation (Viewer& viewer,
                                    int mouse_x,
                                    int from_x,
                                    int mouse_y,
                                    int from_y,
                                    Eigen::RowVector3d pt3D)
{
  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.data().model;
  //project the given point (typically the handle centroid) to get a screen space depth
  Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                      modelview,
                                      viewer.core().proj,
                                      viewer.core().viewport);
  float depth = proj[2];

  double x, y;
  Eigen::Vector3f pos1, pos0;

  //unproject from- and to- points
  x = mouse_x;
  y = viewer.core().viewport(3) - mouse_y;
  pos1 = igl::unproject(Eigen::Vector3f(x,y,depth),
                        modelview,
                        viewer.core().proj,
                        viewer.core().viewport);


  x = from_x;
  y = viewer.core().viewport(3) - from_y;
  pos0 = igl::unproject(Eigen::Vector3f(x,y,depth),
                        modelview,
                        viewer.core().proj,
                        viewer.core().viewport);

  //translation is the vector connecting the two
  Eigen::Vector3f translation = pos1 - pos0;
  return translation;

}


//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector4f computeRotation(Viewer& viewer,
                                int mouse_x,
                                int from_x,
                                int mouse_y,
                                int from_y,
                                Eigen::RowVector3d pt3D)
{

  Eigen::Vector4f rotation;
  rotation.setZero();
  rotation[3] = 1.;

  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.data().model;

  //initialize a trackball around the handle that is being rotated
  //the trackball has (approximately) width w and height h
  double w = viewer.core().viewport[2]/8;
  double h = viewer.core().viewport[3]/8;

  //the mouse motion has to be expressed with respect to its center of mass
  //(i.e. it should approximately fall inside the region of the trackball)

  //project the given point on the handle(centroid)
  Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                      modelview,
                                      viewer.core().proj,
                                      viewer.core().viewport);
  proj[1] = viewer.core().viewport[3] - proj[1];

  //express the mouse points w.r.t the centroid
  from_x -= proj[0]; mouse_x -= proj[0];
  from_y -= proj[1]; mouse_y -= proj[1];

  //shift so that the range is from 0-w and 0-h respectively (similarly to a standard viewport)
  from_x += w/2; mouse_x += w/2;
  from_y += h/2; mouse_y += h/2;

  //get rotation from trackball
  Eigen::Vector4f drot = viewer.core().trackball_angle.coeffs();
  Eigen::Vector4f drot_conj;
  igl::quat_conjugate(drot.data(), drot_conj.data());
  igl::trackball(w, h, float(1.), rotation.data(), from_x, from_y, mouse_x, mouse_y, rotation.data());

  //account for the modelview rotation: prerotate by modelview (place model back to the original
  //unrotated frame), postrotate by inverse modelview
  Eigen::Vector4f out;
  igl::quat_mult(rotation.data(), drot.data(), out.data());
  igl::quat_mult(drot_conj.data(), out.data(), rotation.data());
  return rotation;
}
