#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/dijkstra.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

bool showingUV = false;
bool freeBoundary = true;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;
int distortion_type = 1; // 1 angle, 2 area, 3 length
bool show_distortion = false;
MatrixXd color;
bool arap_initialized = 0;

void Redraw()
{
	viewer.data().clear();

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);

    if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}
	if (show_distortion) {
		viewer.data().show_texture = false;
		viewer.data().set_colors(color);
	}
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
	// ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);
	d.setZero(indices.rows() * 2);
	C.resize(indices.rows() * 2, V.rows() * 2);
	for (int i = 0; i < indices.rows(); ++i)
	{
		d(i) = positions(i, 0);
		d(i + indices.rows()) = positions(i, 1);
	}
	vector<Eigen::Triplet<double> > index_list;
	for (int i = 0; i < indices.rows(); ++i)
	{
		index_list.push_back(Eigen::Triplet<double>(i, indices[i], 1));
		index_list.push_back(Eigen::Triplet<double>(i + indices.rows(), indices[i] + V.rows(), 1));
	}
	C.setFromTriplets(index_list.begin(), index_list.end());
}

void computeParameterization(int type)
{
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;

	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary)
	{
		// The boundary vertices should be fixed to positions on the unit disc. Find these position and
		// save them in the #V x 2 matrix fixed_UV_position.
		igl::boundary_loop(F, fixed_UV_indices);
		igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
		// cout << fixed_UV_indices << endl;
		// cout << "uv indices size: " << endl;
		// cout << V.rows() << "   " << fixed_UV_indices.rows() << endl;
		// cout << fixed_UV_positions.rows() << endl;
	}
	else
	{
		// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
		// Fix two points on the boundary

		set<int> targets;
		vector<vector<int> > VV;
		VectorXd min_distance;
		VectorXi previous;
		double max_distance = 0;
		int v1, v2;
		igl::adjacency_list(F, VV);
		for (int i = 0; i < V.rows(); ++i)
		{
			igl::dijkstra(i, targets, VV, min_distance, previous);
			int max_index = 0;
			double tmp_max = min_distance.maxCoeff(&max_index);
			// cout << tmp_max << endl;
			if (tmp_max > max_distance)
			{
				v1 = i;
				v2 = max_index;
				max_distance = tmp_max;
			}
		}
		fixed_UV_indices.resize(2);
		fixed_UV_positions.resize(2, 2);
		fixed_UV_indices << v1, v2;
		fixed_UV_positions << 0, 0, 1, 1;
	}

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh
		A.resize(V.rows() * 2, V.rows() * 2);
		b.setZero(V.rows() * 2);

		SparseMatrix<double> Adj;
  	    igl::adjacency_matrix(F,Adj);
		// sum each row 
    	SparseVector<double> Asum;
    	igl::sum(Adj, 1, Asum);
    	// Convert row sums into diagonal of sparse matrix
		SparseMatrix<double> Adiag;
		igl::diag(Asum,Adiag);
    	// Build uniform laplacian
		SparseMatrix<double> U;
		U = Adj-Adiag;

		vector<Eigen::Triplet<double> > index_list;	
		for (int i = 0; i < U.outerSize(); ++i) 
		{
			for (SparseMatrix<double>::InnerIterator it(U, i); it; ++it) 
			{
				index_list.push_back(Eigen::Triplet<double> (it.row(), it.col(), it.value()));
				index_list.push_back(Eigen::Triplet<double> (it.row() + V.rows(), it.col() + V.rows(), it.value()));
			}
		}
		A.setFromTriplets(index_list.begin(), index_list.end());
	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~
		A.resize(V.rows() * 2, V.rows() * 2);
		b.setZero(V.rows() * 2);
		SparseMatrix<double> L;
		igl::cotmatrix(V, F, L);

		vector<Eigen::Triplet<double> > index_list;		
		for (int i = 0; i < L.outerSize(); ++i) 
		{
			for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it) 
			{
				index_list.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
				index_list.push_back(Eigen::Triplet<double>(it.row() + V.rows(), it.col() + V.rows(), it.value()));
			}
		}
		A.setFromTriplets(index_list.begin(), index_list.end());
	}

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!
		A.resize(V.rows() * 2, V.rows() * 2);
		b.setZero(2 * V.rows());

		VectorXd area_vec;
		igl::doublearea(V, F, area_vec);
		cout << area_vec.size() << endl;
		SparseMatrix<double> Area_matrix;
		igl::diag(area_vec, Area_matrix);
		// cout << "area matrix size:   " << endl;
		// cout << Area_matrix.rows() << "   " << Area_matrix.cols() << endl;

		SparseMatrix<double> Dx, Dy, DxT, DyT;
		computeSurfaceGradientMatrix(Dx, Dy);
		DxT = Dx.transpose();
		DyT = Dy.transpose();

		SparseMatrix<double> A1, A2, A11, A12, A21, A22;
		
		A11 = (DxT * Area_matrix * Dx + DyT * Area_matrix * Dy);
		A12 = (-DxT * Area_matrix * Dy + DyT * Area_matrix * Dx);
		A21 = (-DyT * Area_matrix * Dx + DxT * Area_matrix * Dy);
		A22 = (DxT * Area_matrix * Dx + DyT * Area_matrix * Dy);
		igl::cat(1, A11, A21, A1);
		// cout << "A1 size:  " << endl;
		// cout << A1.rows() << "   " << A1.cols() << endl;
		igl::cat(1, A12, A22, A2);
		igl::cat(2, A1, A2, A);
	}

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices
		A.resize(V.rows() * 2, V.rows() * 2);
		b.resize(V.rows() * 2);

		SparseMatrix<double> Dx, Dy, DxT, DyT;
		computeSurfaceGradientMatrix(Dx, Dy);
		DxT = Dx.transpose();
		DyT = Dy.transpose();

		VectorXd area_vec;
		igl::doublearea(V, F, area_vec);
		SparseMatrix<double> Area_matrix;
		igl::diag(area_vec, Area_matrix);
		// set A
		SparseMatrix<double> A1, A2, A11, A12, A21, A22;
		A11 = (DxT * Area_matrix * Dx * 0.5 + DyT * Area_matrix * Dy * 0.5);
		A12.resize(V.rows(),V.rows());
		A21.resize(V.rows(),V.rows());
		A22 = A11;
		igl::cat(1, A11, A21, A1);
		igl::cat(1, A12, A22, A2);
		igl::cat(2, A1, A2, A);

		// set b
		VectorXd R11, R12, R21, R22;
		R11.setZero(F.rows());
		R12.setZero(F.rows());
		R21.setZero(F.rows());
		R22.setZero(F.rows());
		MatrixXd DxUV = Dx * UV;
		MatrixXd DyUV = Dy * UV;
		for (int i = 0; i < F.rows(); ++i)
		{
			Matrix2d J, U, Sigma, V, R, VT;
			J << DxUV(i, 0), DyUV(i, 0),
				 DxUV(i, 1), DyUV(i, 1);
			SSVD2x2(J, U, Sigma, V); 
			VT = V.transpose();
			R = U * VT;
			// cout << "R determinant:   " << R.determinant() << endl;
			R11(i) = R(0, 0);
			R12(i) = R(0, 1);
			R21(i) = R(1, 0);
			R22(i) = R(1, 1);
		}

		VectorXd b1, b2;
		b1 = DxT * Area_matrix * R11 * 0.5 + DyT * Area_matrix * R12 * 0.5;
		b2 = DxT * Area_matrix * R21 * 0.5 + DyT * Area_matrix * R22 * 0.5;
		igl::cat(1, b1, b2, b);
	}

	// Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail
	SparseMatrix<double> Ll, Lr, L, zero;
	igl::cat(1, A, C, Ll);
	zero.resize(C.rows(), C.rows());
	SparseMatrix<double> CT = C.transpose();
	igl::cat(1, CT, zero, Lr);
	igl::cat(2, Ll, Lr, L);
	VectorXd r, x;
	r.resize(b.rows() + d.rows());
	x.resize(b.rows() + d.rows());
	igl::cat(1, b, d, r);

	SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
	solver.analyzePattern(L);
	solver.factorize(L);
	x = solver.solve(r);

	// The solver will output a vector
	UV.resize(V.rows(), 2);
	UV.col(0) = x.block(0, 0, V.rows(), 1);
	UV.col(1) = x.block(V.rows(), 0, V.rows(), 1);
}

void compute_distortion()
{
	VectorXd Distortion;
	Distortion.setZero(F.rows());
	SparseMatrix<double> Dx, Dy;
	computeSurfaceGradientMatrix(Dx, Dy);
	MatrixXd DxUV = Dx * UV;
	MatrixXd DyUV = Dy * UV;
	for (int i = 0; i < F.rows(); ++i)
	{
		Matrix2d J;
		J << DxUV(i, 0), DyUV(i, 0),
			 DxUV(i, 1), DyUV(i, 1);
		if (distortion_type == 1)
		{
			Matrix2d JT, I;
			JT = J.transpose();
			I << 1, 0,
				 0, 1;
			Distortion(i) = (J + JT - J.trace() * I).squaredNorm();
		}
		else if (distortion_type == 2)
		{
			Distortion(i) = (J.determinant() - 1) * (J.determinant() - 1);
		}
		else if (distortion_type == 3)
		{
			Matrix2d U, Sigma, V, R, VT;
			SSVD2x2(J, U, Sigma, V); 
			VT = V.transpose();
			R = U * VT;
			// cout << "R determinant:   " << R.determinant() << endl;
			Distortion(i) = (J - R).squaredNorm();
		}
	}
	
	color.resize(F.rows(), 3);
	double min_distortion = Distortion.minCoeff();
	double max_distortion = Distortion.maxCoeff();
	for (int i = 0; i < Distortion.rows(); ++i)
	{
		color(i, 0) = 1;
		color(i, 1) = 1 - (Distortion(i) - min_distortion) / (max_distortion - min_distortion);
		color(i, 2) = color(i, 1);
	}
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
		computeParameterization(key);
		show_distortion = false;
		break;
	case '4':
		if (arap_initialized == 0)
		{
			computeParameterization('3');
			arap_initialized = 1;
		}
		computeParameterization(key);
		show_distortion = false;
		break;
	case '5':
			// Add your code for detecting and displaying flipped triangles in the
			// UV domain here
		compute_distortion();
		show_distortion = true;
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core();
      viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core();
        viewer.core() = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  Redraw();
  viewer.core().align_camera_center(V);
  showingUV = false;

  return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core();
	temp2D = viewer.core();
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
  if(argc != 2) {
    cout << "Usage ex4_bin <mesh.off/obj>" << endl;
    load_mesh("../data/cathead.obj");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmaterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);

			// TODO: Add more parameters to tweak here...
			ImGui::InputInt("Distortion type", &distortion_type, 0, 1);
		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_init = callback_init;

  viewer.launch();
}
