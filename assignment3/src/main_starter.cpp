#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/jet.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/sum.h>
#include <igl/speye.h>
#include <igl/bfs.h>
#include <igl/cotmatrix.h>
#include <igl/principal_curvature.h>
#include <imgui/imgui.h>
/*** insert any libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/doublearea.h>
#include <igl/adjacency_list.h>
#include <set>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #Vx3
Eigen::MatrixXd V;
// Face array, #Fx3
Eigen::MatrixXi F;
//Face normals #Fx3
Eigen::MatrixXd FN;
//Vertex normals #Vx3
Eigen::MatrixXd VN;

// Per-vertex uniform normal array, #Vx3
Eigen::MatrixXd N_uniform;
// Per-vertex area-weighted normal array, #Vx3
Eigen::MatrixXd N_area;
// Per-vertex mean-curvature normal array, #Vx3
Eigen::MatrixXd N_meanCurvature;
// Per-vertex PCA normal array, #Vx3
Eigen::MatrixXd N_PCA;
// Per-vertex quadratic fitted normal array, #Vx3
Eigen::MatrixXd N_quadraticFit;

// Per-vertex mean curvature, #Vx3
Eigen::VectorXd K_mean;
// Per-vertex Gaussian curvature, #Vx3
Eigen::VectorXd K_Gaussian;
// Per-vertex minimal principal curvature, #Vx3
Eigen::VectorXd K_min_principal;
// Per-vertex maximal principal curvature, #Vx3
Eigen::VectorXd K_max_principal;
// Per-vertex color array, #Vx3
Eigen::MatrixXd colors_per_vertex;

// Explicitely smoothed vertex array, #Vx3
Eigen::MatrixXd V_expLap;
// Implicitely smoothed vertex array, #Vx3
Eigen::MatrixXd V_impLap;
// Bilateral smoothed vertex array, #Vx3
Eigen::MatrixXd V_bilateral;

vector<vector<int> > VV;

void uniform_vertex_normal()
{
    N_uniform.setZero(V.rows(), 3);
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j = 0; j < F.cols(); ++j)
        {
            N_uniform.row(F(i, j)) += FN.row(i);
        }
    }
    N_uniform.rowwise().normalize();    
    igl::per_vertex_normals(V, F, VN);
    for (int i = 0; i < N_uniform.rows(); ++i)
    {
        if (N_uniform.row(i).dot(VN.row(i)) < 0)
        {
            N_uniform.row(i) = -N_uniform.row(i);
        }
    }
}

void area_weight_normal()
{
    N_area.setZero(V.rows(), 3);
    Eigen::VectorXd area;
    igl::doublearea(V, F, area);
    area = area.array() / 2;
    
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j = 0; j < F.cols(); ++j)
        {
            N_area.row(F(i, j)) += FN.row(i) * area(i);
        }
    }
    N_area.rowwise().normalize();
    igl::per_vertex_normals(V, F, VN);
    for (int i = 0; i < N_area.rows(); ++i)
    {
        // cout << N_area.row(i).dot(VN.row(i)) << endl;
        if (N_area.row(i).dot(VN.row(i)) < 0)
        {
            N_area.row(i) = -N_area.row(i);
        }
    }
}

void mean_curvature_normal()
{
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);

    N_meanCurvature.setZero(V.rows(), 3);
    // for (int i = 0; i < L.outerSize(); ++i)
    // {
    //     for (Eigen::SparseMatrix<double>::InnerIterator it(L,i); it; ++it)
    //     {
    //         if (it.row() != it.col())
    //             N_meanCurvature.row(it.row()) += it.value() * (V.row(it.row()) - V.row(it.col()));
    //     }
    // }
    N_meanCurvature = L * V;
    N_meanCurvature.rowwise().normalize();

    igl::per_vertex_normals(V, F, VN);
    for (int i = 0; i < N_meanCurvature.rows(); ++i)
    {
        // cout << N_meanCurvature.row(i).dot(VN.row(i)) << endl;
        if (N_meanCurvature.row(i).dot(VN.row(i)) < 0)
        {
            N_meanCurvature.row(i) = -N_meanCurvature.row(i);
        }
    }
}

void get_neighbor(int vertex_id, set<int> &neighbor, int k = 1)
{
    neighbor.insert(vertex_id);
    for (int i = 0; i < VV[vertex_id].size(); ++i)
    {
        neighbor.insert(VV[vertex_id][i]);
    }
    if (k == 2)
    {
        vector<int> tmp;
        for (auto it = neighbor.begin(); it != neighbor.end(); ++it)
        {
            tmp.push_back(*it);
        }
        for (int i = 0; i < tmp.size(); ++i)
        {
            for (int j = 0; j < VV[tmp[i]].size(); ++j)
            {
                neighbor.insert(VV[tmp[i]][j]);
            }
        }
    }
    // cout << "vertex " << vertex_id << "  is connect to:  ";
    // for (auto it = neighbor.begin(); it != neighbor.end(); ++it)
    // {
    //     cout << *it << "  ";
    // }
    // cout << endl;
}

void pca_normal()
{
    igl::adjacency_list(F, VV);
    N_PCA.setZero(V.rows(), 3);
    igl::per_vertex_normals(V, F, VN);
    int k = 1;
    for (int i = 0; i < V.rows(); ++i)
    {
        set<int> neighbor;
        get_neighbor(i, neighbor, k);
        Eigen::MatrixXd C;
        C.setZero(neighbor.size(), 3);
        int idx = 0;
        // cout << "neighbor size. " << neighbor.size() << endl;
        for (auto it = neighbor.begin(); it != neighbor.end(); ++it)
        {
            C.row(idx) << V.row(*it);
            ++idx;
        }
        Eigen::RowVectorXd mean = C.colwise().mean();
        Eigen::MatrixXd aligned = C.rowwise() - mean;
        Eigen::MatrixXd cov = aligned.adjoint() * aligned;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(cov);
        int pos = -1;
        eigensolver.eigenvalues().minCoeff(&pos);
        Eigen::RowVectorXd normal = eigensolver.eigenvectors().col(pos);
        if (normal.dot(VN.row(i)) < 0)
            normal = -normal;
        N_PCA.row(i) << normal;
        // cout << normal.dot(VN.row(i)) << endl;
    }
}

void quadratic_fit_normal()
{
    igl::adjacency_list(F, VV);
    N_quadraticFit.setZero(V.rows(), 3);
    igl::per_vertex_normals(V, F, VN);
    int k = 1;
    for (int i = 0; i < V.rows(); ++i)
    {
        set<int> neighbor;
        get_neighbor(i, neighbor, k);
        Eigen::MatrixXd C;
        C.setZero(neighbor.size(), 3);
        int idx = 0;
        for (auto it = neighbor.begin(); it != neighbor.end(); ++it)
        {
            C.row(idx) << V.row(*it);
            ++idx;
        }
        Eigen::RowVectorXd mean = C.colwise().mean();
        Eigen::MatrixXd aligned = C.rowwise() - mean;
        Eigen::MatrixXd cov = aligned.adjoint() * aligned;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(cov);
        Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
        if (eigenvectors.determinant() < 0)
            eigenvectors = eigenvectors * -1.0;
        int pos = -1;
        eigensolver.eigenvalues().minCoeff(&pos);
        Eigen::RowVectorXd normal = eigenvectors.col(pos);
        // cout << "pos. " << pos << endl; // 0

        // quadratic fitting
        Eigen::MatrixXd basis;
        Eigen::MatrixXd values;
        values.setZero(neighbor.size(), 1);
        basis.resize(neighbor.size(), 5);
        Eigen::RowVectorXd ev1;
        Eigen::RowVectorXd ev2;
        int flag = 0;
        for (int j = 0; j < 3; ++j)
        {
            if (j != pos)
            {
                if (flag == 0)
                {
                    ev1 = eigenvectors.col(j);
                    flag = 1;
                }
                else
                {
                    ev2 = eigenvectors.col(j);
                }
            }
        }
        for (int j = 0; j < C.rows(); ++j)
        {
            double u = (C.row(j) - V.row(i)).dot(ev1);
            double v = (C.row(j) - V.row(i)).dot(ev2);
            double f = (C.row(j) - V.row(i)).dot(normal);
            basis.row(j) << u * u, v * v, v * u, u, v;
            values.row(j) << f;
        }
        Eigen::VectorXd coef = basis.colPivHouseholderQr().solve(values);
        // cout << coef << endl;
        double du = coef(3);
        double dv = coef(4);
        Eigen::RowVector3d pu;
        pu << 1, 0, du;
        Eigen::RowVector3d pv;
        pv << 0, 1, dv;
        Eigen::RowVector3d quadratic_normal = pu.cross(pv);
        quadratic_normal.normalize();
        eigenvectors.col(0) = eigensolver.eigenvectors().col(2);  // change normal to the third column
        eigenvectors.col(2) = normal;
        quadratic_normal = quadratic_normal * eigenvectors.inverse(); // ***inverse

        if (quadratic_normal.dot(VN.row(i)) < 0)
            quadratic_normal = -quadratic_normal;
        // cout << quadratic_normal.dot(VN.row(i)) << endl;
        N_quadraticFit.row(i) << quadratic_normal;
    }
}
void implicit_smooth()
{
    Eigen::SparseMatrix<double> L;
    // Compute Laplace-Beltrami operator: #V by #V
    igl::cotmatrix(V,F,L);

    // Recompute just mass matrix on each step
    Eigen::SparseMatrix<double> M;
    Eigen::MatrixXd U;
    U = V_impLap;
    igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
    // Solve (M-delta*L) U = M*U
    double step_size = 0.01;
    const auto & S = (M - step_size * L);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
    assert(solver.info() == Eigen::Success);
    U = solver.solve(M*U).eval();
    V_impLap = U;
}

void explicit_smooth()
{
    Eigen::SparseMatrix<double> L,M,Minv;
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
    igl::invert_diag(M,Minv);
    Eigen::MatrixXd U;
    U = V_expLap;

    double step_size = 0.00001;
    V_expLap = U + step_size * Minv * L * U;
}

void explicit_smooth_uniform()
{
    Eigen::SparseMatrix<double> L(V.rows(), V.rows());
    // igl::cotmatrix(V,F,L);
    //declare list of non-zero elements (row, column, value)
    std::vector<Eigen::Triplet<double> > tripletList;
    tripletList.push_back(Eigen::Triplet<double>(0, V.rows() - 1, 0.5));
    tripletList.push_back(Eigen::Triplet<double>(0, 0, -1));
    tripletList.push_back(Eigen::Triplet<double>(0, 1, 0.5));
    for (int i = 1; i < V.rows() - 1; ++i)
    {
        tripletList.push_back(Eigen::Triplet<double>(i, i, -1));
        tripletList.push_back(Eigen::Triplet<double>(i, i + 1, 0.5));
        tripletList.push_back(Eigen::Triplet<double>(i, i - 1, 0.5));
    }
    tripletList.push_back(Eigen::Triplet<double>(V.rows() - 1, V.rows() - 1, -1));
    tripletList.push_back(Eigen::Triplet<double>(V.rows() - 1, 0, 0.5));
    tripletList.push_back(Eigen::Triplet<double>(V.rows() - 1, V.rows() - 2, 0.5));
    //construct matrix from the list
    L.setFromTriplets(tripletList.begin(), tripletList.end());

    Eigen::MatrixXd U;
    U = V_expLap;

    double step_size = 0.0001;
    V_expLap = U + step_size * L * U;
}

void bilateral_smoothing();
void area_weight_normal_b();

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing uniform vertex normals here:
        // store in N_uniform

        // Use igl::per_vertex_normals to orient your normals consistently.
        igl::per_face_normals(V, F, FN);
        uniform_vertex_normal();

        // Set the viewer normals.
        // N_uniform.rowwise().normalize();
        viewer.data().set_normals(N_uniform);
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing area-weighted vertex normals here:
        // store in N_area

        // Use igl::per_vertex_normals to orient your normals consistently.
        igl::per_face_normals(V, F, FN);
        area_weight_normal();

        // Set the viewer normals.
        // N_area.rowwise().normalize();
        viewer.data().set_normals(N_area);
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing mean-curvature vertex normals here:
        // store in N_meanCurvature

        // Use igl::per_vertex_normals to orient your normals consistently.
        mean_curvature_normal();

        // Set the viewer normals.
        // N_meanCurvature.rowwise().normalize();
        viewer.data().set_normals(N_meanCurvature);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing PCA vertex normals here:
        // store in N_PCA

        // Use igl::per_vertex_normals to orient your normals consistently.
        pca_normal();

        // Set the viewer normals.
        N_PCA.rowwise().normalize();
        viewer.data().set_normals(N_PCA);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing quadratic fitted vertex normals here:
        // store in N_quadraticFit

        // Use igl::per_vertex_normals to orient your normals consistently.
        quadratic_fit_normal();

        // Set the viewer normals.
        N_quadraticFit.rowwise().normalize();
        viewer.data().set_normals(N_quadraticFit);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete mean curvature:
        // store in K_mean
        K_mean.setZero(V.rows(), 1);

        // For visualization, better to normalize the range of K_mean with the maximal and minimal curvatures.
        // store colors in colors_per_vertex
        // Eigen::MatrixXd HN;
        // Eigen::SparseMatrix<double> L,M,Minv;
        // igl::cotmatrix(V,F,L);
        // igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
        // igl::invert_diag(M,Minv);
        // // Laplace-Beltrami of position
        // HN = -Minv*(L*V);
        // // Extract magnitude as mean curvature
        // Eigen::VectorXd H = HN.rowwise().norm();

        // Compute curvature directions via quadric fitting
        Eigen::MatrixXd PD1,PD2;
        Eigen::VectorXd PV1,PV2;
        igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
        // mean curvature
        K_mean = 0.5*(PV1+PV2);

        double mink = K_mean.minCoeff();
        double maxk = K_mean.maxCoeff();
        for (int i = 0; i < K_mean.rows(); ++i)
        {
            // cout << H(i) << "   " << K_mean(i) << endl;
            igl::jet((K_mean(i) - mink) / (maxk - mink), colors_per_vertex(i, 0), colors_per_vertex(i, 1), colors_per_vertex(i, 2));
        }
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '7') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete Gaussian curvature:
        // store in K_Gaussian
        K_Gaussian.setZero(V.rows(), 1);

        // For visualization, better to normalize the range of K_Gaussian with the maximal and minimal curvatures.
        // store colors in colors_per_vertex

        // Compute integral of Gaussian curvature
        igl::gaussian_curvature(V,F,K_Gaussian);
        // Compute mass matrix
        Eigen::SparseMatrix<double> M,Minv;
        igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
        igl::invert_diag(M,Minv);
        // Divide by area to get integral average
        K_Gaussian = (Minv*K_Gaussian).eval();
        double mink = K_Gaussian.minCoeff();
        double maxk = K_Gaussian.maxCoeff();
        for (int i = 0; i < K_Gaussian.rows(); ++i)
        {
            igl::jet((K_Gaussian(i) - mink) / (maxk - mink), colors_per_vertex(i, 0), colors_per_vertex(i, 1), colors_per_vertex(i, 2));
        }
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '8') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete minimal principal curvature:
        // store in K_min_principal
        K_min_principal.setZero(V.rows(), 1);
       
        // For visualization, better to normalize the range of K_min_principal with the maximal and minimal curvatures.
        // store colors in colors_per_vertex
        Eigen::MatrixXd PD_min,PD2;
        Eigen::VectorXd PV1,PV2;
        igl::principal_curvature(V,F,PD_min,PD2,PV1,PV2);
        K_min_principal = PV1;

        double mink = K_min_principal.minCoeff();
        double maxk = K_min_principal.maxCoeff();
        for (int i = 0; i < K_min_principal.rows(); ++i)
        {
            // cout << PV1(i) << "   " << PV2(i) << endl;
            igl::jet((K_min_principal(i) - mink) / (maxk - mink), colors_per_vertex(i, 0), colors_per_vertex(i, 1), colors_per_vertex(i, 2));
        }
  
        // Uncomment the code below to draw a blue segment parallel to the minimal curvature direction, 
        const double avg = igl::avg_edge_length(V,F);
        Eigen::RowVector3d blue(0.2,0.2,0.8);
        viewer.data().add_edges(V + PD_min*avg, V - PD_min*avg, blue);
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '9') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete maximal principal curvature:
        // store in K_max_principal
        K_max_principal.setZero(V.rows(), 1);

        // For visualization, better to normalize the range of K_max_principal with the maximal and minimal curvatures
        // store colors in colors_per_vertex
        Eigen::MatrixXd PD_min,PD_max;
        Eigen::VectorXd PV1,PV2;
        igl::principal_curvature(V,F,PD_min,PD_max,PV1,PV2);
        K_max_principal = PV2;

        double mink = K_max_principal.minCoeff();
        double maxk = K_max_principal.maxCoeff();
        for (int i = 0; i < K_max_principal.rows(); ++i)
        {
            igl::jet((K_max_principal(i) - mink) / (maxk - mink), colors_per_vertex(i, 0), colors_per_vertex(i, 1), colors_per_vertex(i, 2));
        }

        // Uncomment the code below to draw a red segment parallel to the maximal curvature direction
        
        const double avg = igl::avg_edge_length(V,F);
        Eigen::RowVector3d red(0.8,0.2,0.2);
        viewer.data().add_edges(V + PD_max*avg, V - PD_max*avg, red);
        
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == 'E') {
        // Add your code for computing explicit Laplacian smoothing here:
        // store the smoothed vertices in V_expLap
        V_expLap = V;

        int iter = 1;
        for (int i = 0; i < iter; ++i)
        {
            // explicit_smooth();
            explicit_smooth();
        }

        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_expLap, F);
    }

    if (key == 'D'){
        // Implicit smoothing for comparison
        // store the smoothed vertices in V_impLap
        V_impLap = V;

        int iter = 1;
        for (int i = 0; i < iter; ++i)
        {
            implicit_smooth();
        }

        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_impLap, F);
    }

    if (key == 'B') {
        // Add your code for computing bilateral smoothing here:
        // store the smoothed vertices in V_bilateral
        // be care of the sign mistake in the paper
        // use v' = v - n * (sum / normalizer) to update
        V_bilateral = V;

        int iter = 5;
        for (int i = 0; i < iter; ++i)
        {
            igl::per_face_normals(V_bilateral, F, FN);
            area_weight_normal_b();
            igl::adjacency_list(F, VV);
            bilateral_smoothing();
        }
        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_bilateral, F);
    }


    return true;
}

void area_weight_normal_b()
{
    N_area.setZero(V.rows(), 3);
    Eigen::VectorXd area;
    igl::doublearea(V_bilateral, F, area);
    area = area.array() / 2;
    
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j = 0; j < F.cols(); ++j)
        {
            N_area.row(F(i, j)) += FN.row(i) * area(i);
        }
    }
    N_area.rowwise().normalize();
    igl::per_vertex_normals(V, F, VN);
    for (int i = 0; i < N_area.rows(); ++i)
    {
        // cout << N_area.row(i).dot(VN.row(i)) << endl;
        if (N_area.row(i).dot(VN.row(i)) < 0)
        {
            N_area.row(i) = -N_area.row(i);
        }
    }
}

double get_sigma_s(set<int> neighbor, Eigen::RowVector3d normal, int vertex_id)
{
    double offset_sum = 0;
    vector<double> offsets;
    for (auto it = neighbor.begin(); it != neighbor.end(); ++it)
    {
        Eigen::RowVector3d qi = V_bilateral.row(*it);
        double offset = abs(normal.dot(V_bilateral.row(vertex_id) - qi));
        offsets.push_back(offset);
        offset_sum += offset;
    }
    offset_sum /= neighbor.size();
    double sigma_s = 0;
    for (int i = 0; i < offsets.size(); ++i)
    {
        sigma_s += (offsets[i] - offset_sum) * (offsets[i] - offset_sum);
    }
    sigma_s /= offsets.size();
    cout << "sigma_s:   " << sigma_s << endl;
    return sqrt(sigma_s);
}
void bilateral_smoothing()
{
    // const double avg = igl::avg_edge_length(V_bilateral,F);
    // cout << "avg edge length.  " << avg << endl; // 0.00989599
    cout << "smooth" << endl;
    for (int i = 0; i < V_bilateral.rows(); ++i)
    {
        Eigen::RowVector3d normal = N_area.row(i);
        double sum = 0;
        double normalizer = 0;
        double sigma_c = 0.01;
        // double sigma_s = get_sigma_s(neighbor, normal, i);
        // if (sigma_s < 1e-5)
        //     sigma_s = 1e-5;
        double sigma_s = 1e-3;
        for (int id = 0; id < VV[i].size(); ++id)
        {
            Eigen::RowVector3d qi = V_bilateral.row(VV[i][id]);
            double t = (V_bilateral.row(i) - qi).norm();
            double h = normal.dot(V_bilateral.row(i) - qi);
            double wc = exp(-t * t / (2 * sigma_c * sigma_c));
            double ws = exp(-h * h / (2 * sigma_s * sigma_s));
            sum += (wc * ws) * h;
            normalizer += wc * ws;
        }
        V_bilateral.row(i) = V_bilateral.row(i) - normal * (sum / normalizer);
    }
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    if (filename.substr(filename.length() - 4) == ".off")
    {
        igl::readOFF(filename, V, F);
    }
    else if (filename.substr(filename.length() - 4) == ".obj")
    {
        igl::readOBJ(filename, V, F);
    }
    else
    {
        std::cerr << "Extension unknown (must be '.off' or '.obj')\n";
        return false;
    }
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
        filename = std::string(argv[1]);
    }
    else {
        filename = std::string("../data/bumpy-cube.obj");
    }
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    viewer.launch();
}
