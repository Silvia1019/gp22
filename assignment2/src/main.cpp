#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>

#include <igl/Timer.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Normals evaluated via PCA method, #P x3
Eigen::MatrixXd NP;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 32;

// Parameter: grid resolution
int resolution = 40;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

double epsilon = 0.0;

double enlarge = 1.1;

bool isRotate = 0;

string str;

Eigen::VectorXi neighbors;

vector<vector<vector<vector<int>>>> spatial_index;
vector<vector<vector<vector<int>>>> spatial_index_P;
// Functions
void createGrid();
void evaluateImplicitFunc();
void evaluateImplicitFunc_PolygonSoup();
void getLines();
void pcaNormal();
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);

void rotatePointNormal();
// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid()
{
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d dim = (bb_max - bb_min) * enlarge;
    Eigen::RowVector3d extra = dim - (bb_max - bb_min);

    // Grid spacing
    const double dx = dim[0] / (double)(resolution - 1);
    const double dy = dim[1] / (double)(resolution - 1);
    const double dz = dim[2] / (double)(resolution - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min - extra * 0.5 + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc()
{
    // Sphere center
    auto bb_min = grid_points.colwise().minCoeff().eval();
    auto bb_max = grid_points.colwise().maxCoeff().eval();
    Eigen::RowVector3d center = 0.5 * (bb_min + bb_max);

    double radius = 0.5 * (bb_max - bb_min).minCoeff();

    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);

                // Value at (x,y,z) = implicit function for the sphere
                grid_values[index] = (grid_points.row(index) - center).norm() - radius;
            }
        }
    }
}

void computEpsilon()
{
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d diag = bb_max - bb_min;
    epsilon = 0.01 * diag.norm();
}

void buildSpatialIndex(Eigen::MatrixXd points)
{
    if (epsilon == 0.0)
        computEpsilon();
    cout << "epsilon " << epsilon << endl;

    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d diag = bb_max - bb_min;

    int x_dim = diag(0) / epsilon + 2;
    int y_dim = diag(1) / epsilon + 2;
    int z_dim = diag(2) / epsilon + 2;
    cout << "x: " << x_dim << "  y:  " << y_dim << " z:  " << z_dim << endl;
    spatial_index.resize(x_dim, vector<vector<vector<int>>> (y_dim, vector<vector<int>> (z_dim)));
    cout << spatial_index.size() << "  " << spatial_index[0].size() << "  " << spatial_index[0][0].size() << endl;

    // cout << "points rows  " << points.rows() << endl;
    for (int i = 0; i < points.rows(); ++i)
    {
        Eigen::RowVector3d offset = points.row(i) - bb_min;
        int ix = offset(0) / epsilon;
        int iy = offset(1) / epsilon;
        int iz = offset(2) / epsilon;
        // cout << "ix:  " << ix << "  iy:  " << iy << "  iz:   " << iz << endl;
        
        spatial_index[ix][iy][iz].push_back(i);
    }
}

void buildSpatialIndexP()
{
    if (epsilon == 0.0)
        computEpsilon();
    cout << "epsilon " << epsilon << endl;

    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d diag = bb_max - bb_min;

    int x_dim = diag(0) / epsilon + 2;
    int y_dim = diag(1) / epsilon + 2;
    int z_dim = diag(2) / epsilon + 2;
    cout << "x: " << x_dim << "  y:  " << y_dim << " z:  " << z_dim << endl;
    spatial_index_P.resize(x_dim, vector<vector<vector<int>>> (y_dim, vector<vector<int>> (z_dim)));
    cout << spatial_index_P.size() << "  " << spatial_index_P[0].size() << "  " << spatial_index_P[0][0].size() << endl;

    // cout << "points rows  " << points.rows() << endl;
    for (int i = 0; i < P.rows(); ++i)
    {
        Eigen::RowVector3d offset = P.row(i) - bb_min;
        int ix = offset(0) / epsilon;
        int iy = offset(1) / epsilon;
        int iz = offset(2) / epsilon;
        // cout << "ix:  " << ix << "  iy:  " << iy << "  iz:   " << iz << endl;
        
        spatial_index_P[ix][iy][iz].push_back(i);
    }
}

void getNeighbors(Eigen::RowVector3d point)
{
    neighbors.setZero(constrained_points.rows());
    Eigen::RowVector3d bb_min = P.colwise().minCoeff();
    Eigen::RowVector3d offset = point - bb_min;
    int ix = offset(0) / epsilon;
    int iy = offset(1) / epsilon;
    int iz = offset(2) / epsilon;
    // cout << " ix : " << ix << " " << iy << " " << iz << endl;
    int search = wendlandRadius / epsilon + 1;
    // cout << "search. " << search << endl;
    int index = 0;
    for (int i = ix - search; i <= ix + search; ++i)
        for (int j = iy - search; j <= iy + search; ++j)
            for (int k = iz - search; k <= iz + search; ++k)
            {
                if (i < 0 || j < 0 || k < 0 || i >= spatial_index.size() || j >= spatial_index[0].size() || k >= spatial_index[0][0].size())
                    continue;
                for (int t = 0; t < spatial_index[i][j][k].size(); ++t)
                {
                    if ((constrained_points.row(spatial_index[i][j][k][t]) - point).norm() < wendlandRadius)
                    {
                        neighbors(index) = spatial_index[i][j][k][t];
                        index++;
                    }
                }
            }
    neighbors.conservativeResize(index);
    if (index != neighbors.size())
        cout << "neighbors size:   " << neighbors.size() << "  not equal to index: " << index << endl;
}

void getNeighborsbruteforce(Eigen::RowVector3d point)
{
    neighbors.setZero(constrained_points.rows());
    int index = 0;
    for (int i = 0; i < constrained_points.rows(); ++i)
    {
        if ((constrained_points.row(i) - point).norm() < wendlandRadius)
        {
            neighbors(index) = i;
            index++;
        }
    }
    neighbors.conservativeResize(index);
}

void getNeighborsP(Eigen::RowVector3d point, double r = wendlandRadius)
{
    neighbors.setZero(P.rows());
    Eigen::RowVector3d bb_min = P.colwise().minCoeff();
    Eigen::RowVector3d offset = point - bb_min;
    int ix = offset(0) / epsilon;
    int iy = offset(1) / epsilon;
    int iz = offset(2) / epsilon;
    // cout << " ix : " << ix << " " << iy << " " << iz << endl;
    int search = r / epsilon + 1;
    // cout << "search. " << search << endl;
    int index = 0;
    for (int i = ix - search; i <= ix + search; ++i)
        for (int j = iy - search; j <= iy + search; ++j)
            for (int k = iz - search; k <= iz + search; ++k)
            {
                if (i < 0 || j < 0 || k < 0 || i >= spatial_index_P.size() || j >= spatial_index_P[0].size() || k >= spatial_index_P[0][0].size())
                    continue;
                for (int t = 0; t < spatial_index_P[i][j][k].size(); ++t)
                {
                    if ((P.row(spatial_index_P[i][j][k][t]) - point).norm() < r)
                    {
                        neighbors(index) = spatial_index_P[i][j][k][t];
                        index++;
                    }
                }
            }
    neighbors.conservativeResize(index);
    if (index != neighbors.size())
        cout << "neighbors size:   " << neighbors.size() << "  not equal to index: " << index << endl;
}

double wendlandWeight(double r)
{
    return pow(1 - r / wendlandRadius, 4) * (4 * r / wendlandRadius + 1);
}

void evaluateImplicitFuncMLS()
{
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);

                Eigen::RowVector3d gridpoint = grid_points.row(index);
                getNeighbors(gridpoint);
                if (neighbors.size() == 0)
                {
                    grid_values[index] = 10000;
                }
                else
                {
                    // cout << "neighbors size" << neighbors.size() << endl;
                    Eigen::MatrixXd weight;
                    Eigen::MatrixXd basis;
                    weight.setZero(neighbors.size(), neighbors.size());
                    Eigen::MatrixXd values;
                    values.setZero(neighbors.size(), 1);
                    for (int i = 0; i < neighbors.size(); ++i)
                    {
                        Eigen::RowVector3d neighborpoint = constrained_points.row(neighbors(i));
                        double neighborvalue = constrained_values(neighbors(i));
                        double r = (neighborpoint - gridpoint).norm();

                        weight(i, i) = wendlandWeight(r);
                        values.row(i) << neighborvalue;
                        switch (polyDegree)
                        {
                            case 0:
                                basis.resize(neighbors.size(), 1);
                                basis.row(i) << 1;
                                break;
                            case 1:
                                basis.resize(neighbors.size(), 4);
                                basis.row(i) << 1, neighborpoint(0), neighborpoint(1), neighborpoint(2);
                                break;
                            case 2:
                                basis.resize(neighbors.size(), 10);
                                basis.row(i) << 1, neighborpoint(0), neighborpoint(1), neighborpoint(2), neighborpoint(0) * neighborpoint(1), 
                                                neighborpoint(0) * neighborpoint(2), neighborpoint(1) * neighborpoint(2), neighborpoint(0) * neighborpoint(0),
                                                neighborpoint(1) * neighborpoint(1), neighborpoint(2) * neighborpoint(2);
                                break;
                            default:
                                break;
                        }
                    }
                    Eigen::MatrixXd A = weight * basis;
                    Eigen::MatrixXd b = weight * values;
                    Eigen::VectorXd coef = A.colPivHouseholderQr().solve(b);
                    Eigen::VectorXd base;
                    switch (polyDegree)
                    {
                        case 0:
                            base.resize(1);
                            base << 1;
                            break;
                        case 1:
                            base.resize(4);
                            base << 1, gridpoint(0), gridpoint(1), gridpoint(2);
                            break;
                        case 2:
                            base.resize(10);
                            base << 1, gridpoint(0), gridpoint(1), gridpoint(2), gridpoint(0) * gridpoint(1), 
                                            gridpoint(0) * gridpoint(2), gridpoint(1) * gridpoint(2), gridpoint(0) * gridpoint(0),
                                            gridpoint(1) * gridpoint(1), gridpoint(2) * gridpoint(2);
                            break;
                        default:
                            break;
                    }
                    grid_values[index] = base.dot(coef);
                }
                // cout << grid_values[index] << endl;
            }
        }
    }
}

void evaluateImplicitFunc_PolygonSoup()
{
    // Replace with your code here, for "key == '5'"
    // evaluateImplicitFunc();
    buildSpatialIndexP();

    grid_values.resize(resolution * resolution * resolution);

    // Evaluate signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);

                Eigen::RowVector3d gridpoint = grid_points.row(index);
                getNeighborsP(gridpoint);
                if (neighbors.size() == 0)
                {
                    grid_values[index] = 10000;
                }
                else
                {
                    Eigen::MatrixXd weight;
                    Eigen::MatrixXd basis;
                    weight.setZero(neighbors.size(), neighbors.size());
                    Eigen::MatrixXd values;
                    values.setZero(neighbors.size(), 1);
                    for (int i = 0; i < neighbors.size(); ++i)
                    {
                        Eigen::RowVector3d neighborpoint = P.row(neighbors(i)); // P
                        double neighborvalue = (gridpoint - neighborpoint).dot(N.row(neighbors(i))); // replace with function value
                        double r = (neighborpoint - gridpoint).norm();
                        weight(i, i) = wendlandWeight(r);
                        values.row(i) << neighborvalue;
                        // cout << "weight  " << weight(i, i) << "  value. " << values.row(i) << endl;
                        switch (polyDegree)
                        {
                            case 0:
                                basis.resize(neighbors.size(), 1);
                                basis.row(i) << 1;
                                break;
                            case 1:
                                basis.resize(neighbors.size(), 4);
                                basis.row(i) << 1, neighborpoint(0), neighborpoint(1), neighborpoint(2);
                                break;
                            case 2:
                                basis.resize(neighbors.size(), 10);
                                basis.row(i) << 1, neighborpoint(0), neighborpoint(1), neighborpoint(2), neighborpoint(0) * neighborpoint(1), 
                                                neighborpoint(0) * neighborpoint(2), neighborpoint(1) * neighborpoint(2), neighborpoint(0) * neighborpoint(0),
                                                neighborpoint(1) * neighborpoint(1), neighborpoint(2) * neighborpoint(2);
                                break;
                            default:
                                break;
                        }
                    }
                    Eigen::MatrixXd A = weight * basis;
                    Eigen::MatrixXd b = weight * values;
                    Eigen::VectorXd coef = A.colPivHouseholderQr().solve(b);
                    Eigen::VectorXd base;
                    switch (polyDegree)
                    {
                        case 0:
                            base.resize(1);
                            base << 1;
                            break;
                        case 1:
                            base.resize(4);
                            base << 1, gridpoint(0), gridpoint(1), gridpoint(2);
                            break;
                        case 2:
                            base.resize(10);
                            base << 1, gridpoint(0), gridpoint(1), gridpoint(2), gridpoint(0) * gridpoint(1), 
                                            gridpoint(0) * gridpoint(2), gridpoint(1) * gridpoint(2), gridpoint(0) * gridpoint(0),
                                            gridpoint(1) * gridpoint(1), gridpoint(2) * gridpoint(2);
                            break;
                        default:
                            break;
                    }
                    grid_values[index] = base.dot(coef);
                }
                // cout << grid_values[index] << endl;
            }
        }
    }
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines()
{
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1)
                {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1)
                {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1)
                {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

// Estimation of the normals via PCA.
void pcaNormal()
{
    // NP = -N; // to be replaced with your code

    NP.resize(P.rows(), 3);
    buildSpatialIndexP();
    for (int i = 0; i < P.rows(); ++i)
    {
        double r = wendlandRadius;

        getNeighborsP(P.row(i), r);
        while (neighbors.size() < 3)
        {
            r *= 2;
            getNeighborsP(P.row(i), r);
        }
        Eigen::MatrixXd C;
        C.setZero(neighbors.size(), 3);
        for (int j = 0; j < neighbors.size(); ++j)
        {
            C.row(j) << P.row(neighbors(j));
        }
        Eigen::RowVectorXd mean = C.colwise().mean();
        Eigen::MatrixXd aligned = C.rowwise() - mean;
        Eigen::MatrixXd cov = aligned.adjoint() * aligned;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(cov);
        int pos = -1;
        eigensolver.eigenvalues().minCoeff(&pos);
        Eigen::RowVectorXd normal = eigensolver.eigenvectors().col(pos);
        if (normal.dot(N.row(i)) < 0)
            normal = -normal;
        NP.row(i) = normal;
        // cout << "neighbor size. " << neighbors.size() << "  dot product" << normal.dot(N.row(i)) << endl;
    }
}

void rotatePointNormal()
{
    Eigen::RowVectorXd mean = P.colwise().mean();
    Eigen::MatrixXd aligned = P.rowwise() - mean;
    Eigen::MatrixXd cov = aligned.adjoint() * aligned;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(cov);
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
    // cout << "eigenvector:  " << eigenvectors.determinant() << endl;
    if (eigenvectors.determinant() < 0)
        eigenvectors = eigenvectors * -1.0;
    cout << eigenvectors << endl;
    cout << eigenvectors.rightCols(3) << endl;
    Eigen::Vector3d eigenVal = eigensolver.eigenvalues();
    // cout << "eigenvalue:   " << eigenVal << endl;
    P = aligned * eigenvectors;
    N = N * eigenvectors;
}

bool isNearest(int index, Eigen::RowVector3d point)
{
    double dist = (point - P.row(index)).squaredNorm();
    for (int i = 0; i < P.rows(); ++i)
    {
        if ((point - P.row(i)).squaredNorm() < dist)
        {
            return false;
        }
    }
    return true;
}
bool isNearestAccelaration(int index, Eigen::RowVector3d point)
{
    double dist = (point - P.row(index)).squaredNorm();
    Eigen::RowVector3d bb_min = P.colwise().minCoeff();
    Eigen::RowVector3d offset = P.row(index) - bb_min;
    int ix = offset(0) / epsilon;
    int iy = offset(1) / epsilon;
    int iz = offset(2) / epsilon;
    for (int i = ix - 1; i <= ix + 1; ++i)
        for (int j = iy - 1; j <= iy + 1; ++j)
            for (int k = iz - 1; k <= iz + 1; ++k)
            {
                if (i < 0 || j < 0 || k < 0 || i >= spatial_index_P.size() || j >= spatial_index_P[0].size() || k >= spatial_index_P[0][0].size())
                    continue;
                for (int t = 0; t < spatial_index_P[i][j][k].size(); ++t)
                {
                    if (spatial_index_P[i][j][k][t] != index && (P.row(spatial_index_P[i][j][k][t]) - point).squaredNorm() < dist)
                    {
                        return false;
                    }
                }
            }
    return true;
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == '1')
    {
        // Show imported points
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
    }

    if (key == '2')
    {
        // Show all constraints
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        // Add your code for computing auxiliary constraint points here
        // Add code for displaying all points, as above
        constrained_points.resize(P.rows() * 3, 3);
        constrained_values.resize(P.rows() * 3, 1);

        if (epsilon == 0.0)
            computEpsilon();
        double eps_plus = epsilon;
        double eps_minus = epsilon;
        cout << "epsilon:   " << epsilon << endl;
        cout << "P rows:  " << P.rows() << endl;

        igl::Timer timer;
        timer.start();

        buildSpatialIndexP();
        for (int i = 0; i < P.rows(); i++)
        {
            eps_plus = epsilon;
            Eigen::RowVector3d p_plus = P.row(i) + eps_plus * N.row(i);
            // if (isNearest(i, p_plus) != isNearestAccelaration(i, p_plus))
                // cout << "accelaration wrong plus" << endl;
            while(!isNearestAccelaration(i, p_plus))
            {
                eps_plus *= 0.5;
                p_plus = P.row(i) + eps_plus * N.row(i);
            }
            eps_minus = epsilon;
            Eigen::RowVector3d p_minus = P.row(i) - eps_minus * N.row(i);
            // if (isNearest(i, p_minus) != isNearestAccelaration(i, p_minus))
                // cout << "accelaration wrong minus" << endl;
            while(!isNearestAccelaration(i, p_minus))
            {
                eps_minus *= 0.5;
                p_minus = P.row(i) - eps_minus * N.row(i);
            }
            constrained_points.row(i) << P.row(i);
            constrained_points.row(i + P.rows()) << p_plus;
            constrained_points.row(i + P.rows() * 2) << p_minus;

            constrained_values.row(i) << 0;
            constrained_values.row(i + P.rows()) << eps_plus;
            constrained_values.row(i + P.rows() * 2) << -eps_minus;

        }
        timer.stop();
        cout << "constraints take time: " << timer.getElapsedTime() << endl;

        buildSpatialIndex(constrained_points);
        viewer.data().point_size = 11;

		viewer.data().add_points(constrained_points.block(0,0,P.rows(),3), Eigen::RowVector3d(0, 0, 1));
		viewer.data().add_points(constrained_points.block(P.rows(), 0, P.rows(), 3), Eigen::RowVector3d(1, 0, 0));
		viewer.data().add_points(constrained_points.block(2 * P.rows(), 0, P.rows(), 3), Eigen::RowVector3d(0, 1, 0));
    }

    if (key == '3')
    {
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        // Add code for creating a grid
        // Add your code for evaluating the implicit function at the grid points
        // Add code for displaying points and lines
        // You can use the following example:

        /*** begin: sphere example, replace (at least partially) with your code ***/
        // Make grid
        createGrid();

        igl::Timer timermls;
        timermls.start();

        // Evaluate implicit function
        evaluateImplicitFuncMLS();

        timermls.stop();
        cout << "MLS take time: " << timermls.getElapsedTime() << endl;


        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i)
        {
            double value = grid_values(i);
            if (value < 0)
            {
                grid_colors(i, 1) = 1;
            }
            else
            {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/
    }

    if (key == '4')
    {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0))
        {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
        if (V.rows() == 0)
        {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);

        cout << "str. " << str << endl;
        string path = "../res/" + str + ".off";
		igl::writeOFF(path, V, F);
    }

    if (key == '5')
    {
        // Use the structure for key=='3' but replace the function evaluateImplicitFunc();
        // with a function performing the approximation of the implicit surface from polygon soup
        // Ref: Chen Shen, James F. O’Brien, and Jonathan Richard Shewchuk. Interpolating and approximating implicit surfaces from polygon soup.

        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function --> Function to be modified here
        evaluateImplicitFunc_PolygonSoup();

        // get grid lines
        getLines();

        // Display the reconstruction
        callback_key_down(viewer, '4', modifiers);
    }

    if (key == '6' || key == '7' || key == '8')
    {
        // Implement PCA Normal Estimation --> Function to be modified here
        pcaNormal();

        // To use the normals estimated via PCA instead of the input normals and then restaurate the input normals
        Eigen::MatrixXd N_tmp = N;
        N = NP;

        switch (key)
        {
        case '6':
            callback_key_down(viewer, '2', modifiers);
            break;
        case '7':
            callback_key_down(viewer, '3', modifiers);
            break;
        case '8':
            callback_key_down(viewer, '3', modifiers);
            callback_key_down(viewer, '4', modifiers);
            break;
        default:
            break;
        }

        // Restore input normals
        N = N_tmp;
    }

    return true;
}

bool callback_load_mesh(Viewer &viewer, string filename)
{
    igl::readOFF(filename, P, F, N);
    callback_key_down(viewer, '1', 0);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../data/sphere.off", P, F, N);
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1], P, F, N);
    }

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::InputInt("Resolution", &resolution, 0, 0);
            if (ImGui::Button("Reset Grid", ImVec2(-1, 0)))
            {
                std::cout << "ResetGrid\n";
                // Recreate the grid
                createGrid();
                // Switch view to show the grid
                callback_key_down(viewer, '3', 0);
            }
            ImGui::InputInt("polyDegree", &polyDegree, 0, 0);
            ImGui::InputDouble("radius", &wendlandRadius, 0, 0);
            ImGui::InputText("input text", str, IM_ARRAYSIZE(&str));
            // TODO: Add more parameters to tweak here...
        }
    };

    if (isRotate)
    {
        rotatePointNormal();
    }

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
