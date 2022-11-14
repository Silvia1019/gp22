#include <iostream>
#include <igl/readOFF.h>
#include <igl/readTGF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/sum.h>
#include <igl/speye.h>
#include <igl/cotmatrix.h>
#include <imgui/imgui.h>
/*** insert any libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/column_to_quats.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/readDMAT.h>
#include <igl/deform_skeleton.h>
#include <hedra/line_cylinders.h>
#include <igl/slice_into.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/lbs_matrix.h>
#include <igl/dqs.h>
#include <igl/grad.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #Vx3; unposed vertex in task 8
Eigen::MatrixXd V, V_unpose;
// Joints, # x 3
Eigen::MatrixXd Bone;
// Face array, #Fx3
Eigen::MatrixXi F;
// Skeleton edge array, #Ex2
Eigen::MatrixXi E;

// # parents of point handle indices
Eigen::VectorXi P;

Eigen::MatrixXd RM; // RM: Rotation Matrix from .dmat
Eigen::MatrixXi H_ref; // reference handles, #V (4780) x 1 for hand
vector<vector<int>> H; // computed handle list

typedef std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;
RotationList pose; // quaternions for rotations of last frame
const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);

double anim_t = 1.0;
double anim_t_dir = -0.03;
int frame = 0;
int frame_dir = 1; // for animation of each frame
int skeleton_id = 1; // returned from viewer.append_mesh()
bool skeleton_added = false;

int animation = 2; // task 3: 1 for rotation matrix, 2 for quarternion
int selected = 3; // selected skeleton index, for task 4
int task = 0;

Eigen::MatrixXd W, U, M; // W: harmonic weight, U: updated V, M: returned from lbs_matrix
Eigen::VectorXi handles; // #V * 1
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver;

// for task 7 
Eigen::SparseMatrix<double> G, D_area, G_c, G_f;
Eigen::SparseMatrix<double, RowMajor> GG;
Eigen::VectorXi vc, vf;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> deform_solver;
RotationList Qf;

// for task 8
int context_frame = 0;
int context_frame_dir = 1;
Eigen::MatrixXd AF; // all frames
vector<Eigen::MatrixXd> Sigma(4);
vector<Eigen::MatrixXd> Pj(4);
vector<vector<double>> c(4, vector<double>(4, 0.0));
vector<double> aj(4, 0);

void compute_handle() 
{
    // igl::readDMAT("../data/hand/hand-handles.dmat", H_ref);
    H.clear();
    H.resize(E.rows());
    handles.setZero(V.rows());
    // handle selection
    for (int j = 0; j < V.rows(); ++j)
    {
        handles(j) = -1;
    }
    for (int i = 0; i < E.rows(); ++i) // Geometric selection
    {
        double min_distance = 10000;
        Eigen::RowVector3d middle_point = (Bone.row(E(i, 0)) + Bone.row(E(i, 1))) / 2;
        for (int j = 0; j < V.rows(); ++j)
        {
            double dis = (middle_point - V.row(j)).squaredNorm();
            min_distance = min(min_distance, dis);
        }
        for (int j = 0; j < V.rows(); ++j)
        {
            double dis = (middle_point - V.row(j)).squaredNorm();
            double thres = 2.5 * min_distance;
            if (handles(j) == -1 && dis <= thres)
            {
                H[i].push_back(j);
                handles(j) = i;
            }
        }
    }
}

void compute_context_handle()
{
    igl::readDMAT("../data/context-aware/handles.dmat", H_ref);
    H.clear();
    H.resize(E.rows());
    handles.setZero(V.rows());
    // use referenced handles
    for (int i = 0; i < H_ref.rows(); ++i)
    {
        if (H_ref(i, 0) != -1)
        {
            H[H_ref(i, 0)].push_back(i);
        }
        handles(i) = H_ref(i, 0);
    }
}

void visualize_handle(Viewer& viewer)
{
    // Per-vertex color array, #Vx3
    Eigen::MatrixXd colors_per_vertex;
    colors_per_vertex.setZero(V.rows(),3);
    for (int i = 0; i < V.rows(); ++i)
    {
        colors_per_vertex.row(i) = Eigen::RowVector3d(0,0,205./255.);
    }
    for (int i = 0; i < H[selected].size(); ++i)
    {
        colors_per_vertex.row(H[selected][i]) = Eigen::RowVector3d(0.9,0.9,0);;
    }
    viewer.data().set_colors(colors_per_vertex);
}

void compute_weight(Viewer& viewer)
{
    // sovle Lwk = 0, similar to assignment 5
    SparseMatrix<double> L, A, B, lhs;
    igl::cotmatrix(V, F, L);
    W.setZero(V.rows(), E.rows());
    Eigen::VectorXi Non_handles, selected_handles; // index for handles and free vertice
    int handle_num = 0;
    for (int j = 0; j < V.rows(); ++j)
    {
        if (handles(j) != -1)
            handle_num++;
    }
    Non_handles.setZero(V.rows() - handle_num);
    selected_handles.setZero(handle_num);
    int non_handle_index = 0, handle_index = 0;
    for (int j = 0; j < V.rows(); ++j)
    {
        if (handles(j) == -1)
        {
            Non_handles(non_handle_index) = j;
            non_handle_index++;
        }
        else
        {
            selected_handles(handle_index) = j;
            handle_index++;
        }
    }
    igl::slice(L, Non_handles, Non_handles, lhs);
    igl::slice(L, Non_handles, selected_handles, A);
    solver.compute(lhs);
    for (int i = 0; i < E.rows(); ++i)
    {
        Eigen::VectorXd c, rhs, x, d;
        c.setZero(V.rows());
        for (int j = 0; j < H[i].size(); ++j) // wk = 1 for Hk, others = 0
        {
            c(H[i][j]) = 1;
        }
        igl::slice(c, selected_handles, 1, d);
        rhs = -A * d;
        x = solver.solve(rhs);
        igl::slice_into(x, Non_handles, 1, c);
        W.col(i) << c;
    }
    igl::normalize_row_sums(W,W);
    viewer.data().set_data(W.col(selected));
}

void compute_face_quaternion(RotationList &vQ)
{
    RotationList logQ(vQ.size());
    for (int i = 0; i < vQ.size(); ++i)
    {
        double exp_w = vQ[i].norm();
        double w = log(exp_w);
        double a = acos(vQ[i].w() / exp_w);
        if (vQ[i].vec().norm() == 0)
        {
            logQ[i] = Eigen::Quaterniond(w, 0.0, 0.0, 0.0);
            continue;
        }
        logQ[i].w() = w;
        logQ[i].vec() = vQ[i].vec() / vQ[i].vec().norm() * a;
    }
    Qf.resize(F.rows());
    for (int i = 0; i < Qf.size(); ++i)
    {
        Eigen::Quaterniond sum(0.0, 0.0, 0.0, 0.0);
        for (int j = 0; j < vQ.size(); ++j) // 1-k
        {
            double hkf = (W(F(i, 0), j) + W(F(i, 1), j) + W(F(i, 2), j)) / 3.0;
            sum.w() = sum.w() + hkf * logQ[j].w();
            sum.vec() += hkf * logQ[j].vec();
        }
        double v_norm = sum.vec().norm();
        double exp_a = exp(sum.w());
        if (v_norm == 0)
        {
            Qf[i] = Eigen::Quaterniond(exp_a * cos(v_norm), 0.0, 0.0, 0.0);
            continue;
        }
        Qf[i].w() = exp_a * cos(v_norm);
        Qf[i].vec() = exp_a * sin(v_norm) * sum.vec() / v_norm;
    }
}

void poission_stitching(Viewer& viewer)
{
    // compute rhs Q
    Eigen::MatrixXd G_q; // face deformation gradient matrix
    G_q.setZero(3 * F.rows(), 3);
    for (int i = 0; i < F.rows(); ++i)
    {
        G_q.block(i * 3, 0, 3, 3) = Qf[i].toRotationMatrix().transpose();
    }
    // solve linear system
    Eigen::MatrixXd Rhs, R_f, V_c, V_f2;
    Rhs = GG.transpose() * D_area * G_q;
    igl::slice(Rhs, vf, 1, R_f);
    igl::slice(V, vc, 1, V_c); // vi^l = vi^0
    R_f = R_f - G_c * V_c;
    V_f2 = deform_solver.solve(R_f);

    igl::slice_into(V_f2, vf, 1, U);
    igl::slice_into(V_c, vc, 1, U);
    viewer.data().set_vertices(U);
}

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
    // refer to tutorial 403
    // cout << "pre_draw" << endl;
    if(viewer.core().is_animating)
    {
        // Interpolate pose and identity
        RotationList anim_pose(pose.size());
        if (animation == 2)
        { 
            for(int e = 0;e < pose.size(); e++)
            {
                anim_pose[e] = pose[e].slerp(anim_t, Quaterniond::Identity());
            }
        }
        else if (animation == 1)
        {
            int start_position = frame * pose.size() * 3;
            for(int e = 0;e < pose.size(); e++)
            {
                Eigen::Matrix3d rotation;
                rotation = RM.block(start_position + 3 * e, 0, 3, 3);
                anim_pose[e] = rotation; // convert rotation matrix to quaternions
            }
        }
        // Propagate relative rotations via FK to retrieve absolute transformations
        RotationList vQ;
        vector<Vector3d> vT;
        igl::forward_kinematics(Bone, E, P, anim_pose, vQ, vT);
        const int dim = Bone.cols();
        MatrixXd T(E.rows() * (dim + 1), dim);
        for(int e = 0;e < E.rows(); e++)
        {
            Affine3d a = Affine3d::Identity();
            a.translate(vT[e]);
            a.rotate(vQ[e]);
            T.block(e*(dim+1),0,dim+1,dim) =
                a.matrix().transpose().block(0,0,dim+1,dim);
        }
        if (task == 3)
        {
            // Also deform skeleton edges
            MatrixXd CT;
            MatrixXi BET;
            igl::deform_skeleton(Bone, E, T, CT, BET);
            viewer.data().clear();
            viewer.data_list[skeleton_id].set_points(CT, Eigen::RowVector3d(1,0,0));
            
            // viewer.data().set_edges(CT, BET, sea_green);
            Eigen::MatrixXd faceCenter1, faceCenter2, normalColors, VNormals, CNormals;
            Eigen::MatrixXi TNormals;
            normalColors.resize(E.rows(), 3);
            faceCenter1.resize(E.rows(), 3);
            faceCenter2.resize(E.rows(), 3);
            for (int i = 0; i < normalColors.rows(); ++i)
            {
                normalColors.row(i) << sea_green;
                faceCenter1.row(i) << CT.row(BET(i, 0));
                faceCenter2.row(i) << CT.row(BET(i, 1));
            }
            hedra::line_cylinders(faceCenter1, faceCenter2, 0.02, normalColors, 10, VNormals, TNormals, CNormals);
            viewer.data_list[skeleton_id].set_mesh(VNormals, TNormals);
            viewer.data_list[skeleton_id].set_colors(CNormals);
            viewer.data_list[skeleton_id].show_lines=false;
        }
        else if (task == 5)
        {
            U = M*T;
            viewer.data().set_vertices(U);
        }
        else if (task == 6)
        {
            igl::dqs(V,W,vQ,vT,U);
            viewer.data().set_vertices(U);
        }
        else if (task == 7)
        {
            compute_face_quaternion(vQ);
            poission_stitching(viewer);
        }

        viewer.selected_data_index=0;
        viewer.data().compute_normals();
        if (animation == 2)
        {
            anim_t += anim_t_dir;
            anim_t_dir *= (anim_t>=1.0 || anim_t<=0.0?-1.0:1.0);
        }
        else if (animation == 1)
        {
            frame += frame_dir;
            frame_dir *= (frame >= 33 || frame <= 0?-1:1);
        }
    }
    return false;
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

void preprocess_hand(Viewer& viewer)
{
    load_mesh(viewer, "../data/hand/hand.off", V, F);
    igl::readTGF("../data/hand/hand.tgf", Bone, E);
    igl::directed_edge_parents(E, P); // P: -1, 0-20

    Eigen::MatrixXd Q;
    igl::readDMAT("../data/hand/hand-pose_matrot.dmat", RM); // row: 3 x 20 x 34
    igl::readDMAT("../data/hand/hand-pose_quat.dmat", Q); // row: 80
    igl::column_to_quats(Q, pose); // Read pose as matrix of quaternions per row
    assert(pose.size() == E.rows());
}

void compute_unpose(Eigen::MatrixXd &P_j)
{
    // compute rotation and translation as before
    RotationList anim_pose(E.rows());
    for(int e = 0;e < E.rows(); e++)
    {
        Eigen::Matrix3d rotation;
        rotation = P_j.block(3 * e, 0, 3, 3);
        anim_pose[e] = rotation;
    }        
    RotationList vQ;
    vector<Vector3d> vT;
    igl::forward_kinematics(Bone, E, P, anim_pose, vQ, vT);
    const int dim = Bone.cols();
    MatrixXd T(E.rows() * (dim + 1), dim);
    for(int e = 0;e < E.rows(); e++)
    {
        Affine3d a = Affine3d::Identity();
        a.translate(vT[e]);
        a.rotate(vQ[e]);
        T.block(e*(dim+1),0,dim+1,dim) =
            a.matrix().transpose().block(0,0,dim+1,dim);
    }

    // compute V_unpose
    V_unpose.setZero(V.rows(), V.cols());
    for (int i = 0; i < V.rows(); ++i)
    {
        Eigen::MatrixXd M, A, MT;
        M.setZero(4, 3);
        A.setZero(3, 3);
        for (int j = 0; j < W.cols(); ++j)
        {
            Eigen::MatrixXd R_j;
            R_j = T.block(4 * j, 0, 4, 3);
            M = M + R_j * W(i, j);
        }
        Eigen::Vector3d v_i, v_i_unpose;
        v_i << V(i, 0), V(i, 1), V(i, 2); // v_i = M.transpose() * (vi_unpose, 1)^T;
        MT = M.transpose();
        A << MT.col(0), MT.col(1), MT.col(2);
        v_i_unpose = A.colPivHouseholderQr().solve(v_i - MT.col(3));
        V_unpose.row(i) << v_i_unpose(0), v_i_unpose(1), v_i_unpose(2);
    }
}

void compute_weight_a(Eigen::MatrixXd &PP)
{
    double sum = 0;
    for (int i = 0; i < 4; ++i)
    {
        aj[i] = 0;
        for (int j = 0; j < 4; ++j)
        {
            aj[i] += exp(-(Pj[j] - PP).squaredNorm()) * c[i][j];
        }
        sum += aj[i];
    }
    for (int i = 0; i < 4; ++i)
    {
        aj[i] /= sum;
    }
}

bool pre_draw_context(igl::opengl::glfw::Viewer & viewer)
{
    if(viewer.core().is_animating)
    {
        int start_position = context_frame * E.rows() * 3;
        // cout << start_position << endl;
        Eigen::MatrixXd PP;
        PP = AF.block(start_position, 0, E.rows() * 3, 3);
        RotationList anim_pose(E.rows());
        for(int e = 0;e < E.rows(); e++)
        {
            Eigen::Matrix3d rotation;
            rotation = PP.block(3 * e, 0, 3, 3);
            anim_pose[e] = rotation;
        }
        RotationList vQ;
        vector<Vector3d> vT;
        igl::forward_kinematics(Bone, E, P, anim_pose, vQ, vT);
        const int dim = Bone.cols();
        MatrixXd T(E.rows() * (dim + 1), dim);
        for(int e = 0;e < E.rows(); e++)
        {
            Affine3d a = Affine3d::Identity();
            a.translate(vT[e]);
            a.rotate(vQ[e]);
            T.block(e*(dim+1),0,dim+1,dim) =
                a.matrix().transpose().block(0,0,dim+1,dim);
        }
        compute_weight_a(PP);
        for (int i = 0; i < V.rows(); ++i)
        {
            Eigen::MatrixXd M;
            Eigen::VectorXd v_i, v_0, s_i;
            M.setZero(4, 3);
            for (int j = 0; j < W.cols(); ++j)
            {
                Eigen::MatrixXd R_j;
                R_j = T.block(4 * j, 0, 4, 3);
                M = M + R_j * W(i, j);
            }
            v_0.setZero(3);
            v_0 << V(i, 0), V(i, 1), V(i, 2);
            for (int j = 0; j < 4; ++j)
            {
                s_i.setZero(3);
                s_i << Sigma[j](i, 0), Sigma[j](i, 1), Sigma[j](i, 2);
                v_0 += aj[j] * s_i;
            }
            v_i.setZero(4);
            v_i << v_0(0), v_0(1), v_0(2), 1;
            v_i = M.transpose() * v_i;
            U.row(i) << v_i(0), v_i(1), v_i(2);
        }
        viewer.data().set_vertices(U);
        viewer.selected_data_index=0;
        viewer.data().compute_normals();
        context_frame += context_frame_dir;
        context_frame_dir *= (context_frame >= 393 || context_frame <= 0?-1:1);
    }
    return false;
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '3') {
        viewer.data().clear();
        task = 3;
        viewer.callback_pre_draw = &pre_draw;
        viewer.core().is_animating = true;

        preprocess_hand(viewer);
        if (skeleton_added == false)
        {
            skeleton_added = true;
            skeleton_id = viewer.append_mesh();
        }
    }
    if (key == '4') {
        task = 4;
        viewer.core().is_animating = false;
        preprocess_hand(viewer);
        if (skeleton_added == true)
        {
            viewer.data_list[skeleton_id].clear();
        }
        // task 4
        compute_handle();
        // visualize_handle(viewer);
        compute_weight(viewer);
    }
    if (key == '5') {
        task = 5;
        preprocess_hand(viewer);
        if (skeleton_added == true)
        {
            viewer.data_list[skeleton_id].clear();
        }

        compute_handle();
        compute_weight(viewer);
        // precompute linear blend skinning matrix
        igl::lbs_matrix(V,W,M);
        U = V;
        viewer.data().clear();
        viewer.data().set_mesh(U, F);
        viewer.callback_pre_draw = &pre_draw;
        viewer.core().is_animating = true;
    }
    if (key == '6') {
        task = 6;
        preprocess_hand(viewer);
        if (skeleton_added == true)
        {
            viewer.data_list[skeleton_id].clear();
        }

        compute_handle();
        compute_weight(viewer);
        U = V;
        viewer.data().clear();
        viewer.data().set_mesh(U, F);
        viewer.callback_pre_draw = &pre_draw;
        viewer.core().is_animating = true;
    }
    if (key == '7')
    {
        task = 7;
        preprocess_hand(viewer);
        if (skeleton_added == true)
        {
            viewer.data_list[skeleton_id].clear();
        }

        compute_handle();
        compute_weight(viewer);
        U = V;
        viewer.data().clear();
        viewer.data().set_mesh(U, F);
        viewer.callback_pre_draw = &pre_draw;
        viewer.core().is_animating = true;

        // compute D
        Eigen::VectorXd area_vec;
        igl::doublearea(V, F, area_vec);
        VectorXd area_vec_3(area_vec.size() * 3);
        for (int i = 0; i < area_vec.size(); ++i)
        {
        area_vec(i) /= 2;
        area_vec_3(3 * i) = area_vec(i);
        area_vec_3(3 * i + 1) = area_vec(i);
        area_vec_3(3 * i + 2) = area_vec(i);
        }
        igl::diag(area_vec_3, D_area);

        // compute G
        igl::grad(V, F, G);
        // shuffle G
        GG.resize(G.rows(), G.cols());
        for (int i = 0; i < F.rows(); ++i)
        {
            GG.row(3 * i) = G.row(i);
            GG.row(3 * i + 1) = G.row(i + F.rows());
            GG.row(3 * i + 2) = G.row(i + F.rows() * 2);
        }    

        // boundary
        vc.setZero(H[0].size());
        vf.setZero(V.rows() - H[0].size());
        int vf_index = 0;
        for (int i = 0; i < H[0].size(); ++i) // compute H1 vertice list
        {
            vc(i) = H[0][i];
        }
        for (int i = 0; i < V.rows(); ++i) // compute free vertice list
        {
            if (handles(i) != 0)
            {
                vf(vf_index) = i;
                vf_index++;
            }
        }           
        Eigen::SparseMatrix<double> Lhs;            
        Lhs = GG.transpose() * D_area * GG;   
        igl::slice(Lhs, vf, vc, G_c);
        igl::slice(Lhs, vf, vf, G_f);
        deform_solver.compute(G_f);
    }
    if (key == '9')
    {
        task = 8;
        Eigen::MatrixXd V_0;
        load_mesh(viewer, "../data/context-aware/reference.obj", V_0, F);
        igl::readDMAT("../data/context-aware/all_frames.dmat", AF); // 394 x 84 
        for (int i = 0; i < 4; ++i) // 8.1.2
        {
            string path = "../data/context-aware/eg" + to_string(i);
            // cout << path << endl;
            igl::readOBJ(path + ".obj", V, F);

            Eigen::MatrixXd P_j;
            igl::readDMAT(path + ".dmat", P_j);

            igl::readTGF("../data/context-aware/skeleton.tgf", Bone, E);
            igl::directed_edge_parents(E, P); 
            compute_context_handle();
            compute_weight(viewer);

            compute_unpose(P_j);
            Sigma[i].resize(V_unpose.rows(), V_unpose.cols());
            Sigma[i] = V_unpose - V_0;
        }
        for (int i = 0; i < 4; ++i)
        {
            string path = "../data/context-aware/eg" + to_string(i);
            igl::readDMAT(path + ".dmat", Pj[i]);
        }
        for (int i = 0; i < 4; ++i) // compute ci,j
        {
            Eigen::Matrix4d A;
            for (int j = 0; j < 4; ++j)
            {
                for (int k = 0; k < 4; ++k)
                {
                    A(j, k) = exp(-(Pj[j] - Pj[k]).squaredNorm());
                }
            }
            Eigen::Vector4d b, x;
            b.setZero(4);
            b(i) = 1;
            x = A.colPivHouseholderQr().solve(b);
            c[i][0] = x(0);
            c[i][1] = x(1);
            c[i][2] = x(2);
            c[i][3] = x(3);
        }

        load_mesh(viewer, "../data/context-aware/reference.obj", V, F);
        U = V;
        viewer.data().clear();
        viewer.data().set_mesh(U, F);
        viewer.callback_pre_draw = &pre_draw_context;
        viewer.core().is_animating = true;
    }
    if (key == '8')
    {
        task = 8;
        if (skeleton_added == true)
        {
            viewer.data_list[skeleton_id].clear();
        }
        viewer.core().is_animating = false;
        load_mesh(viewer, "../data/context-aware/reference.obj", V, F); // #V 18453, #F 36902
        igl::readOBJ("../data/context-aware/eg3.obj", V, F);

        Eigen::MatrixXd P_j;
        igl::readDMAT("../data/context-aware/eg3.dmat", P_j);
        // cout << P_j.rows() << "   " << P_j.cols() << endl; // 84 x 3 (rotation matrix)

        igl::readTGF("../data/context-aware/skeleton.tgf", Bone, E);
        igl::directed_edge_parents(E, P); 
        compute_context_handle();
        compute_weight(viewer);

        compute_unpose(P_j);
        viewer.data().clear();
        viewer.data().set_mesh(V_unpose,F);
    }
    if (key == '.') // for task 4
    {
        selected++;
    }
    if (key == ',')
    {
        selected--;
    }

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
        filename = std::string("../data/hand/hand.off");
    }
    load_mesh(viewer,filename,V,F);

    viewer.data().line_width = 1;
    viewer.callback_pre_draw = &pre_draw;
    viewer.core().is_animating = false;

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();
	};
    cout << "Press 3-8 for corresponding task" << endl;
    cout << "When in task 4, press '.' to increase the index of selected joint, press ',' to decrease" << endl;
    
    viewer.launch();
}
