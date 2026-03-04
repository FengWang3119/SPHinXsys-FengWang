#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <numeric> 
#include <string>

using namespace std;
using Vec = std::vector<double>;

// -------- RANS coefficients --------
double std_kw_beta_star_ = 0.09;
double std_kw_sigma_star_ = 0.6;
double std_kw_alpha_ = 0.52;
double std_kw_sigma_ = 0.5;
double std_kw_f_beta_ = 1.0;
double std_kw_beta_0_ = 0.0708;
double std_kw_sigma_do_ = 0.125;
double std_kw_C_lim_ = 0.875;
double std_kw_beta_i_ = 0.075;
double sigma_d_value = 1.0 / 8.0;

double std_kw_beta_ = std_kw_beta_0_ * std_kw_f_beta_;
double std_kw_beta_star_25_ = pow(std_kw_beta_star_, 0.25);
double std_kw_beta_star_5_ = pow(std_kw_beta_star_, 0.5);

double tiny = 1.0e-6;

double compute_flow_rate_simpson_safe(const std::vector<double>& U, double h)
{
    int N = U.size();

    if (N < 3)
    {
        throw std::runtime_error(
            "Simpson rule requires at least 3 points."
        );
    }

    if (N % 2 == 0)
    {
        // 前 N-1 点用 Simpson
        double sum = U[0] + U[N - 2];

        for (int i = 1; i < N - 2; ++i)
        {
            if (i % 2 == 0)
                sum += 2.0 * U[i];
            else
                sum += 4.0 * U[i];
        }

        double Q_simpson = sum * h / 3.0;

        // 最后一个区间用梯形
        double Q_trap = 0.5 * h * (U[N - 2] + U[N - 1]);

        return Q_simpson + Q_trap;
    }
    else
    {
        double sum = U[0] + U[N - 1];

        for (int i = 1; i < N - 1; ++i)
        {
            if (i % 2 == 0)
                sum += 2.0 * U[i];
            else
                sum += 4.0 * U[i];
        }

        return sum * h / 3.0;
    }
}

// ================= TDMA =================
// Solve a tridiagonal system: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
// a[0] must be 0, c[n-1] will be ignored
Vec tdma(const Vec& a, const Vec& b, const Vec& c, const Vec& d) {
    int n = d.size();
    if (a.size() != n || b.size() != n || c.size() != n) {
        throw std::invalid_argument("TDMA: vector size mismatch!");
    }

    Vec cp(n, 0.0); // modified upper diagonal
    Vec dp(n, 0.0); // modified right-hand side
    Vec x(n, 0.0);  // solution vector

    // ---------------- Step 0: first row ----------------
    if (std::abs(b[0]) < 1e-14) throw std::runtime_error("TDMA: b[0] too small!");
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    // ---------------- Forward sweep ----------------
    for (int i = 1; i < n; i++) 
    {
        double denom = b[i] - a[i] * cp[i - 1];
        if (std::abs(denom) < 1e-14) {
            throw std::runtime_error("TDMA: denom too small at row ");
        }
        cp[i] = (i < n - 1) ? c[i] / denom : 0.0;      // last row has no right neighbor
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }

    // ---------------- Back substitution ----------------
    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    return x;
}


void solve_1D_half_channel() 
{

    // -------- Input parameters --------
    double utau_init = 6.37309e-02;
    double H = 2.0;
    double delta = H / 2.0;
    double nu = 3.5e-4;
    double U_avg = 1.0;

    double u_init = 1.0;
    double k_init = 1.0e-5;
    double turbu_omega_init = 2.056;

    double convergence_criteria_outer = 1.0e-6;
        
    double flow_rate_target = U_avg * delta;
    double utau = utau_init;

    //Note: Node arrangement, solving for half channel
    int ny = 16;
    double hy = delta / ny;
    double y_p = 0.5 * hy;
    //================== Input index ==================
    int NF = 2 * ny;
    int index = 2;

    // -------- RANS coefficients --------
    double std_kw_beta_star_ = 0.09;
    double std_kw_sigma_star_ = 0.6;
    double std_kw_alpha_ = 0.52;
    double std_kw_sigma_ = 0.5;
    double std_kw_f_beta_ = 1.0;
    double std_kw_beta_0_ = 0.0708;
    double std_kw_sigma_do_ = 0.125;
    double std_kw_C_lim_ = 0.875;
    double std_kw_beta_i_ = 0.075;
    double sigma_d_value = 1.0 / 8.0;

    double std_kw_beta_ = std_kw_beta_0_ * std_kw_f_beta_;
    double std_kw_beta_star_25_ = pow(std_kw_beta_star_,0.25);
    double std_kw_beta_star_5_ = pow(std_kw_beta_star_,0.5);

    double tiny = 1.0e-6;

    printf("hy=%f\n", hy);
    printf("yp=%f\n", y_p);

    double yplus = utau * y_p / nu;
    printf("yplus=%f\n", yplus);

    double u_p = utau * yplus;
    printf("u_p=%f\n",u_p);
    double turbu_omega_p = 6.0 * nu / (std_kw_beta_i_ * y_p * y_p);

    Vec y(ny);  //down half nodes, actually solved nodes
    for (int i = 0; i < ny; ++i) {
        y[i] = y_p + i * hy;
    }
    Vec y_whole(2*ny);
    int Ny_whole = ny * 2;
    for (int i = 0; i < Ny_whole; ++i) {
        y_whole[i] = y_p + i * hy;
    }
    std::cout << "The number of fluid particle along the Y direction is "
        << Ny_whole << std::endl;

    std::cout << "y_whole = ";
    for (const auto& v : y_whole) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    std::cout << "press to continue" << std::endl;
    std::cin.get();

    Vec u_init_value(ny, u_init);
    Vec k_init_value(ny, k_init);
    Vec turbu_omega_init_value(ny, turbu_omega_init);
    // u_init_value
    std::cout << "u_init_value:" << std::endl;
    for (const auto& v : u_init_value) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    // k_init_value
    std::cout << "\n k_init_value:" << std::endl;
    for (const auto& v : k_init_value) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    // turbu_omega_init_value
    std::cout << "\n turbu_omega_init_value:" << std::endl;
    for (const auto& v : turbu_omega_init_value) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    std::cout << "\n Press Enter to continue..." << std::endl;
    std::cin.get();

    std::vector<double> phi_current(0, 0.0);
    std::vector<double> phi_solved(3 * ny, 0.0);
    int num_iter_out = 0;

    // concatenate
    phi_current.insert(phi_current.end(), u_init_value.begin(), u_init_value.end());
    phi_current.insert(phi_current.end(), k_init_value.begin(), k_init_value.end());
    phi_current.insert(phi_current.end(), turbu_omega_init_value.begin(), turbu_omega_init_value.end());

    // copy
    phi_solved = phi_current;

    // print
    std::cout << "Initial value = ";
    for (const auto& v : phi_current) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    double differ = 1.0;
    while (differ > convergence_criteria_outer)
    {
        std::cout << "Entering while, differ = " << differ
            << ", convergence_criteria_outer = " << convergence_criteria_outer << std::endl;
        // ---update the star values---
        int n_start = 0;
        std::vector<double> u_star(
            phi_solved.begin() + n_start,
            phi_solved.begin() + n_start + ny
        );

        n_start = ny;
        std::vector<double> k_star(
            phi_solved.begin() + n_start,
            phi_solved.begin() + n_start + ny
        );

        n_start = 2 * ny;
        std::vector<double> turbu_omega_star(
            phi_solved.begin() + n_start,
            phi_solved.begin() + n_start + ny
        );

        // initialisation 
        std::vector<double> dkdy(ny, 0.0);
        std::vector<double> dwdy(ny, 0.0);
        // backward diff（from i=1 ）
        for (int i = 1; i < ny; ++i) {
            dkdy[i] = (k_star[i] - k_star[i - 1]) / hy;
            dwdy[i] = (turbu_omega_star[i] - turbu_omega_star[i - 1]) / hy;
        }
        // B.C.
        dkdy[0] = 0.0;
        dwdy[0] = 0.0;

        std::vector<double> dudy_discretized(ny, 0.0);
        for (int i = 1; i < ny; ++i) {
            dudy_discretized[i] = (u_star[i] - u_star[i - 1]) / hy;
        }
        dudy_discretized[0] = 0.0;

        std::vector<double> turbu_omega_tilde(ny);
        std::vector<double> nut_star(ny);

        for (int i = 0; i < ny; ++i) {
            turbu_omega_tilde[i] = std::max(
                turbu_omega_star[i],
                std_kw_C_lim_ * dudy_discretized[i] / std_kw_beta_star_5_
            );
            nut_star[i] = k_star[i] / (turbu_omega_tilde[i] + tiny);
        }

        std::vector<double> dudy_star(ny);
        for (int i = 0; i < ny; ++i) {
            dudy_star[i] = utau * utau * (1.0 - y[i] / delta) / (nu + nut_star[i]);
        }

        // ---update the linearized source term variables---
        // For velocity
        std::vector<double> Sc_u(ny);
        std::vector<double> Sp_u(ny, 0.0);

        for (int i = 0; i < ny; ++i) {
            Sc_u[i] = -1.0 / (nu + nut_star[i]);
            // Sp_u already is 0
        }

        // diffusion coefficient for k
        std::vector<double> diffusion_coefficient_k(ny);
        for (int i = 0; i < ny; ++i) {
            diffusion_coefficient_k[i] = nu + std_kw_sigma_star_ * k_star[i] / (turbu_omega_star[i] + tiny);
        }

        // 后向差分
        std::vector<double> dDkdy(ny, 0.0);
        for (int i = 1; i < ny; ++i) {
            dDkdy[i] = (diffusion_coefficient_k[i] - diffusion_coefficient_k[i - 1]) / hy;
        }

        // 边界
        dDkdy[0] = 0.0;

        // part_extra_viscous_term_k = dDkdy * dkdy
        std::vector<double> part_extra_viscous_term_k(ny);
        std::vector<double> Sc_k(ny);
        std::vector<double> Sp_k(ny);

        for (int i = 0; i < ny; ++i) {
            part_extra_viscous_term_k[i] = dDkdy[i] * dkdy[i];
            Sc_k[i] = (nut_star[i] * dudy_star[i] * dudy_star[i] + part_extra_viscous_term_k[i])
                / diffusion_coefficient_k[i];
            Sp_k[i] = std_kw_beta_star_ * turbu_omega_star[i] / diffusion_coefficient_k[i];
        }

        // diffusion coefficient for turbu_omega
        std::vector<double> diffusion_coefficient_turbu_omega(ny);
        for (int i = 0; i < ny; ++i) {
            diffusion_coefficient_turbu_omega[i] = nu + std_kw_sigma_ * k_star[i] / (turbu_omega_star[i] + tiny);
        }

        // 后向差分
        std::vector<double> dDwdy(ny, 0.0);
        for (int i = 1; i < ny; ++i) {
            dDwdy[i] = (diffusion_coefficient_turbu_omega[i] - diffusion_coefficient_turbu_omega[i - 1]) / hy;
        }
        dDwdy[0] = 0.0;

        // part_extra_viscous_term_w = dDwdy * dwdy
        std::vector<double> part_extra_viscous_term_w(ny);
        for (int i = 0; i < ny; ++i) {
            part_extra_viscous_term_w[i] = dDwdy[i] * dwdy[i];
        }

        // part_production = alpha*omega/k * nut * (dudy)^2
        std::vector<double> part_production(ny);
        for (int i = 0; i < ny; ++i) {
            part_production[i] = std_kw_alpha_ * turbu_omega_star[i] / (k_star[i] + tiny)
                * nut_star[i] * dudy_star[i] * dudy_star[i];
        }

        // grad_prod = dkdy * dwdy
        std::vector<double> grad_prod(ny);
        std::vector<double> part_cross_diffusion(ny);
        for (int i = 0; i < ny; ++i) {
            grad_prod[i] = dkdy[i] * dwdy[i];
            double sigma_d = (grad_prod[i] > 0.0) ? sigma_d_value : 0.0;
            part_cross_diffusion[i] = sigma_d / (turbu_omega_star[i] + tiny) * grad_prod[i];
        }

        // Sc and Sp
        std::vector<double> Sc_turbu_omega(ny);
        std::vector<double> Sp_turbu_omega(ny);
        for (int i = 0; i < ny; ++i) {
            Sc_turbu_omega[i] = (part_production[i] + part_cross_diffusion[i] + part_extra_viscous_term_w[i])
                / diffusion_coefficient_turbu_omega[i];
            Sp_turbu_omega[i] = std_kw_beta_ * turbu_omega_star[i] / diffusion_coefficient_turbu_omega[i];
        }

        // 初始化
        std::vector<double> a_u(ny, 0.0);
        std::vector<double> b_u(ny, 0.0);
        std::vector<double> c_u(ny, 0.0);
        std::vector<double> d_u(ny, 0.0);

        // 循环赋值
        for (int i = 1; i < ny; ++i) {
            a_u[i] = -1.0;
            b_u[i] = 1.0;
            c_u[i] = 0.0;
            d_u[i] = utau * utau * (1.0 - y[i] / delta) * hy * Sc_u[i] * (-1.0);
        }

        // 边界，第一个节点
        a_u[0] = 0.0;
        b_u[0] = 1.0;
        c_u[0] = 0.0;
        d_u[0] = u_p;

        // TDMA 求解
        std::vector<double> U_new = tdma(a_u, b_u, c_u, d_u);

        // 初始化
        std::vector<double> a_k(ny, 0.0);
        std::vector<double> b_k(ny, 0.0);
        std::vector<double> c_k(ny, 0.0);
        std::vector<double> d_k(ny, 0.0);

        // 内部节点 i = 1 ... ny-2
        for (int i = 1; i < ny - 1; ++i) {
            a_k[i] = 1.0;
            b_k[i] = -2.0 - Sp_k[i] * hy * hy;
            c_k[i] = 1.0;
            d_k[i] = -1.0 * hy * hy * Sc_k[i];
        }

        // 边界，第一个节点 grad k = 0
        a_k[0] = 0.0;
        b_k[0] = -1.0 - Sp_k[0] * hy * hy;
        c_k[0] = 1.0;
        d_k[0] = -1.0 * hy * hy * Sc_k[0];

        // 边界，最后一个节点 grad k = 0
        int last = ny - 1;
        a_k[last] = 1.0;
        b_k[last] = -1.0 - Sp_k[last] * hy * hy;
        c_k[last] = 0.0;
        d_k[last] = -1.0 * hy * hy * Sc_k[last];

        // TDMA 求解
        std::vector<double> K_new = tdma(a_k, b_k, c_k, d_k);

        // 初始化
        std::vector<double> a_w(ny, 0.0);
        std::vector<double> b_w(ny, 0.0);
        std::vector<double> c_w(ny, 0.0);
        std::vector<double> d_w(ny, 0.0);

        // 内部节点 i = 1 ... ny-2
        for (int i = 1; i < ny - 1; ++i) {
            a_w[i] = 1.0;
            b_w[i] = -2.0 - Sp_turbu_omega[i] * hy * hy;
            c_w[i] = 1.0;
            d_w[i] = -1.0 * hy * hy * Sc_turbu_omega[i];
        }

        // 边界，第一个节点 w = Wp
        a_w[0] = 0.0;
        b_w[0] = 1.0;
        c_w[0] = 0.0;
        d_w[0] = turbu_omega_p;

        // 边界，最后一个节点 grad w = 0
        last = ny - 1;
        a_w[last] = 1.0;
        b_w[last] = -1.0 - Sp_turbu_omega[last] * hy * hy;
        c_w[last] = 0.0;
        d_w[last] = -1.0 * hy * hy * Sc_turbu_omega[last];

        // TDMA 求解
        std::vector<double> Turbu_omega_new = tdma(a_w, b_w, c_w, d_w);

        // 更新 phi_solved
        n_start = 0;
        for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = U_new[i];

        n_start += ny;
        for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = K_new[i];

        n_start += ny;
        for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = Turbu_omega_new[i];

        // 打印 phi_solved
        std::cout << "phi_solved = ";
        for (const auto& v : phi_solved) std::cout << v << " ";
        std::cout << std::endl;

        // --- Check flow rate ---
        double flow_rate_current = accumulate(U_new.begin(), U_new.end(), 0.0) * hy;

        // --- flow control ---
        double ratio = flow_rate_target / (flow_rate_current + 1e-12);
        double alpha = 0.03;  

        double utau_new = utau * std::sqrt(ratio);
        utau = (1.0 - alpha) * utau + alpha * utau_new;

        // 打印
        std::cout << "flow_rate_current = " << flow_rate_current
            << ", target = " << flow_rate_target << std::endl;
        std::cout << "updated utau = " << utau << std::endl;

        // --- Check convergence ---
        differ = 0.0;
        for (size_t i = 0; i < phi_current.size(); ++i) {
            double diff = phi_current[i] - phi_solved[i];
            differ += diff * diff;
        }
        differ = std::sqrt(differ);

        std::cout << "differ: " << differ << std::endl;

        // 更新 phi_current
        phi_current = phi_solved;

        // 更新迭代次数
        num_iter_out += 1;

        // 打印
        std::cout << "num_iter_out = " << num_iter_out << std::endl;
        std::cout << "------------" << std::endl;

    }
    std::cout << "******Converge******" << std::endl;

    std::cout << "******The results are: ******" << std::endl;

    int n_start = 0;
    std::cout << "U = ";
    for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
    std::cout << std::endl;

    n_start += ny;
    std::cout << "K = ";
    for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
    std::cout << std::endl;

    n_start += ny;
    std::cout << "Turbu_omega = ";
    for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
    std::cout << std::endl;

    // ================== 提取解 ==================
    n_start = 0;

    std::vector<double> U(ny);
    for (int i = 0; i < ny; ++i) U[i] = phi_solved[n_start + i];

    n_start += ny;
    std::vector<double> K(ny);
    for (int i = 0; i < ny; ++i) K[i] = phi_solved[n_start + i];

    n_start += ny;
    std::vector<double> OMEGA(ny);
    for (int i = 0; i < ny; ++i) OMEGA[i] = phi_solved[n_start + i];

    // ================== 使用已有的 y ==================
    std::vector<double> Y = y;

    // ================== 计算湍流粘性系数 ==================
    double eps = 1e-12;
    std::vector<double> NUT(ny);
    for (int i = 0; i < ny; ++i) {
        NUT[i] = K[i] / (OMEGA[i] + eps);
    }

    // ================== 输出 Tecplot 文件 ==================
    std::string header_line = "ZONE T=\"SPH(1D)45 NF="
        + std::to_string(NF)
        + " ("
        + std::to_string(index)
        + ")\"";
    std::string filename = "pipe_kw_nf" + std::to_string(NF) + "_" + std::to_string(index) + ".dat";

    std::ofstream fout(filename);
    fout << "$VARIABLES = \"Y\", \"U\", \"K\", \"OMEGA\", \"NUT(k/omega)\"\n";
    fout << header_line << "\n";

    fout << std::scientific << std::setprecision(8);

    for (int i = 0; i < ny; ++i) {
        fout << Y[i] << " "
            << U[i] << " "
            << K[i] << " "
            << OMEGA[i] << " "
            << NUT[i] << "\n";
    }

    fout.close();

    std::cout << "Tecplot 文件已生成: " << filename << std::endl;
}

void solve_1D_sublayer(double u_p_outer, double k_p_outer, double w_p_outer, double vel_grad_p_outer, double nut_p_outer, double h_sublayer, double Q_target)
{
    //------------------------------------------------↓ Input parameters ↓------------------------------------------------
    double utau_init = 6.37309e-02;
    double H = 2.0;
    double height_sublayer = h_sublayer;
    double nu = 3.5e-4;
    double U_avg = 1.0;

    double u_init = u_p_outer;
    double k_init = k_p_outer;
    double turbu_omega_init = w_p_outer;

    double convergence_criteria_outer = 1.0e-6;

    double flow_rate_target = Q_target;
    double utau = utau_init;
    //------------------------------------------------↑ Input parameters ↑------------------------------------------------

    //------------------------------------------------↓ Node arrangement, for sublayer ↓------------------------------------------------
    int ny = 5;
    double hy = height_sublayer / (double(ny) + 0.5); // distance from node U to P_outer is hy, hence with a 0.5
    double y_p = 0.5 * hy;
    Vec y(ny);  //computational nodes
    for (int i = 0; i < ny; ++i) 
    {
        y[i] = y_p + i * hy;
    }
    //------------------------------------------------↑ Node arrangement, for sublayer ↑------------------------------------------------
    
    //------------------------------------------------↓ Input index ↓------------------------------------------------
    int NF = 2 * ny;
    int index = 28;
    //------------------------------------------------↑ Input index ↑------------------------------------------------

    //------------------------------------------------↓ Calculate P value ↓------------------------------------------------
    double yplus = utau * y_p / nu;
    double u_p = utau * yplus;
    double turbu_omega_p = 6.0 * nu / (std_kw_beta_i_ * y_p * y_p);
    //------------------------------------------------↑ Calculate P value ↑------------------------------------------------

    printf("hy=%f\n", hy);
    printf("yp=%f\n", y_p);
    printf("yplus=%f\n", yplus);
    printf("u_p=%f\n", u_p);
    
    Vec y_whole(2 * ny);
    int Ny_whole = ny * 2;
    for (int i = 0; i < Ny_whole; ++i) 
    {
        y_whole[i] = y_p + i * hy;
    }
    std::cout << "NF= " << Ny_whole << std::endl;
    std::cout << "ny= " << ny << std::endl;
    std::cout << "y = ";
    for (const auto& v : y) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    std::cout << "press to continue" << std::endl;
    std::cin.get();

    //------------------------------------------------↓ Construct initial value ↓------------------------------------------------
    Vec u_init_value(ny, u_init);
    Vec k_init_value(ny, k_init);
    Vec turbu_omega_init_value(ny, turbu_omega_init);
    //------------------------------------------------↑ Construct initial value ↑------------------------------------------------
    
    std::cout << "u_init_value:" << std::endl;
    for (const auto& v : u_init_value) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    std::cout << "\n k_init_value:" << std::endl;
    for (const auto& v : k_init_value) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    std::cout << "\n turbu_omega_init_value:" << std::endl;
    for (const auto& v : turbu_omega_init_value) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    std::cout << "\n Press Enter to continue..." << std::endl;
    std::cin.get();

    //------------------------------------------------↓ Construct solution vector ↓------------------------------------------------
    std::vector<double> phi_current(0, 0.0);
    std::vector<double> phi_solved(3 * ny, 0.0);
    phi_current.insert(phi_current.end(), u_init_value.begin(), u_init_value.end());
    phi_current.insert(phi_current.end(), k_init_value.begin(), k_init_value.end());
    phi_current.insert(phi_current.end(), turbu_omega_init_value.begin(), turbu_omega_init_value.end());
    phi_solved = phi_current;
    //------------------------------------------------↑ Construct solution vector ↑------------------------------------------------

    std::cout << "Initial value = ";
    for (const auto& v : phi_current) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    double differ = 1.0; // Should have a value 
    int num_iter_out = 0;
    int n_start = 0;
    int last = 0;
    int first_index = 0;
    double flow_rate_current = 0.0;
    while (differ > convergence_criteria_outer)
    {
        //------------------------------------------------↓ Update the star values ↓------------------------------------------------
        n_start = 0;
        std::vector<double> u_star(
            phi_solved.begin() + n_start,
            phi_solved.begin() + n_start + ny
        );

        n_start = ny;
        std::vector<double> k_star(
            phi_solved.begin() + n_start,
            phi_solved.begin() + n_start + ny
        );

        n_start = 2 * ny;
        std::vector<double> turbu_omega_star(
            phi_solved.begin() + n_start,
            phi_solved.begin() + n_start + ny
        );
        //------------------------------------------------↑ Update the star values ↑------------------------------------------------

        //------------------------------------------------↓ Calculate Dk, Dw, C_su ↓------------------------------------------------
        std::vector<double> diffusion_coefficient_k(ny);
        std::vector<double> diffusion_coefficient_turbu_omega(ny);
        std::vector<double> C_su(ny);
        double tau_over_rho_outer = (nu + nut_p_outer) * vel_grad_p_outer;
        for (int i = 0; i < ny; ++i) {
            diffusion_coefficient_k[i] = nu + std_kw_sigma_star_ * k_star[i] / (turbu_omega_star[i] + tiny);
            diffusion_coefficient_turbu_omega[i] = nu + std_kw_sigma_ * k_star[i] / (turbu_omega_star[i] + tiny);
            C_su[i] = (utau * utau - tau_over_rho_outer) * (1.0 - y[i] / height_sublayer) + tau_over_rho_outer;
        }
        //------------------------------------------------↑ Calculate Dk, Dw  ↑------------------------------------------------

        //------------------------------------------------↓ Calculate gradients of u, k, omega, Dk, Dw ↓------------------------------------------------
        std::vector<double> dudy_discretized(ny, 0.0);
        std::vector<double> dkdy(ny, 0.0);
        std::vector<double> dwdy(ny, 0.0);
        // central diff（from i=1 to i=ny-2 ）
        for (int i = 1; i < ny - 1; ++i)
        {
            dudy_discretized[i] = (u_star[i + 1] - u_star[i - 1]) / (2.0 * hy);
            dkdy[i] = (k_star[i + 1] - k_star[i - 1]) / (2.0 * hy);
            dwdy[i] = (turbu_omega_star[i + 1] - turbu_omega_star[i - 1]) / (2.0 * hy);
        }
        // B.C. near wall (i=0, central difference)
        dudy_discretized[0] = (u_star[1] + u_star[0]) / (2.0 * hy); // mirror B.C. u_star[-1] = -u_star[0]
        dkdy[0] = (k_star[1] - k_star[0]) / (2.0 * hy); // zero gradient, k_star[-1] = k_star[0]
        dwdy[0] = (turbu_omega_star[1] - turbu_omega_star[0]) / (2.0 * hy); //zero gradient, turbu_omega_star[-1] = turbu_omega[0]
        // B.C. near P_outer (i=ny-1, central difference)
        dudy_discretized[ny - 1] = (u_p_outer - u_star[ny - 2]) / (2.0 * hy);
        dkdy[ny - 1] = (k_p_outer - k_star[ny - 2]) / (2.0 * hy);
        dwdy[ny - 1] = (w_p_outer - turbu_omega_star[ny - 2]) / (2.0 * hy);
        //------------------------------------------------↑ Calculate gradients of u, k, omega, Dk, Dw ↑------------------------------------------------

        //------------------------------------------------↓ Calculate nut_star ↓------------------------------------------------
        std::vector<double> turbu_omega_tilde(ny);
        std::vector<double> nut_star(ny);
        for (int i = 0; i < ny; ++i) 
        {
            turbu_omega_tilde[i] = std::max(
                turbu_omega_star[i],
                std_kw_C_lim_ * dudy_discretized[i] / std_kw_beta_star_5_
            );
            nut_star[i] = k_star[i] / (turbu_omega_tilde[i] + tiny);
        }
        //------------------------------------------------↑ Calculate nut_star ↑------------------------------------------------

        //------------------------------------------------↓ Calculate analytical gradient of u ↓------------------------------------------------
        std::vector<double> dudy_star(ny);
        for (int i = 0; i < ny; ++i) 
        {
            dudy_star[i] = C_su[i] / (nu + nut_star[i]);
        }
        //------------------------------------------------↑ Calculate analytical gradient of u ↑------------------------------------------------
         
        //------------------------------------------------↓ Update linearized source terms Sc Sp ↓------------------------------------------------
        // 
        //-------------------------------------↓ For turbulent specific dissipation ↓-------------------------------------
        // calcualte part_cross_diffusion 
        std::vector<double> grad_prod(ny);
        std::vector<double> part_cross_diffusion(ny);
        for (int i = 0; i < ny; ++i) {
            grad_prod[i] = dkdy[i] * dwdy[i];
            double sigma_d = (grad_prod[i] > 0.0) ? sigma_d_value : 0.0;
            part_cross_diffusion[i] = sigma_d / (turbu_omega_star[i] + tiny) * grad_prod[i];
        }
        //-------------------------------------↑ For turbulent specific dissipation ↑-------------------------------------
        // 
        //------------------------------------------------↑ Update linearized source terms Sc Sp ↑------------------------------------------------

        //------------------------------------------------↓ Start solution using TDMA ↓------------------------------------------------
        // 
        //-------------------------------------↓ For velocity ↓-------------------------------------
        std::vector<double> a_u(ny, 0.0);
        std::vector<double> b_u(ny, 0.0);
        std::vector<double> c_u(ny, 0.0);
        std::vector<double> d_u(ny, 0.0);
        // inner node
        for (int i = 1; i < ny - 1; ++i) 
        {
            double nu_eff_i = nu + nut_star[i];
            double nu_eff_i_plus = nu + nut_star[i+1];
            double nu_eff_i_minus = nu + nut_star[i-1];
            double nu_eff_i_plus_half = 2.0 * nu_eff_i_plus * nu_eff_i / (nu_eff_i_plus + nu_eff_i);
            double nu_eff_i_minus_half = 2.0 * nu_eff_i_minus * nu_eff_i / (nu_eff_i_minus + nu_eff_i);
            a_u[i] = -nu_eff_i_minus_half;
            b_u[i] = (nu_eff_i_plus_half + nu_eff_i_minus_half);
            c_u[i] = -nu_eff_i_plus_half;
            d_u[i] = (utau * utau - tau_over_rho_outer) / height_sublayer * hy * hy;
        }
        // first node
        a_u[0] = 0.0;
        b_u[0] = 1.0;
        c_u[0] = 0.0;
        d_u[0] = u_p;
        // last node
        last = ny - 1;
        double nu_eff_last = nu + nut_star[last];
        double nu_eff_last_plus = nu + nut_p_outer; // B.C.
        double nu_eff_last_minus = nu + nut_star[last - 1];
        double nu_eff_last_plus_half = 2.0 * nu_eff_last_plus * nu_eff_last / (nu_eff_last_plus + nu_eff_last);
        double nu_eff_last_minus_half = 2.0 * nu_eff_last_minus * nu_eff_last / (nu_eff_last_minus + nu_eff_last);
        a_u[last] = -nu_eff_last_minus_half;
        b_u[last] = (nu_eff_last_plus_half + nu_eff_last_minus_half);
        c_u[last] = 0.0;
        d_u[last] = (utau * utau - tau_over_rho_outer) / height_sublayer * hy * hy + nu_eff_last_plus_half * u_p_outer;
        // solving
        std::vector<double> U_new = tdma(a_u, b_u, c_u, d_u);
        //-------------------------------------↑ For velocity ↑-------------------------------------
        
        //-------------------------------------↓ For turbulent kinetic energy ↓-------------------------------------
        std::vector<double> a_k(ny, 0.0);
        std::vector<double> b_k(ny, 0.0);
        std::vector<double> c_k(ny, 0.0);
        std::vector<double> d_k(ny, 0.0);
        // inner node
        for (int i = 1; i < ny - 1; ++i) {
            double Dk_i = diffusion_coefficient_k[i];
            double Dk_i_plus = diffusion_coefficient_k[i + 1];
            double Dk_i_minus = diffusion_coefficient_k[i - 1];
            double Dk_i_plus_half = 2.0 * Dk_i_plus * Dk_i / (Dk_i_plus + Dk_i);
            double Dk_i_minus_half = 2.0 * Dk_i_minus * Dk_i / (Dk_i_minus + Dk_i);
            a_k[i] = -1.0 * Dk_i_minus_half;
            b_k[i] = Dk_i_plus_half + Dk_i_minus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[i];
            c_k[i] = -1.0 * Dk_i_plus_half;
            d_k[i] = hy * hy * nut_star[i] * dudy_star[i] * dudy_star[i];
        }
        // first node, grad k = 0
        double Dk_first = diffusion_coefficient_k[0];
        double Dk_first_plus = diffusion_coefficient_k[0 + 1];
        double Dk_first_plus_half = 2.0 * Dk_first_plus * Dk_first / (Dk_first_plus + Dk_first);
        a_k[0] = 0.0;
        b_k[0] = Dk_first_plus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[0];
        c_k[0] = -1.0 * Dk_first_plus_half;
        d_k[0] = hy * hy * nut_star[0] * dudy_star[0] * dudy_star[0];
        // last node
        last = ny - 1;
        double Dk_last = diffusion_coefficient_k[last];
        double Dk_last_plus = nu + std_kw_sigma_star_ * k_p_outer / (w_p_outer + tiny);
        double Dk_last_minus = diffusion_coefficient_k[last - 1];
        double Dk_last_plus_half = 2.0 * Dk_last_plus * Dk_last / (Dk_last_plus + Dk_last);
        double Dk_last_minus_half = 2.0 * Dk_last_minus * Dk_last / (Dk_last_minus + Dk_last);
        a_k[last] = -1.0 * Dk_last_minus_half;
        b_k[last] = Dk_last_plus_half + Dk_last_minus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[last];
        c_k[last] = 0.0;
        d_k[last] = hy * hy * nut_star[last] * dudy_star[last] * dudy_star[last] + Dk_last_plus_half * k_p_outer;
        // solving
        std::vector<double> K_new = tdma(a_k, b_k, c_k, d_k);
        // avoid negative value, K_new = max(K_new, k_min)
        double k_min = 1e-10;
        for (int i = 0; i < ny; ++i) {
            K_new[i] = std::max(K_new[i], k_min);
        }
        //-------------------------------------↑ For turbulent kinetic energy ↑-------------------------------------

        //-------------------------------------↓ For turbulent specific dissipation ↓-------------------------------------
        std::vector<double> a_w(ny, 0.0);
        std::vector<double> b_w(ny, 0.0);
        std::vector<double> c_w(ny, 0.0);
        std::vector<double> d_w(ny, 0.0);
        // inner node i = 1 ... ny-2
        for (int i = 1; i < ny - 1; ++i) {
            double Dw_i = diffusion_coefficient_turbu_omega[i];
            double Dw_i_plus = diffusion_coefficient_turbu_omega[i + 1];
            double Dw_i_minus = diffusion_coefficient_turbu_omega[i - 1];
            double Dw_i_plus_half = 2.0 * Dw_i_plus * Dw_i / (Dw_i_plus + Dw_i);
            double Dw_i_minus_half = 2.0 * Dw_i_minus * Dw_i / (Dw_i_minus + Dw_i);
            a_w[i] = -1.0 * Dw_i_minus_half;
            b_w[i] = Dw_i_plus_half + Dw_i_minus_half + hy * hy * std_kw_beta_ * turbu_omega_star[i];
            c_w[i] = -1.0 * Dw_i_plus_half;
            double part_production = std_kw_alpha_ * turbu_omega_star[i] / (k_star[i] + tiny) * nut_star[i] * dudy_star[i] * dudy_star[i];
            d_w[i] = hy * hy * (part_production + part_cross_diffusion[i]);
        }
        // first node, w = Wp
        a_w[0] = 0.0;
        b_w[0] = 1.0;
        c_w[0] = 0.0;
        d_w[0] = turbu_omega_p;
        // last node
        last = ny - 1;
        double Dw_last = diffusion_coefficient_turbu_omega[last];
        double Dw_last_plus = nu + std_kw_sigma_ * k_p_outer / (w_p_outer + tiny);
        double Dw_last_minus = diffusion_coefficient_turbu_omega[last - 1];
        double Dw_last_plus_half = 2.0 * Dw_last_plus * Dw_last / (Dw_last_plus + Dw_last);
        double Dw_last_minus_half = 2.0 * Dw_last_minus * Dw_last / (Dw_last_minus + Dw_last);
        a_w[last] = -1.0 * Dw_last_minus_half;
        b_w[last] = Dw_last_plus_half + Dw_last_minus_half + hy * hy * std_kw_beta_ * turbu_omega_star[last];
        c_w[last] = 0.0;
        double part_production_last = std_kw_alpha_ * turbu_omega_star[last] / (k_star[last] + tiny) * nut_star[last] * dudy_star[last] * dudy_star[last];
        d_w[last] = hy * hy * (part_production_last + part_cross_diffusion[last]) + Dw_last_plus_half * w_p_outer;
        // solving
        std::vector<double> Turbu_omega_new = tdma(a_w, b_w, c_w, d_w);
        // avoid negative value, Turbu_omega_new = max(Turbu_omega_new, omega_min)
        double omega_min = 1e-10;
        for (int i = 0; i < ny; ++i) {
            Turbu_omega_new[i] = std::max(Turbu_omega_new[i], omega_min);
        }
        //-------------------------------------↑ For turbulent specific dissipation ↑-------------------------------------
        // 
        //------------------------------------------------↑ Start solution using TDMA ↑------------------------------------------------
        
        //------------------------------------------------↓ update phi_solved with under-relaxation ↓------------------------------------------------
        n_start = 0;
        double relax_u = 0.9;
        for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_u) * u_star[i] + relax_u * U_new[i];
        n_start += ny;
        double relax_k = 0.9;
        for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_k) * k_star[i] + relax_k * K_new[i];
        n_start += ny;
        double relax_w = 0.9;
        for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_w) * turbu_omega_star[i] + relax_w * Turbu_omega_new[i];
        //------------------------------------------------↑ update phi_solved with under-relaxation ↑------------------------------------------------
        
        //std::cout << "phi_solved = ";
        //for (const auto& v : phi_solved) std::cout << v << " ";
        //std::cout << std::endl;

        //------------------------------------------------↓ Check and update flow rate ↓------------------------------------------------
        flow_rate_current = accumulate(U_new.begin(), U_new.end(), 0.0) * hy;
        flow_rate_current += u_p_outer * hy * 0.5;

        // flow control
        double ratio = flow_rate_target / (flow_rate_current + 1e-12);
        double alpha = 0.03;
        double utau_new = utau * std::sqrt(ratio);
        utau = (1.0 - alpha) * utau + alpha * utau_new;
        /*double error = U_new[ny-1] - u_p_outer;
        utau -= 0.1 * error;*/
        //------------------------------------------------↑ Check and update flow rate ↑------------------------------------------------

        std::cout << "flow_rate_current = " << flow_rate_current
            << ", target = " << flow_rate_target << std::endl;
        std::cout << "updated utau = " << utau << std::endl;

        //------------------------------------------------↓ Calculate residue ↓------------------------------------------------
        differ = 0.0;
        for (size_t i = 0; i < phi_current.size(); ++i) {
            double diff = phi_current[i] - phi_solved[i];
            differ += diff * diff;
        }
        differ = std::sqrt(differ);
        //------------------------------------------------↑ Calculate residue ↑------------------------------------------------
        std::cout << "differ: " << differ << std::endl;

        //------------------------------------------------↓ Update ↓------------------------------------------------
        phi_current = phi_solved;
        num_iter_out += 1;
        //------------------------------------------------↑ Update ↑------------------------------------------------

        std::cout << "num_iter_out = " << num_iter_out << std::endl;
        std::cout << "------------" << std::endl;

    }
    std::cout << "******Converge******" << std::endl;

    std::cout << "******The results are: ******" << std::endl;

    n_start = 0;
    std::cout << "U = ";
    for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
    std::cout << std::endl;

    n_start += ny;
    std::cout << "K = ";
    for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
    std::cout << std::endl;

    n_start += ny;
    std::cout << "Turbu_omega = ";
    for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
    std::cout << std::endl;

    // ================== 提取解 ==================
    n_start = 0;

    std::vector<double> U(ny);
    for (int i = 0; i < ny; ++i) U[i] = phi_solved[n_start + i];

    n_start += ny;
    std::vector<double> K(ny);
    for (int i = 0; i < ny; ++i) K[i] = phi_solved[n_start + i];

    n_start += ny;
    std::vector<double> OMEGA(ny);
    for (int i = 0; i < ny; ++i) OMEGA[i] = phi_solved[n_start + i];

    // ================== 使用已有的 y ==================
    std::vector<double> Y = y;

    // ================== 计算湍流粘性系数 ==================
    double eps = 1e-12;
    std::vector<double> NUT(ny);
    for (int i = 0; i < ny; ++i) {
        NUT[i] = K[i] / (OMEGA[i] + eps);
    }

    // ================== 输出 Tecplot 文件 ==================
    std::string header_line = "ZONE T=\"SPH(1D)45 NF="
        + std::to_string(NF)
        + " ("
        + std::to_string(index)
        + ")\"";
    std::string filename = "pipe_kw_nf" + std::to_string(NF) + "_" + std::to_string(index) + ".dat";

    std::ofstream fout(filename);
    fout << "$VARIABLES = \"Y\", \"U\", \"K\", \"OMEGA\", \"NUT(k/omega)\"\n";
    fout << "$friction velocity = " << utau <<"\n";
    fout << "$current flow rate = " << flow_rate_current << "\n";
    fout << "$num_iter_out = " << num_iter_out << "\n";
    fout << header_line << "\n";

    fout << std::scientific << std::setprecision(8);

    for (int i = 0; i < ny; ++i) {
        fout << Y[i] << " "
            << U[i] << " "
            << K[i] << " "
            << OMEGA[i] << " "
            << NUT[i] << "\n";
    }

    fout.close();

    std::cout << "Tecplot 文件已生成: " << filename << std::endl;
}

// ================= Simpson 测试 =================
void test_simpson()
{
    double y_max = 1.0;
    int ny = 5;
    double hy = y_max / (double(ny) + 0.5); // distance from node U to P_outer is hy, hence with a 0.5
    double y_p = 0.5 * hy;
    Vec y(ny);  //computational nodes
    vector<double> U(ny);
    for (int i = 0; i < ny; ++i)
    {
        y[i] = y_p + i * hy;
        U[i] = y[i] * y[i];
    }
    double U_max = 1.0 * 1.0;

    // Simpson 积分
    double Q_simpson = compute_flow_rate_simpson_safe(U, hy);

    // 矩形法（你现在用的）
    double Q_rect = accumulate(U.begin(), U.end(), 0.0) * hy;

    // 精确解
    double Q_exact = 1.0 / 3.0 - 0.5 * hy * U_max;

    cout << "Exact integral     = " << Q_exact << endl;
    cout << "Simpson integral   = " << Q_simpson << endl;
    cout << "Rectangle integral = " << Q_rect << endl;

    cout << "Simpson error      = " << fabs(Q_simpson - Q_exact) << endl;
    cout << "Rectangle error    = " << fabs(Q_rect - Q_exact) << endl;
}

int main()
{
    //solve_1D_half_channel();
    double U_p = 3.674356e-01;
    double K_p = 2.145868e-03;
    double W_p = 4.179796e+01;
    double NUT_p = 5.166097e-05;
    double vel_grad_p = 1.021562e+01;
    double distance_to_wall = 0.03125;
    double target_flow_rate = 5.871126e-03;
    solve_1D_sublayer(U_p, 
        K_p, W_p, vel_grad_p, NUT_p, 
        distance_to_wall, target_flow_rate);
    //test_simpson();
    return 0;
}