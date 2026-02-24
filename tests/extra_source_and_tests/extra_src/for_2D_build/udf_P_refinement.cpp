//#pragma once
#include "udf_P_refinement.hpp"
namespace SPH
{
//=================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
namespace udf
{
//=================================================================================================//
    P_refinement::
        P_refinement(SPHBody& sph_body)
        : LocalDynamics(sph_body),
        num_sub_node_(5),
        //
        is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
        y_p_(particles_->getVariableDataByName<Real>("Y_P")),
        vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
        turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
        turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
        rho_(particles_->getVariableDataByName<Real>("Density")),
        viscosity_(DynamicCast<Viscosity>(this, particles_->getBaseMaterial())), 
        mu_(viscosity_.ReferenceViscosity()) {}
    
    //=================================================================================================//
    void P_refinement::update(size_t index_i, Real dt)
    {
        if (is_near_wall_P1_[index_i] == 1)
        {
            Real nu = mu_ / rho_[index_i];
            Real u_init = vel_[index_i][xAxis];  //** Temporary treatment *
            Real k_init = turbu_k_[index_i];
            Real omega_init = turbu_omega_[index_i];
            Real distance_to_wall = y_p_[index_i];
            Real hy = distance_to_wall / (Real(num_sub_node_) + 0.5);
            Real y_p_sub = 0.5 * hy;
            std::vector<Real> y(num_sub_node_);
            for (int i = 0; i < num_sub_node_; ++i) 
            {
                y[i] = y_p_sub + i * hy;
            }
            std::vector<Real> u_init_value(num_sub_node_, u_init);
            std::vector<Real> k_init_value(num_sub_node_, k_init);
            std::vector<Real> turbu_omega_init_value(num_sub_node_, omega_init);
            
            std::vector<double> phi_current(0, 0.0);
            std::vector<double> phi_solved(3 * num_sub_node_, 0.0);
            int num_iter_out = 0;
            phi_current.insert(phi_current.end(), u_init_value.begin(), u_init_value.end());
            phi_current.insert(phi_current.end(), k_init_value.begin(), k_init_value.end());
            phi_current.insert(phi_current.end(), turbu_omega_init_value.begin(), turbu_omega_init_value.end());
            phi_solved = phi_current;

            Real convergence_criteria_outer = 1.0e-6;
            double tiny = 1.0e-6;
            double utau_init = 6.37309e-02;
            double utau = utau_init;
            double yplus_sub = utau * y_p_sub / nu;
            double u_p = utau * yplus_sub;
            double turbu_omega_p = 6.0 * nu / (std_kw_beta_i_ * y_p_sub * y_p_sub);

            double differ = 1.0;
            int n_start = 0;
            while (differ > convergence_criteria_outer)
            {
                // ---update the star values---
                n_start = 0;
                std::vector<double> u_star(
                    phi_solved.begin() + n_start,
                    phi_solved.begin() + n_start + num_sub_node_
                );

                n_start = num_sub_node_;
                std::vector<double> k_star(
                    phi_solved.begin() + n_start,
                    phi_solved.begin() + n_start + num_sub_node_
                );

                n_start = 2 * num_sub_node_;
                std::vector<double> turbu_omega_star(
                    phi_solved.begin() + n_start,
                    phi_solved.begin() + n_start + num_sub_node_
                );

                // initialisation 
                std::vector<double> dkdy(num_sub_node_, 0.0);
                std::vector<double> dwdy(num_sub_node_, 0.0);
                // backward diff（from i=1 ）
                for (int i = 1; i < num_sub_node_; ++i) {
                    dkdy[i] = (k_star[i] - k_star[i - 1]) / hy;
                    dwdy[i] = (turbu_omega_star[i] - turbu_omega_star[i - 1]) / hy;
                }
                // B.C.
                dkdy[0] = 0.0;
                dwdy[0] = 0.0;

                std::vector<double> dudy_discretized(num_sub_node_, 0.0);
                for (int i = 1; i < num_sub_node_; ++i) {
                    dudy_discretized[i] = (u_star[i] - u_star[i - 1]) / hy;
                }
                dudy_discretized[0] = 0.0;

                std::vector<double> turbu_omega_tilde(num_sub_node_);
                std::vector<double> nut_star(num_sub_node_);

                for (int i = 0; i < num_sub_node_; ++i) {
                    turbu_omega_tilde[i] = std::max(
                        turbu_omega_star[i],
                        std_kw_C_lim_ * dudy_discretized[i] / std_kw_beta_star_5_
                    );
                    nut_star[i] = k_star[i] / (turbu_omega_tilde[i] + tiny);
                }

                std::vector<double> dudy_star(num_sub_node_);
                for (int i = 0; i < num_sub_node_; ++i) {
                    dudy_star[i] = utau * utau * (1.0 - y[i] / delta) / (nu + nut_star[i]);
                }

                // ---update the linearized source term variables---
                // For velocity
                std::vector<double> Sc_u(num_sub_node_);
                std::vector<double> Sp_u(num_sub_node_, 0.0);

                for (int i = 0; i < num_sub_node_; ++i) {
                    Sc_u[i] = -1.0 / (nu + nut_star[i]);
                    // Sp_u already is 0
                }

                // diffusion coefficient for k
                std::vector<double> diffusion_coefficient_k(num_sub_node_);
                for (int i = 0; i < num_sub_node_; ++i) {
                    diffusion_coefficient_k[i] = nu + std_kw_sigma_star_ * k_star[i] / (turbu_omega_star[i] + tiny);
                }

                // 后向差分
                std::vector<double> dDkdy(num_sub_node_, 0.0);
                for (int i = 1; i < num_sub_node_; ++i) {
                    dDkdy[i] = (diffusion_coefficient_k[i] - diffusion_coefficient_k[i - 1]) / hy;
                }

                // 边界
                dDkdy[0] = 0.0;

                // part_extra_viscous_term_k = dDkdy * dkdy
                std::vector<double> part_extra_viscous_term_k(num_sub_node_);
                std::vector<double> Sc_k(num_sub_node_);
                std::vector<double> Sp_k(num_sub_node_);

                for (int i = 0; i < num_sub_node_; ++i) {
                    part_extra_viscous_term_k[i] = dDkdy[i] * dkdy[i];
                    Sc_k[i] = (nut_star[i] * dudy_star[i] * dudy_star[i] + part_extra_viscous_term_k[i])
                        / diffusion_coefficient_k[i];
                    Sp_k[i] = std_kw_beta_star_ * turbu_omega_star[i] / diffusion_coefficient_k[i];
                }

                // diffusion coefficient for turbu_omega
                std::vector<double> diffusion_coefficient_turbu_omega(num_sub_node_);
                for (int i = 0; i < num_sub_node_; ++i) {
                    diffusion_coefficient_turbu_omega[i] = nu + std_kw_sigma_ * k_star[i] / (turbu_omega_star[i] + tiny);
                }

                // 后向差分
                std::vector<double> dDwdy(num_sub_node_, 0.0);
                for (int i = 1; i < num_sub_node_; ++i) {
                    dDwdy[i] = (diffusion_coefficient_turbu_omega[i] - diffusion_coefficient_turbu_omega[i - 1]) / hy;
                }
                dDwdy[0] = 0.0;

                // part_extra_viscous_term_w = dDwdy * dwdy
                std::vector<double> part_extra_viscous_term_w(num_sub_node_);
                for (int i = 0; i < num_sub_node_; ++i) {
                    part_extra_viscous_term_w[i] = dDwdy[i] * dwdy[i];
                }

                // part_production = alpha*omega/k * nut * (dudy)^2
                std::vector<double> part_production(num_sub_node_);
                for (int i = 0; i < num_sub_node_; ++i) {
                    part_production[i] = std_kw_alpha_ * turbu_omega_star[i] / (k_star[i] + tiny)
                        * nut_star[i] * dudy_star[i] * dudy_star[i];
                }

                // grad_prod = dkdy * dwdy
                std::vector<double> grad_prod(num_sub_node_);
                std::vector<double> part_cross_diffusion(num_sub_node_);
                for (int i = 0; i < num_sub_node_; ++i) {
                    grad_prod[i] = dkdy[i] * dwdy[i];
                    double sigma_d = (grad_prod[i] > 0.0) ? std_kw_sigma_do_ : 0.0;
                    part_cross_diffusion[i] = sigma_d / (turbu_omega_star[i] + tiny) * grad_prod[i];
                }

                // Sc and Sp
                std::vector<double> Sc_turbu_omega(num_sub_node_);
                std::vector<double> Sp_turbu_omega(num_sub_node_);
                for (int i = 0; i < num_sub_node_; ++i) {
                    Sc_turbu_omega[i] = (part_production[i] + part_cross_diffusion[i] + part_extra_viscous_term_w[i])
                        / diffusion_coefficient_turbu_omega[i];
                    Sp_turbu_omega[i] = std_kw_beta_ * turbu_omega_star[i] / diffusion_coefficient_turbu_omega[i];
                }

                // 初始化
                std::vector<double> a_u(num_sub_node_, 0.0);
                std::vector<double> b_u(num_sub_node_, 0.0);
                std::vector<double> c_u(num_sub_node_, 0.0);
                std::vector<double> d_u(num_sub_node_, 0.0);

                // 循环赋值
                for (int i = 1; i < num_sub_node_; ++i) {
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
                std::vector<double> a_k(num_sub_node_, 0.0);
                std::vector<double> b_k(num_sub_node_, 0.0);
                std::vector<double> c_k(num_sub_node_, 0.0);
                std::vector<double> d_k(num_sub_node_, 0.0);

                // 内部节点 i = 1 ... num_sub_node_-2
                for (int i = 1; i < num_sub_node_ - 1; ++i) {
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
                int last = num_sub_node_ - 1;
                a_k[last] = 1.0;
                b_k[last] = -1.0 - Sp_k[last] * hy * hy;
                c_k[last] = 0.0;
                d_k[last] = -1.0 * hy * hy * Sc_k[last];

                // TDMA 求解
                std::vector<double> K_new = tdma(a_k, b_k, c_k, d_k);

                // 初始化
                std::vector<double> a_w(num_sub_node_, 0.0);
                std::vector<double> b_w(num_sub_node_, 0.0);
                std::vector<double> c_w(num_sub_node_, 0.0);
                std::vector<double> d_w(num_sub_node_, 0.0);

                // 内部节点 i = 1 ... num_sub_node_-2
                for (int i = 1; i < num_sub_node_ - 1; ++i) {
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
                last = num_sub_node_ - 1;
                a_w[last] = 1.0;
                b_w[last] = -1.0 - Sp_turbu_omega[last] * hy * hy;
                c_w[last] = 0.0;
                d_w[last] = -1.0 * hy * hy * Sc_turbu_omega[last];

                // TDMA 求解
                std::vector<double> Turbu_omega_new = tdma(a_w, b_w, c_w, d_w);

                // 更新 phi_solved
                n_start = 0;
                for (int i = 0; i < num_sub_node_; ++i) phi_solved[n_start + i] = U_new[i];

                n_start += num_sub_node_;
                for (int i = 0; i < num_sub_node_; ++i) phi_solved[n_start + i] = K_new[i];

                n_start += num_sub_node_;
                for (int i = 0; i < num_sub_node_; ++i) phi_solved[n_start + i] = Turbu_omega_new[i];

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

        }
    }
//=================================================================================================//
} // namespace udf
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//