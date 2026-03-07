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
    P_refinement_GetVelocityGradientInnerOnlyP::
        P_refinement_GetVelocityGradientInnerOnlyP(BaseInnerRelation& inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
        velocity_gradient_inner_only_P_(particles_->registerStateVariableData<Matd>("VelocityGradientInnerOnlyP")),
        turbu_B_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
        Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
        vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
        is_near_wall_P1_(this->particles_->template getVariableDataByName<int>("IsNearWallP1")) {}
    //=================================================================================================//
    void P_refinement_GetVelocityGradientInnerOnlyP::interaction(size_t index_i, Real dt)
    {
        velocity_gradient_inner_only_P_[index_i] = Matd::Zero();
        if (is_near_wall_P1_[index_i] == 1)
        {
            Vecd vel_i = vel_[index_i];
            const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
                velocity_gradient_inner_only_P_[index_i] += -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
            }
        }
    }
    //=================================================================================================//
    void P_refinement_GetVelocityGradientInnerOnlyP::update(size_t index_i, Real dt)
    {
        if (is_near_wall_P1_[index_i] == 1)
        {
            velocity_gradient_inner_only_P_[index_i] *= turbu_B_[index_i];
        }
    }
//=================================================================================================//
    P_refinement::
        P_refinement(SPHBody& sph_body)
        : LocalDynamics(sph_body),
        num_sub_node_(5),
        friction_velocity_from_sublayer_(particles_->registerStateVariableData<Real>("FrictionVelocityFromSublayer")),
        //
        is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
        y_p_(particles_->getVariableDataByName<Real>("Y_P")),
        vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
        turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
        turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
        rho_(particles_->getVariableDataByName<Real>("Density")),
        viscosity_(DynamicCast<Viscosity>(this, particles_->getBaseMaterial())), 
        mu_(viscosity_.ReferenceViscosity()),
        velocity_gradient_inner_only_P_(particles_->getVariableDataByName<Matd>("VelocityGradientInnerOnlyP")),
        turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
        wall_shear_stress_(particles_->getVariableDataByName<Real>("WallShearStress")),
        e_nearest_normal_(particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector")),
        fluid_particle_spacing_(sph_body.getSPHAdaptation().ReferenceSpacing()) 
    {
        particles_->addVariableToWrite<Real>("FrictionVelocityFromSublayer");
    }
    //=================================================================================================//
    void P_refinement::update(size_t index_i, Real dt)
    {
        friction_velocity_from_sublayer_[index_i] = 0.0;
        if (is_near_wall_P1_[index_i] == 1)
        {
            Vecd normal = e_nearest_normal_[index_i];

            Real nu = mu_ / rho_[index_i];
            Real u_outer = obtainTangentialComponent(vel_[index_i], normal);
            Real k_outer = turbu_k_[index_i];
            Real omega_outer = turbu_omega_[index_i];
            
            Vecd dudn_vector = velocity_gradient_inner_only_P_[index_i] * normal;
            Real dudn = obtainTangentialComponent(dudn_vector, normal);
            
            Real nut_outer = turbu_mu_[index_i] / rho_[index_i];
            Real distance_to_wall = y_p_[index_i];
            Real friction_vel_magnitude = std::sqrt(wall_shear_stress_[index_i] / rho_[index_i]);

            Real u_ps = u_outer + dudn * 0.5 * fluid_particle_spacing_;
            Real flow_rate_half = (u_outer + u_ps) * fluid_particle_spacing_ / 4.0;
            Real flow_rate_whole = u_outer * fluid_particle_spacing_;
            Real flow_rate_local = flow_rate_whole - flow_rate_half;
            
            friction_velocity_from_sublayer_[index_i] = solve_1D_sublayer(nu, u_outer, k_outer, omega_outer, std::abs(dudn),
                nut_outer, distance_to_wall, friction_vel_magnitude, std::abs(flow_rate_local));
        }
    }
    //=================================================================================================//
    double P_refinement::solve_1D_sublayer(double kinematic_viscosity, double u_p_outer, double k_p_outer, double w_p_outer,
        double vel_grad_p_outer, double nut_p_outer, double h_sublayer, double utau_outer, double Q_target)
    {
        //------------------------------------------------¡ý Input parameters ¡ý------------------------------------------------
        double utau_init = utau_outer;
        double height_sublayer = h_sublayer;
        double nu = kinematic_viscosity;

        double u_init = u_p_outer;
        double k_init = k_p_outer;
        double turbu_omega_init = w_p_outer;

        double convergence_criteria_outer = 1.0e-6;
        double tiny = 1.0e-6;
        double relax_u = 0.9;
        double relax_k = 0.9;
        double relax_w = 0.9;
        double alpha = 0.03;

        double flow_rate_target = Q_target;
        double utau = utau_init;
        //------------------------------------------------¡ü Input parameters ¡ü------------------------------------------------

        //------------------------------------------------¡ý Node arrangement, for sublayer ¡ý------------------------------------------------
        int ny = num_sub_node_;
        double hy = height_sublayer / (double(ny) + 0.5); // distance from node U to P_outer is hy, hence with a 0.5
        double y_p = 0.5 * hy;
        std::vector<double> y(ny);  //computational nodes
        for (int i = 0; i < ny; ++i)
        {
            y[i] = y_p + i * hy;
        }
        //------------------------------------------------¡ü Node arrangement, for sublayer ¡ü------------------------------------------------

        //------------------------------------------------¡ý Calculate P value ¡ý------------------------------------------------
        double yplus = utau * y_p / nu;
        double u_p = utau * yplus;
        double turbu_omega_p = 6.0 * nu / (std_kw_beta_i_ * y_p * y_p);
        //------------------------------------------------¡ü Calculate P value ¡ü------------------------------------------------

        //printf("hy=%f\n", hy);
        //printf("yp=%f\n", y_p);
        //printf("yplus=%f\n", yplus);
        //printf("u_p=%f\n", u_p);
        //std::cout << "ny= " << ny << std::endl;
        //std::cout << "y = ";
        //for (const auto& v : y) {
        //    std::cout << v << " ";
        //}
        //std::cout << std::endl;
        //std::cout << "press to continue" << std::endl;
        //std::cin.get();

        //------------------------------------------------¡ý Construct initial value ¡ý------------------------------------------------
        std::vector<double> u_init_value(ny, u_init);
        std::vector<double> k_init_value(ny, k_init);
        std::vector<double> turbu_omega_init_value(ny, turbu_omega_init);
        //------------------------------------------------¡ü Construct initial value ¡ü------------------------------------------------

        //------------------------------------------------¡ý Construct solution vector ¡ý------------------------------------------------
        std::vector<double> phi_current(0, 0.0);
        std::vector<double> phi_solved(3 * ny, 0.0);
        phi_current.insert(phi_current.end(), u_init_value.begin(), u_init_value.end());
        phi_current.insert(phi_current.end(), k_init_value.begin(), k_init_value.end());
        phi_current.insert(phi_current.end(), turbu_omega_init_value.begin(), turbu_omega_init_value.end());
        phi_solved = phi_current;
        //------------------------------------------------¡ü Construct solution vector ¡ü------------------------------------------------

        double differ = 1.0; // Should have a value 
        int num_iter_out = 0;
        int n_start = 0;
        int last = 0;
        int first_index = 0;
        double flow_rate_current = 0.0;
        while (differ > convergence_criteria_outer)
        {
            //------------------------------------------------¡ý Update the star values ¡ý------------------------------------------------
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
            //------------------------------------------------¡ü Update the star values ¡ü------------------------------------------------

            //------------------------------------------------¡ý Calculate Dk, Dw, C_su ¡ý------------------------------------------------
            std::vector<double> diffusion_coefficient_k(ny);
            std::vector<double> diffusion_coefficient_turbu_omega(ny);
            std::vector<double> C_su(ny);
            double tau_over_rho_outer = (nu + nut_p_outer) * vel_grad_p_outer;
            for (int i = 0; i < ny; ++i) {
                diffusion_coefficient_k[i] = nu + std_kw_sigma_star_ * k_star[i] / (turbu_omega_star[i] + tiny);
                diffusion_coefficient_turbu_omega[i] = nu + std_kw_sigma_ * k_star[i] / (turbu_omega_star[i] + tiny);
                C_su[i] = (utau * utau - tau_over_rho_outer) * (1.0 - y[i] / height_sublayer) + tau_over_rho_outer;
            }
            //------------------------------------------------¡ü Calculate Dk, Dw  ¡ü------------------------------------------------

            //------------------------------------------------¡ý Calculate gradients of u, k, omega, Dk, Dw ¡ý------------------------------------------------
            std::vector<double> dudy_discretized(ny, 0.0);
            std::vector<double> dkdy(ny, 0.0);
            std::vector<double> dwdy(ny, 0.0);
            // central diff£¨from i=1 to i=ny-2 £©
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
            //------------------------------------------------¡ü Calculate gradients of u, k, omega, Dk, Dw ¡ü------------------------------------------------

            //------------------------------------------------¡ý Calculate nut_star ¡ý------------------------------------------------
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
            //------------------------------------------------¡ü Calculate nut_star ¡ü------------------------------------------------

            //------------------------------------------------¡ý Calculate analytical gradient of u ¡ý------------------------------------------------
            std::vector<double> dudy_star(ny);
            for (int i = 0; i < ny; ++i)
            {
                dudy_star[i] = C_su[i] / (nu + nut_star[i]);
            }
            //------------------------------------------------¡ü Calculate analytical gradient of u ¡ü------------------------------------------------

            //------------------------------------------------¡ý Update linearized source terms Sc Sp ¡ý------------------------------------------------
            // 
            //-------------------------------------¡ý For turbulent specific dissipation ¡ý-------------------------------------
            // calcualte part_cross_diffusion 
            std::vector<double> grad_prod(ny);
            std::vector<double> part_cross_diffusion(ny);
            for (int i = 0; i < ny; ++i) {
                grad_prod[i] = dkdy[i] * dwdy[i];
                double sigma_d = (grad_prod[i] > 0.0) ? std_kw_sigma_do_ : 0.0;
                part_cross_diffusion[i] = sigma_d / (turbu_omega_star[i] + tiny) * grad_prod[i];
            }
            //-------------------------------------¡ü For turbulent specific dissipation ¡ü-------------------------------------
            // 
            //------------------------------------------------¡ü Update linearized source terms Sc Sp ¡ü------------------------------------------------

            //------------------------------------------------¡ý Start solution using TDMA ¡ý------------------------------------------------
            // 
            //-------------------------------------¡ý For velocity ¡ý-------------------------------------
            std::vector<double> a_u(ny, 0.0);
            std::vector<double> b_u(ny, 0.0);
            std::vector<double> c_u(ny, 0.0);
            std::vector<double> d_u(ny, 0.0);
            // inner node
            for (int i = 1; i < ny - 1; ++i)
            {
                double nu_eff_i = nu + nut_star[i];
                double nu_eff_i_plus = nu + nut_star[i + 1];
                double nu_eff_i_minus = nu + nut_star[i - 1];
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
            //-------------------------------------¡ü For velocity ¡ü-------------------------------------

            //-------------------------------------¡ý For turbulent kinetic energy ¡ý-------------------------------------
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
            //-------------------------------------¡ü For turbulent kinetic energy ¡ü-------------------------------------

            //-------------------------------------¡ý For turbulent specific dissipation ¡ý-------------------------------------
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
            //-------------------------------------¡ü For turbulent specific dissipation ¡ü-------------------------------------
            // 
            //------------------------------------------------¡ü Start solution using TDMA ¡ü------------------------------------------------

            //------------------------------------------------¡ý update phi_solved with under-relaxation ¡ý------------------------------------------------
            n_start = 0;
            for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_u) * u_star[i] + relax_u * U_new[i];
            n_start += ny;
            for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_k) * k_star[i] + relax_k * K_new[i];
            n_start += ny;
            for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_w) * turbu_omega_star[i] + relax_w * Turbu_omega_new[i];
            //------------------------------------------------¡ü update phi_solved with under-relaxation ¡ü------------------------------------------------

            //std::cout << "phi_solved = ";
            //for (const auto& v : phi_solved) std::cout << v << " ";
            //std::cout << std::endl;

            //------------------------------------------------¡ý Check and update flow rate ¡ý------------------------------------------------
            flow_rate_current = accumulate(U_new.begin(), U_new.end(), 0.0) * hy;
            flow_rate_current += u_p_outer * hy * 0.5;

            // flow control
            double ratio = flow_rate_target / (flow_rate_current + 1e-12);
            double utau_new = utau * std::sqrt(ratio);
            utau = (1.0 - alpha) * utau + alpha * utau_new;
            /*double error = U_new[ny-1] - u_p_outer;
            utau -= 0.1 * error;*/
            //------------------------------------------------¡ü Check and update flow rate ¡ü------------------------------------------------

            //std::cout << "flow_rate_current = " << flow_rate_current
            //    << ", target = " << flow_rate_target << std::endl;
            //std::cout << "updated utau = " << utau << std::endl;

            //------------------------------------------------¡ý Calculate residue ¡ý------------------------------------------------
            differ = 0.0;
            for (size_t i = 0; i < phi_current.size(); ++i) {
                double diff = phi_current[i] - phi_solved[i];
                differ += diff * diff;
            }
            differ = std::sqrt(differ);
            //------------------------------------------------¡ü Calculate residue ¡ü------------------------------------------------
            //std::cout << "differ: " << differ << std::endl;

            //------------------------------------------------¡ý Update ¡ý------------------------------------------------
            phi_current = phi_solved;
            num_iter_out += 1;
            //------------------------------------------------¡ü Update ¡ü------------------------------------------------

            //std::cout << "num_iter_out = " << num_iter_out << std::endl;
            //std::cout << "------------" << std::endl;

        }
        //std::cout << "******Converge******" << std::endl;

        //std::cout << "******The results are: ******" << std::endl;

        //n_start = 0;
        //std::cout << "U = ";
        //for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
        //std::cout << std::endl;

        //n_start += ny;
        //std::cout << "K = ";
        //for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
        //std::cout << std::endl;

        //n_start += ny;
        //std::cout << "Turbu_omega = ";
        //for (int i = 0; i < ny; ++i) std::cout << phi_solved[n_start + i] << " ";
        //std::cout << std::endl;

        // ================== Extract Solution ==================
        //n_start = 0;

        //std::vector<double> U(ny);
        //for (int i = 0; i < ny; ++i) U[i] = phi_solved[n_start + i];

        //n_start += ny;
        //std::vector<double> K(ny);
        //for (int i = 0; i < ny; ++i) K[i] = phi_solved[n_start + i];

        //n_start += ny;
        //std::vector<double> OMEGA(ny);
        //for (int i = 0; i < ny; ++i) OMEGA[i] = phi_solved[n_start + i];

        //// ================== use existing y ==================
        //std::vector<double> Y = y;

        //// ================== calculate nut ==================
        //double eps = 1e-12;
        //std::vector<double> NUT(ny);
        //for (int i = 0; i < ny; ++i) {
        //    NUT[i] = K[i] / (OMEGA[i] + eps);
        //}
        return utau;
    }
    //=================================================================================================//
    // ================= TDMA =================
    // Solve a tridiagonal system: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
    // a[0] must be 0, c[n-1] will be ignored
    std::vector<double> P_refinement::tdma(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) 
    {
        int n = d.size();
        if (a.size() != n || b.size() != n || c.size() != n) {
            throw std::invalid_argument("TDMA: vector size mismatch!");
        }

        std::vector<double> cp(n, 0.0); // modified upper diagonal
        std::vector<double> dp(n, 0.0); // modified right-hand side
        std::vector<double> x(n, 0.0);  // solution vector

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
//=================================================================================================//
} // namespace udf
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//