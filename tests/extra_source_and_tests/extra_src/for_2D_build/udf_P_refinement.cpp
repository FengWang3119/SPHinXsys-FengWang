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
    P_refinement_GetVelocityGradient<Inner<>>::
        P_refinement_GetVelocityGradient(BaseInnerRelation& inner_relation)
        : P_refinement_GetVelocityGradient<DataDelegateInner>(inner_relation),
        turbu_B_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
        B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
    {
        particles_->addVariableToWrite<Matd>("VelocityGradientInnerOnlyP");
    }
    //=================================================================================================//
    void P_refinement_GetVelocityGradient<Inner<>>::interaction(size_t index_i, Real dt)
    {
        velocity_gradient_only_P_[index_i] = Matd::Zero();
        k_gradient_only_P_[index_i] = Vecd::Zero();
        omega_gradient_only_P_[index_i] = Vecd::Zero();
        if (is_near_wall_P1_[index_i] == 1)
        {
            Vecd vel_i = vel_[index_i];
            Real k_i = turbu_k_[index_i];
            Real w_i = turbu_omega_[index_i];
            const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
                velocity_gradient_only_P_[index_i] += -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
                k_gradient_only_P_[index_i] += -1.0 * (k_i - turbu_k_[index_j]) * nablaW_ijV_j;
                omega_gradient_only_P_[index_i] += -1.0 * (w_i - turbu_omega_[index_j]) * nablaW_ijV_j;
            }
        }
    }
    //=================================================================================================//
    void P_refinement_GetVelocityGradient<Inner<>>::update(size_t index_i, Real dt)
    {
        if (is_near_wall_P1_[index_i] == 1)
        {
            //velocity_gradient_only_P_[index_i] *= turbu_B_[index_i];
            velocity_gradient_only_P_[index_i] *= B_[index_i];
        }
    }
    //=================================================================================================//
    P_refinement_GetVelocityGradient<Contact<Wall>>::P_refinement_GetVelocityGradient(BaseContactRelation& contact_relation)
        : InteractionWithWall<P_refinement_GetVelocityGradient>(contact_relation) {}
    //=================================================================================================//
    void P_refinement_GetVelocityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
    {
        if (is_near_wall_P1_[index_i] == 1)
        {
            Matd vel_grad = Matd::Zero();
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Vecd* vel_ave_k = wall_vel_ave_[k];
                Real* Vol_k = wall_Vol_[k];
                Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
                    //vel_grad += -1.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
                    vel_grad += -1.0 * 2.0 * (vel_[index_i] - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
                }
            }
            velocity_gradient_only_P_[index_i] += vel_grad;
        }
    }
//=================================================================================================//
    P_refinement::
        P_refinement(SPHBody& sph_body)
        : LocalDynamics(sph_body),
        num_sub_node_(5), // ** Needs tobe modified in sublayer function, as well *
        friction_velocity_from_sublayer_(particles_->registerStateVariableData<Real>("FrictionVelocityFromSublayer")),
        target_flow_rate_in_sublayer_(particles_->registerStateVariableData<Real>("TargetFlowRateInSublayer")),
        vel_ps_magnitude_(particles_->registerStateVariableData<Real>("VelPS")),
        dudn_for_local_flow_rate_(particles_->registerStateVariableData<Real>("dudnForLocalFlowRate")),
        utau_node_(particles_->registerStateVariableData<Real>("utauNode")),
        node_value_(particles_->registerStateVariableData<Vec6d>("NodeValue")),
        dUdn_P_sublayer_magnitude_(particles_->registerStateVariableData<Real>("dUdnFromSublayerMagnitude")),
        dUdn_P_sublayer_(particles_->registerStateVariableData<Matd>("dUdnFromSublayer")),
        vel_nodeO_(particles_->registerStateVariableData<Real>("VelNodeO")),
        vel_nodeUM_(particles_->registerStateVariableData<Real>("VelNodeUM")),
        //
        is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
        y_p_(particles_->getVariableDataByName<Real>("Y_P")),
        vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
        turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
        turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
        rho_(particles_->getVariableDataByName<Real>("Density")),
        viscosity_(DynamicCast<Viscosity>(this, particles_->getBaseMaterial())), 
        mu_(viscosity_.ReferenceViscosity()),
        velocity_gradient_only_P_(particles_->getVariableDataByName<Matd>("VelocityGradientInnerOnlyP")),
        turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
        distance_to_dummy_interface_(particles_->getVariableDataByName<Real>("DistanceToDummyInterface")),
        wall_shear_stress_(particles_->getVariableDataByName<Real>("WallShearStress")),
        e_nearest_normal_(particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector")),
        fluid_particle_spacing_(sph_body.getSPHAdaptation().ReferenceSpacing()),
        physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")),
        velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
        k_gradient_only_P_(particles_->getVariableDataByName<Vecd>("TurbulentKineticEnergyGradientOnlyP")),
        omega_gradient_only_P_(particles_->getVariableDataByName<Vecd>("TurbulentSpecificDissipationGradientOnlyP"))
    {
        particles_->addVariableToWrite<Real>("FrictionVelocityFromSublayer");
        particles_->addVariableToWrite<Real>("TargetFlowRateInSublayer");
        particles_->addVariableToWrite<Real>("VelPS");
        particles_->addVariableToWrite<Real>("dudnForLocalFlowRate");
        particles_->addVariableToWrite<Real>("utauNode");
        particles_->addVariableToWrite<Real>("DistanceToDummyInterface");
        particles_->addVariableToWrite<Vec6d>("NodeValue");
        particles_->addVariableToWrite<Real>("dUdnFromSublayerMagnitude");
        particles_->addVariableToWrite<Matd>("dUdnFromSublayer");
        particles_->addVariableToWrite<Real>("VelNodeO");
        particles_->addVariableToWrite<Real>("VelNodeUM");
    }
    //=================================================================================================//
    void P_refinement::update(size_t index_i, Real dt)
    {
        Real sum_node_vel_difference = 0.0;
        Real vel_nodeO_i_prior = 0.0;
        if (is_near_wall_P1_[index_i] == 1)
        {
            Vec6d node_value_i_prior = node_value_[index_i];
            for (int j = 0; j < num_sub_node_; ++j)
            {
                Real vel_difference = node_value_i_prior[j + 1] - node_value_i_prior[j];
                sum_node_vel_difference += vel_difference;
            }
            vel_nodeO_i_prior = vel_nodeO_[index_i];
        }

        friction_velocity_from_sublayer_[index_i] = 0.0;
        target_flow_rate_in_sublayer_[index_i] = 0.0;
        vel_ps_magnitude_[index_i] = 0.0;
        dudn_for_local_flow_rate_[index_i] = 0.0;
        utau_node_[index_i] = 0.0;
        node_value_[index_i] = Vec6d::Zero();
        dUdn_P_sublayer_magnitude_[index_i] = 0.0;
        dUdn_P_sublayer_[index_i] = Matd::Zero();
        vel_nodeO_[index_i] = 0.0;
        vel_nodeUM_[index_i] = 0.0;
        double U_nodeO = 0.0;
        double U_nodeUM = 0.0;
        if (is_near_wall_P1_[index_i] == 1)
        {
            //** Define outside values, 3 grad values, 1 local flowrate, 5 initial values *
            Real dudn_outer = 0.0;
            Real dkdn_outer = 0.0;
            Real dwdn_outer = 0.0;
            Real flow_rate_local = 0.0;
            Real u_outer = 0.0;
            Real k_outer = 0.0;
            Real omega_outer = 0.0;
            Real nut_outer = 0.0;
            Real friction_vel_magnitude_outer = 0.0;

            Vecd normal = e_nearest_normal_[index_i];
            Real nu = mu_ / rho_[index_i];
            
            u_outer = obtainTangentialComponent(vel_[index_i], normal);
            k_outer = turbu_k_[index_i];
            omega_outer = turbu_omega_[index_i];
            nut_outer = turbu_mu_[index_i] / rho_[index_i];
            friction_vel_magnitude_outer = std::sqrt(wall_shear_stress_[index_i] / rho_[index_i]);

            //** Two ways to determine distance to wall *
            Real distance_to_wall = y_p_[index_i];
            //Real distance_to_wall = distance_to_dummy_interface_[index_i];
            
            Vecd velocity_gradient_only_P_normal = velocity_gradient_only_P_[index_i] * normal;
            Real dudn_from_SPH = obtainTangentialComponent(velocity_gradient_only_P_normal, normal);
            //** Two ways to determine dudn_outer *
            //** If use full from SPH *
            dudn_outer = dudn_from_SPH;
            //** If use weighting combination *
            //Real weight_SPH = fluid_particle_spacing_;
            //Real sum_weight_sublayer = distance_to_wall; //** Assume uniform division *
            //dudn_outer = (dudn_from_SPH * weight_SPH + sum_node_vel_difference) / (weight_SPH + sum_weight_sublayer); 

            dkdn_outer = k_gradient_only_P_[index_i].dot(normal); //** Currently, only inner contribution is considered *
            dwdn_outer = omega_gradient_only_P_[index_i].dot(normal);
            //------------------------------------------------¡ý For test 1D analytical ¡ý------------------------------------------------
            //
            //-------------------------¡ý If input fix analytical value ¡ý-------------------------
            if(0)
            {
                std::cout << "Fixed input value test starts." << std::endl;
                //** Define outside values *
                flow_rate_local = 3.781607e-3;
                dudn_outer = 1.127180e+1;
                k_outer = 1.118813e-3;
                omega_outer = 5.998469e+1;
                nut_outer = 1.885024e-5;
                u_outer = 3.000607e-1;
                friction_vel_magnitude_outer = 6.37309e-02;
                dkdn_outer = 1.391004e-01;
                dwdn_outer = -4.190042e+03;
                nu = 3.5e-4;
                U_nodeO = 0.0;
                U_nodeUM = 0.0;
                //** Remember to activate output function inside *
                node_value_[index_i] = solve_1D_sublayer(nu, u_outer, k_outer, omega_outer, std::abs(dudn_outer),
                    nut_outer, distance_to_wall, friction_vel_magnitude_outer, std::abs(flow_rate_local), dkdn_outer, dwdn_outer, U_nodeO, U_nodeUM);
                std::cout << "Fixed input value test ends, stop here." << std::endl;
                std::cin.get();
            }
            //-------------------------¡ü If input fix analytical value ¡ü-------------------------

            //-------------------------¡ý If dynamic test ¡ý-------------------------
            if (0)
            {
                std::cout << "Dynamic test starts." << std::endl;
                //------¡ý Mimic SPH average value ¡ý------
                Real analytical_k_S = 8.260674e-03;
                Real analytical_k_P = 1.118813e-3;
                Real analytical_k_grad_P_inner = (analytical_k_S - analytical_k_P) / fluid_particle_spacing_;

                Real analytical_w_S = 1.242432e+01;
                Real analytical_w_P = 5.998469e+1;
                Real analytical_w_grad_P_inner = (analytical_w_S - analytical_w_P) / fluid_particle_spacing_;

                Real analytical_vel_S = 6.487377e-01;
                Real analytical_vel_P = 3.000607e-1;
                Real analytical_vel_grad_P_inner = (analytical_vel_S - analytical_vel_P) / fluid_particle_spacing_;

                analytical_vel_grad_P_inner = 1.127180e+1; //** Overwrite *

                Real analytical_flow_rate_whole_PS_to_Wall = 1.433123e-2;
                //------¡ü Mimic SPH average value ¡ü------
                // 
                //------¡ý Start testing ¡ý------
                Real flow_rate_local_prior = 0.0;
                Real residue = 1.0e3;
                vel_nodeO_i_prior = analytical_vel_P;
                
                //** Fix 5+1 initial values *
                k_outer = 1.118813e-3;
                omega_outer = 5.998469e+1;
                nut_outer = 1.885024e-5;
                u_outer = 3.000607e-1;
                nu = 3.5e-4;
                friction_vel_magnitude_outer = 6.37309e-02;

                //** Transfer 3 gradient values *
                dudn_outer = analytical_vel_grad_P_inner;
                dkdn_outer = analytical_k_grad_P_inner;
                dwdn_outer = analytical_w_grad_P_inner;

                while (residue > 1.0e-6)
                {
                    flow_rate_local = get_loacal_flow_rate(analytical_flow_rate_whole_PS_to_Wall, dudn_outer, vel_nodeO_i_prior, fluid_particle_spacing_);

                    U_nodeO = 0.0;
                    U_nodeUM = 0.0;
                    node_value_[index_i] = solve_1D_sublayer(nu, u_outer, k_outer, omega_outer, std::abs(dudn_outer),
                        nut_outer, distance_to_wall, friction_vel_magnitude_outer, std::abs(flow_rate_local), dkdn_outer, dwdn_outer, U_nodeO, U_nodeUM);

                    residue = std::abs(flow_rate_local - flow_rate_local_prior);
                    std::cout << "residue =" << residue << std::endl;
                    flow_rate_local_prior = flow_rate_local;
                    Real relax_factor = 0.3;
                    vel_nodeO_i_prior = (1.0 - relax_factor) * vel_nodeO_i_prior + relax_factor * U_nodeO;
                }

                std::cout << "flow rate local converge!." << std::endl;
                std::cout << "node_value_[index_i][0]=" << node_value_[index_i][0] << std::endl;
                std::cout << "node_value_[index_i][1]=" << node_value_[index_i][1] << std::endl;
                std::cout << "node_value_[index_i][2]=" << node_value_[index_i][2] << std::endl;
                std::cout << "node_value_[index_i][3]=" << node_value_[index_i][3] << std::endl;
                std::cout << "node_value_[index_i][4]=" << node_value_[index_i][4] << std::endl;
                std::cout << "node_value_[index_i][5]=" << node_value_[index_i][5] << std::endl;
                std::cout << "U_nodeO=" << U_nodeO << std::endl;
                std::cout << "U_nodeUM=" << U_nodeUM << std::endl;

                writeTecplotFromVec6d(
                    node_value_[index_i],
                    U_nodeO,
                    U_nodeUM,
                    distance_to_wall,
                    num_sub_node_,
                    40,
                    80
                );

                std::cout << "Dynamic test ends, stop here." << std::endl;
                std::cin.get();
                //------¡ü Start testing ¡ü------
            }
            //-------------------------¡ü If dynamic test ¡ü-------------------------
            // 
            //------------------------------------------------¡ü For test 1D analytical ¡ü------------------------------------------------
            
            Real flow_rate_local_prior = 0.0;
            Real residue = 1.0e3;
            while (residue > 1.0e-6)
            {
                flow_rate_local = get_loacal_flow_rate(u_outer * fluid_particle_spacing_, dudn_outer, vel_nodeO_i_prior, fluid_particle_spacing_); //** This is for better testing *
                U_nodeO = 0.0;
                U_nodeUM = 0.0;
                node_value_[index_i] = solve_1D_sublayer(nu, u_outer, k_outer, omega_outer, std::abs(dudn_outer),
                    nut_outer, distance_to_wall, friction_vel_magnitude_outer, std::abs(flow_rate_local), dkdn_outer, dwdn_outer, U_nodeO, U_nodeUM);

                residue = std::abs(flow_rate_local - flow_rate_local_prior);
                //std::cout << "residue =" << residue << std::endl;
                flow_rate_local_prior = flow_rate_local;
                Real relax_factor = 0.3;
                vel_nodeO_i_prior = (1.0 - relax_factor) * vel_nodeO_i_prior + relax_factor * U_nodeO;
            }

            //** Check results *
            if (!std::isfinite(node_value_[index_i][0]))
            {
                if (*physical_time_ < start_time_laminar_)
                {
                    std::cout << "Warning: initial NaN in sublayer solver\n";
                    node_value_[index_i][0] = 0.0;
                }
                else
                {
                    std::cout << "Warning: unexpected NaN in sublayer solver\n";
                    
                    std::cout << "x is NaN\n";
                    std::cout << "u_outer=" << u_outer << std::endl;
                    std::cout << "k_outer=" << k_outer << std::endl;
                    std::cout << "omega_outer=" << omega_outer << std::endl;
                    std::cout << "dudn_outer=" << dudn_outer << std::endl;
                    std::cout << "nut_outer=" << nut_outer << std::endl;
                    std::cout << "distance_to_wall=" << distance_to_wall << std::endl;
                    std::cout << "friction_vel_magnitude_outer=" << friction_vel_magnitude_outer << std::endl;
                    std::cout << "flow_rate_local=" << flow_rate_local << std::endl;
                    std::cout << "node_value_[index_i][0] = " << node_value_[index_i][0] << std::endl;

                    node_value_[index_i][0] = 0.0;
                }
            }

            //** Extract results *
            friction_velocity_from_sublayer_[index_i] = node_value_[index_i][0];
            vel_nodeO_[index_i] = U_nodeO;
            vel_nodeUM_[index_i] = U_nodeUM;

            //** For testing *
            target_flow_rate_in_sublayer_[index_i] = flow_rate_local;
            vel_ps_magnitude_[index_i] = U_nodeO + dudn_outer * 0.5 * fluid_particle_spacing_;
            dudn_for_local_flow_rate_[index_i] = dudn_outer;
            utau_node_[index_i] = friction_vel_magnitude_outer;


            //** If calculate local gradient for correcting SPH solver *
            if(0)
            { 
                Vecd vel_tangential = vel_[index_i] - vel_[index_i].dot(normal) * normal;
                Real tangential_velocity_P_magnitude = vel_tangential.norm();
                Vecd tangential = vel_tangential / (tangential_velocity_P_magnitude + TinyReal);
                Real tangential_velocity_node_U = node_value_[index_i][5]; //** Temporary treatment *
                Real dist_nodeU_P = distance_to_wall / (double(num_sub_node_) + 0.5);
                Real dUdn_P_sublayer_magnitude = std::abs(tangential_velocity_P_magnitude - tangential_velocity_node_U) / (dist_nodeU_P + TinyReal);
                //Real dUdn_P_sublayer_magnitude = std::abs(node_value_[index_i][5] - node_value_[index_i][1]) / (4.0 * dist_nodeU_P + TinyReal);
                Matd dUdn_P_sublayer = dUdn_P_sublayer_magnitude * (tangential * normal.transpose());
                dUdn_P_sublayer_magnitude_[index_i] = dUdn_P_sublayer_magnitude;
                dUdn_P_sublayer_[index_i] = dUdn_P_sublayer;
            }
        }
    }
    //=================================================================================================//
    Vec6d P_refinement::solve_1D_sublayer(double kinematic_viscosity, double u_p_outer, double k_p_outer, 
        double w_p_outer, double vel_grad_p_outer, double nut_p_outer, double h_sublayer, double utau_outer, 
        double Q_target, double k_grad_p_outer, double w_grad_p_outer, double& vel_nodeO, double& vel_nodeUM)
    {
        //------------------------------------------------¡ý Input parameters ¡ý------------------------------------------------
        constexpr int ny = 5; // Manually determine

        if (num_sub_node_ != ny)
        {
            std::cout << "Node number mismatch! Stop here." << std::endl;
            std::cin.get();
        }

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
        double alpha = 0.9;

        double flow_rate_target = Q_target;
        double utau = utau_init;
        //------------------------------------------------¡ü Input parameters ¡ü------------------------------------------------

        //------------------------------------------------¡ý Node arrangement, for sublayer ¡ý------------------------------------------------
        double hy = height_sublayer / double(ny); // distance from node U to P_outer is hy/2
        double y_p = 0.5 * hy;
        double y[ny];  //computational nodes
        for (int i = 0; i < ny; ++i)
        {
            y[i] = y_p + i * hy;
        }
        //------------------------------------------------¡ü Node arrangement, for sublayer ¡ü------------------------------------------------

        //------------------------------------------------¡ý Calculate nodeP value ¡ý------------------------------------------------
        //** These values need tobe updated in each iteration *
        double yplus = utau * y_p / nu; 
        double u_p = utau * yplus;
        double turbu_omega_p = 6.0 * nu / (std_kw_beta_i_ * y_p * y_p);
        //------------------------------------------------¡ü Calculate nodeP value ¡ü------------------------------------------------

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
        double u_init_value[ny];
        double k_init_value[ny];
        double turbu_omega_init_value[ny];
        std::fill_n(u_init_value, ny, u_init);
        std::fill_n(k_init_value, ny, k_init);
        std::fill_n(turbu_omega_init_value, ny, turbu_omega_init);
        //------------------------------------------------¡ü Construct initial value ¡ü------------------------------------------------

        //------------------------------------------------¡ý Calculate nodeO and nodeUM ¡ý------------------------------------------------
        //** These values need tobe updated in each iteration *
        double u_nodeUM = std::max((u_init_value[ny - 1] + vel_grad_p_outer * hy), tiny);
        double k_nodeUM = std::max((k_init_value[ny - 1] + k_grad_p_outer * hy), tiny);
        double w_nodeUM = std::max((turbu_omega_init_value[ny - 1] + w_grad_p_outer * hy), tiny);
        double w_tilde_nodeUM = std::max(w_nodeUM, std_kw_C_lim_ * vel_grad_p_outer / std_kw_beta_star_5_); //** For calculating limiter of nut, assume vel_grad_p_outer = vel_grad_ndoeUM*
        double nut_nodeUM = k_nodeUM / (w_tilde_nodeUM + tiny);
        double u_nodeO = std::max((u_init_value[ny - 1] + vel_grad_p_outer * 0.5 * hy), tiny);
        double k_nodeO = std::max((k_init_value[ny - 1] + k_grad_p_outer * 0.5 * hy), tiny);
        double w_nodeO = std::max((turbu_omega_init_value[ny - 1] + w_grad_p_outer * 0.5 * hy), tiny);
        double w_tilde_nodeO = std::max(w_nodeO, std_kw_C_lim_ * vel_grad_p_outer / std_kw_beta_star_5_);
        double nut_nodeO = k_nodeO / (w_tilde_nodeO + tiny);
        //------------------------------------------------¡ü Calculate nodeO and nodeUM ¡ü------------------------------------------------

        //------------------------------------------------¡ý Construct solution vector ¡ý------------------------------------------------
        double phi_current[3 * ny];
        double phi_solved[3 * ny];
        for (int j = 0; j < ny; ++j)
        {
            phi_current[j] = u_init_value[j];
            phi_current[j + ny] = k_init_value[j];
            phi_current[j + 2 * ny] = turbu_omega_init_value[j];
        }
        std::copy_n(phi_current, 3 * ny, phi_solved);
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
            double* u_star = phi_solved;
            double* k_star = phi_solved + ny;
            double* turbu_omega_star = phi_solved + 2 * ny;
            //------------------------------------------------¡ü Update the star values ¡ü------------------------------------------------

            //------------------------------------------------¡ý Update boundary values ¡ý------------------------------------------------
            if(num_iter_out!=0) 
            {
            //-------------------------¡ý Update P value ¡ý-------------------------
            yplus = utau * y_p / nu;
            u_p = utau * yplus;
            turbu_omega_p = 6.0 * nu / (std_kw_beta_i_ * y_p * y_p);
            //-------------------------¡ü Update P value ¡ü-------------------------
            //
            //-------------------------¡ý Calculate nodeO and nodeUM ¡ý-------------------------
            u_nodeUM = std::max((u_star[ny - 1] + vel_grad_p_outer * hy), tiny);
            //u_nodeUM = 3.277959e-1;
            k_nodeUM = std::max((k_star[ny - 1] + k_grad_p_outer * hy), tiny);
            //k_nodeUM = 1.495287e-3;
            w_nodeUM = std::max((turbu_omega_star[ny - 1] + w_grad_p_outer * hy), tiny);
            w_tilde_nodeUM = std::max(w_nodeUM, std_kw_C_lim_ * vel_grad_p_outer / std_kw_beta_star_5_); //** For calculating limiter of nut, assume vel_grad_p_outer = vel_grad_ndoeUM*
            nut_nodeUM = k_nodeUM / (w_tilde_nodeUM + tiny);
            u_nodeO = std::max((u_star[ny - 1] + vel_grad_p_outer * 0.5 * hy), tiny);
            k_nodeO = std::max((k_star[ny - 1] + k_grad_p_outer * 0.5 * hy), tiny);
            w_nodeO = std::max((turbu_omega_star[ny - 1] + w_grad_p_outer * 0.5 * hy), tiny);
            w_tilde_nodeO = std::max(w_nodeO, std_kw_C_lim_ * vel_grad_p_outer / std_kw_beta_star_5_);
            nut_nodeO = k_nodeO / (w_tilde_nodeO + tiny);
            //-------------------------¡ü Calculate nodeO and nodeUM ¡ü-------------------------
            }
            //------------------------------------------------¡ü Update boundary values ¡ü------------------------------------------------

            //------------------------------------------------¡ý Calculate Dk, Dw, C_su ¡ý------------------------------------------------
            double diffusion_coefficient_k[ny]{};
            double diffusion_coefficient_turbu_omega[ny]{};
            double C_su[ny]{};
            double tau_over_rho_outer = (nu + nut_nodeO) * vel_grad_p_outer;  //** Here, nodeO value is used not nodeUM, since this is from anlytical *
            for (int i = 0; i < ny; ++i) {
                diffusion_coefficient_k[i] = nu + std_kw_sigma_star_ * k_star[i] / (turbu_omega_star[i] + tiny);
                diffusion_coefficient_turbu_omega[i] = nu + std_kw_sigma_ * k_star[i] / (turbu_omega_star[i] + tiny);
                C_su[i] = (utau * utau - tau_over_rho_outer) * (1.0 - y[i] / height_sublayer) + tau_over_rho_outer;
            }
            //------------------------------------------------¡ü Calculate Dk, Dw  ¡ü------------------------------------------------

            //------------------------------------------------¡ý Calculate gradients of u, k, omega, Dk, Dw ¡ý------------------------------------------------
            double dudy_discretized[ny]{};
            double dkdy[ny]{};
            double dwdy[ny]{};
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
            dudy_discretized[ny - 1] = (u_nodeUM - u_star[ny - 2]) / (2.0 * hy);
            dkdy[ny - 1] = (k_nodeUM - k_star[ny - 2]) / (2.0 * hy);
            dwdy[ny - 1] = (w_nodeUM - turbu_omega_star[ny - 2]) / (2.0 * hy);
            //------------------------------------------------¡ü Calculate gradients of u, k, omega, Dk, Dw ¡ü------------------------------------------------

            //------------------------------------------------¡ý Calculate nut_star ¡ý------------------------------------------------
            double turbu_omega_tilde[ny]{};
            double nut_star[ny]{};
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
            double dudy_star[ny]{};
            for (int i = 0; i < ny; ++i)
            {
                dudy_star[i] = C_su[i] / (nu + nut_star[i]);
            }
            //------------------------------------------------¡ü Calculate analytical gradient of u ¡ü------------------------------------------------

            //------------------------------------------------¡ý Update linearized source terms Sc Sp ¡ý------------------------------------------------
            // 
            //-------------------------------------¡ý For turbulent specific dissipation ¡ý-------------------------------------
            // calcualte part_cross_diffusion 
            double grad_prod[ny]{};
            double part_cross_diffusion[ny]{};
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
            double a_u[ny]{};
            double b_u[ny]{};
            double c_u[ny]{};
            double d_u[ny]{};
            // inner node
            for (int i = 1; i < ny - 1; ++i)
            {
                double nu_eff_i = nu + nut_star[i];
                double nu_eff_i_plus = nu + nut_star[i + 1];
                double nu_eff_i_minus = nu + nut_star[i - 1];
                double nu_eff_i_plus_half = 2.0 * nu_eff_i_plus * nu_eff_i / std::max((nu_eff_i_plus + nu_eff_i), tiny);
                double nu_eff_i_minus_half = 2.0 * nu_eff_i_minus * nu_eff_i / std::max((nu_eff_i_minus + nu_eff_i), tiny);
                a_u[i] = -nu_eff_i_minus_half;
                b_u[i] = std::max((nu_eff_i_plus_half + nu_eff_i_minus_half), tiny);
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
            double nu_eff_last_plus = nu + nut_nodeUM; // B.C.
            double nu_eff_last_minus = nu + nut_star[last - 1];
            double nu_eff_last_plus_half = 2.0 * nu_eff_last_plus * nu_eff_last / std::max((nu_eff_last_plus + nu_eff_last), tiny);
            double nu_eff_last_minus_half = 2.0 * nu_eff_last_minus * nu_eff_last / std::max((nu_eff_last_minus + nu_eff_last), tiny);
            a_u[last] = -nu_eff_last_minus_half;
            b_u[last] = std::max((nu_eff_last_plus_half + nu_eff_last_minus_half), tiny);
            c_u[last] = 0.0;
            d_u[last] = (utau * utau - tau_over_rho_outer) / height_sublayer * hy * hy + nu_eff_last_plus_half * u_nodeUM;
            // solving
            double U_new[ny]{};
            tdma5(a_u, b_u, c_u, d_u, U_new);
            //-------------------------------------¡ü For velocity ¡ü-------------------------------------

            //-------------------------------------¡ý For turbulent kinetic energy ¡ý-------------------------------------
            double a_k[ny]{};
            double b_k[ny]{};
            double c_k[ny]{};
            double d_k[ny]{};
            // inner node
            for (int i = 1; i < ny - 1; ++i) {
                double Dk_i = diffusion_coefficient_k[i];
                double Dk_i_plus = diffusion_coefficient_k[i + 1];
                double Dk_i_minus = diffusion_coefficient_k[i - 1];
                double Dk_i_plus_half = 2.0 * Dk_i_plus * Dk_i / std::max((Dk_i_plus + Dk_i), tiny);
                double Dk_i_minus_half = 2.0 * Dk_i_minus * Dk_i / std::max((Dk_i_minus + Dk_i), tiny);
                a_k[i] = -1.0 * Dk_i_minus_half;
                b_k[i] = Dk_i_plus_half + Dk_i_minus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[i];
                c_k[i] = -1.0 * Dk_i_plus_half;
                d_k[i] = hy * hy * nut_star[i] * dudy_star[i] * dudy_star[i];
            }
            // first node, grad k = 0
            double Dk_first = diffusion_coefficient_k[0];
            double Dk_first_plus = diffusion_coefficient_k[0 + 1];
            double Dk_first_plus_half = 2.0 * Dk_first_plus * Dk_first / std::max((Dk_first_plus + Dk_first), tiny);
            a_k[0] = 0.0;
            b_k[0] = Dk_first_plus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[0];
            c_k[0] = -1.0 * Dk_first_plus_half;
            d_k[0] = hy * hy * nut_star[0] * dudy_star[0] * dudy_star[0];
            // last node
            last = ny - 1;
            double Dk_last = diffusion_coefficient_k[last];
            double Dk_last_plus = nu + std_kw_sigma_star_ * k_nodeUM / (w_nodeUM + tiny);
            double Dk_last_minus = diffusion_coefficient_k[last - 1];
            double Dk_last_plus_half = 2.0 * Dk_last_plus * Dk_last / std::max((Dk_last_plus + Dk_last), tiny);
            double Dk_last_minus_half = 2.0 * Dk_last_minus * Dk_last / std::max((Dk_last_minus + Dk_last), tiny);
            a_k[last] = -1.0 * Dk_last_minus_half;
            b_k[last] = Dk_last_plus_half + Dk_last_minus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[last];
            c_k[last] = 0.0;
            d_k[last] = hy * hy * nut_star[last] * dudy_star[last] * dudy_star[last] + Dk_last_plus_half * k_nodeUM;
            // solving
            double K_new[ny]{};
            tdma5(a_k, b_k, c_k, d_k, K_new);
            // avoid negative value, K_new = max(K_new, k_min)
            double k_min = 1e-10;
            for (int i = 0; i < ny; ++i) {
                K_new[i] = std::max(K_new[i], k_min);
            }
            //-------------------------------------¡ü For turbulent kinetic energy ¡ü-------------------------------------

            //-------------------------------------¡ý For turbulent specific dissipation ¡ý-------------------------------------
            double a_w[ny]{};
            double b_w[ny]{};
            double c_w[ny]{};
            double d_w[ny]{};
            // inner node i = 1 ... ny-2
            for (int i = 1; i < ny - 1; ++i) {
                double Dw_i = diffusion_coefficient_turbu_omega[i];
                double Dw_i_plus = diffusion_coefficient_turbu_omega[i + 1];
                double Dw_i_minus = diffusion_coefficient_turbu_omega[i - 1];
                double Dw_i_plus_half = 2.0 * Dw_i_plus * Dw_i / std::max((Dw_i_plus + Dw_i), tiny);
                double Dw_i_minus_half = 2.0 * Dw_i_minus * Dw_i / std::max((Dw_i_minus + Dw_i), tiny);
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
            double Dw_last_plus = nu + std_kw_sigma_ * k_nodeUM / (w_nodeUM + tiny);
            double Dw_last_minus = diffusion_coefficient_turbu_omega[last - 1];
            double Dw_last_plus_half = 2.0 * Dw_last_plus * Dw_last / std::max((Dw_last_plus + Dw_last), tiny);
            double Dw_last_minus_half = 2.0 * Dw_last_minus * Dw_last / std::max((Dw_last_minus + Dw_last), tiny);
            a_w[last] = -1.0 * Dw_last_minus_half;
            b_w[last] = Dw_last_plus_half + Dw_last_minus_half + hy * hy * std_kw_beta_ * turbu_omega_star[last];
            c_w[last] = 0.0;
            double part_production_last = std_kw_alpha_ * turbu_omega_star[last] / (k_star[last] + tiny) * nut_star[last] * dudy_star[last] * dudy_star[last];
            d_w[last] = hy * hy * (part_production_last + part_cross_diffusion[last]) + Dw_last_plus_half * w_nodeUM;
            // solving
            double Turbu_omega_new[ny]{};
            tdma5(a_w, b_w, c_w, d_w, Turbu_omega_new);
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
            flow_rate_current = std::accumulate(U_new, U_new + ny, 0.0) * hy;
            //flow_rate_current -= U_new[ny-1] * hy * 0.5;

            // flow control
            double ratio = flow_rate_target / (flow_rate_current + 1e-12);
            
            // smooth limiter
            ratio = 0.5 * (ratio + std::sqrt(ratio * ratio + 1e-12));
            // optional upper bound
            ratio = std::min(ratio, 10.0);
            
            double utau_new = utau * std::sqrt(ratio);
            utau = (1.0 - alpha) * utau + alpha * utau_new;
            /*double error = U_new[ny-1] - u_p_outer;
            utau -= 0.1 * error;*/
            if (!std::isfinite(utau)) 
            {
                utau = utau_init;   // fallback
            }
            utau = std::max(utau, tiny);
            //------------------------------------------------¡ü Check and update flow rate ¡ü------------------------------------------------

            //std::cout << "flow_rate_current = " << flow_rate_current
            //    << ", target = " << flow_rate_target << std::endl;
            //std::cout << "updated utau = " << utau << std::endl;

            //------------------------------------------------¡ý Calculate residue ¡ý------------------------------------------------
            differ = 0.0;
            for (int i = 0; i < 3 * ny; ++i) {
                double diff = phi_current[i] - phi_solved[i];
                differ += diff * diff;
            }
            differ = std::sqrt(differ);
            //------------------------------------------------¡ü Calculate residue ¡ü------------------------------------------------
            //std::cout << "differ: " << differ << std::endl;

            //------------------------------------------------¡ý Update ¡ý------------------------------------------------
            for (int i = 0; i < 3 * ny; ++i) {
                phi_current[i] = phi_solved[i];
            }
            num_iter_out += 1;
            //------------------------------------------------¡ü Update ¡ü------------------------------------------------

            //std::cout << "num_iter_out = " << num_iter_out << std::endl;
            //std::cout << "------------" << std::endl;

        }

        ///*
        //** This is a temporary treatment *
        if (ny != 5)
        {
            std::cout << "ny is not 5, currently not allowed! Stop here." << std::endl;
            std::cin.get();
        }
        Vec6d results = Vec6d::Zero();
        results[0] = utau;
        for (int i = 0; i < ny; ++i) {
            results[i+1] = phi_solved[i];
        }
        vel_nodeO = u_nodeO;
        vel_nodeUM = u_nodeUM;
        return results;
        //*/

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

        // ================== Extract Solution ==================
        n_start = 0;

        std::vector<double> U(ny);
        for (int i = 0; i < ny; ++i) U[i] = phi_solved[n_start + i];

        n_start += ny;
        std::vector<double> K(ny);
        for (int i = 0; i < ny; ++i) K[i] = phi_solved[n_start + i];

        n_start += ny;
        std::vector<double> OMEGA(ny);
        for (int i = 0; i < ny; ++i) OMEGA[i] = phi_solved[n_start + i];

        // ================== use existing y ==================
        //std::vector<double> Y = y;

        // ================== calculate nut ==================
        double eps = 1e-12;
        std::vector<double> NUT(ny);
        for (int i = 0; i < ny; ++i) {
            NUT[i] = K[i] / (OMEGA[i] + eps);
        }
        
        // ================== Êä³ö Tecplot ÎÄ¼þ ==================
        int NF = 40;
        int index = 76;
        std::string header_line = "ZONE T=\"SPH(1D)46 NF="
            + std::to_string(NF)
            + " ("
            + std::to_string(index)
            + ")\"";
        std::string filename = "pipe_kw_nf" + std::to_string(NF) + "_" + std::to_string(index) + ".dat";

        std::ofstream fout(filename);
        fout << "$VARIABLES = \"Y\", \"U\", \"K\", \"OMEGA\", \"NUT(k/omega)\"\n";
        fout << "$friction velocity = " << utau << "\n";
        fout << "$current flow rate = " << flow_rate_current << "\n";
        fout << "$num_iter_out = " << num_iter_out << "\n";
        fout << header_line << "\n";

        fout << std::scientific << std::setprecision(8);

        for (int i = 0; i < ny; ++i) {
            fout << y[i] << " "
                << U[i] << " "
                << K[i] << " "
                << OMEGA[i] << " "
                << NUT[i] << "\n";
        }

        fout.close();

        std::cout << "Tecplot has been created: " << filename << std::endl;

        std::cin.get();
    }
    //=================================================================================================//
    // ================= TDMA =================
    // Solve a tridiagonal system: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
    // a[0] must be 0, c[n-1] will be ignored
    void P_refinement::tdma(int N, const double* a, const double* b, const double* c, const double* d, double* x)
    {
        std::vector<double> cp(N, 0.0);
        std::vector<double> dp(N, 0.0);

        if (std::abs(b[0]) < 1e-14)
            throw std::runtime_error("TDMA: b[0] too small!");

        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];

        for (int i = 1; i < N; ++i)
        {
            double denom = b[i] - a[i] * cp[i - 1];
            if (std::abs(denom) < 1e-14)
                throw std::runtime_error("TDMA: denom too small");

            cp[i] = (i == N - 1) ? 0.0 : c[i] / denom;
            dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
        }

        x[N - 1] = dp[N - 1];
        for (int i = N - 2; i >= 0; --i)
        {
            x[i] = dp[i] - cp[i] * x[i + 1];
        }
    }
//=================================================================================================//
    // ================= TDMA5 =================
    // Solve a tridiagonal system of size 5:
    // a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
    // a[0] must be 0, c[4] will be ignored
    void P_refinement::tdma5(const double a[5], const double b[5], const double c[5], const double d[5], double x[5])
    {
        if (num_sub_node_ != 5)
        {
            std::cout << "TDMA5: Node number mismatch! Stop here." << std::endl;
            std::cin.get();
        }

        double cp[5]{ 0.0 }; // modified upper diagonal
        double dp[5]{ 0.0 }; // modified right-hand side

        // ---------------- Step 0: first row ----------------
        if (std::abs(b[0]) < 1e-14) throw std::runtime_error("TDMA5: b[0] too small!");
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];

        // ---------------- Forward sweep ----------------
        // i = 1
        {
            double denom = b[1] - a[1] * cp[0];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA5: denom too small at row 1");
            cp[1] = c[1] / denom;
            dp[1] = (d[1] - a[1] * dp[0]) / denom;
        }

        // i = 2
        {
            double denom = b[2] - a[2] * cp[1];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA5: denom too small at row 2");
            cp[2] = c[2] / denom;
            dp[2] = (d[2] - a[2] * dp[1]) / denom;
        }

        // i = 3
        {
            double denom = b[3] - a[3] * cp[2];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA5: denom too small at row 3");
            cp[3] = c[3] / denom;
            dp[3] = (d[3] - a[3] * dp[2]) / denom;
        }

        // i = 4
        {
            double denom = b[4] - a[4] * cp[3];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA5: denom too small at row 4");
            cp[4] = 0.0; // last row has no right neighbor
            dp[4] = (d[4] - a[4] * dp[3]) / denom;
        }

        // ---------------- Back substitution ----------------
        x[4] = dp[4];
        x[3] = dp[3] - cp[3] * x[4];
        x[2] = dp[2] - cp[2] * x[3];
        x[1] = dp[1] - cp[1] * x[2];
        x[0] = dp[0] - cp[0] * x[1];
    }
    // ================= TDMA10 =================
    // Solve a tridiagonal system of size 10:
    // a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
    // a[0] must be 0, c[9] will be ignored
    void P_refinement::tdma10(const double a[10], const double b[10], const double c[10], const double d[10], double x[10])
    {
        if (num_sub_node_ != 10)
        {
            std::cout << "TDMA10: Node number mismatch! Stop here." << std::endl;
            std::cin.get();
        }

        double cp[10]{ 0.0 }; // modified upper diagonal
        double dp[10]{ 0.0 }; // modified RHS

        // ---------------- Step 0 ----------------
        if (std::abs(b[0]) < 1e-14) throw std::runtime_error("TDMA10: b[0] too small!");
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];

        // ---------------- Forward sweep ----------------
        // i = 1
        {
            double denom = b[1] - a[1] * cp[0];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 1");
            cp[1] = c[1] / denom;
            dp[1] = (d[1] - a[1] * dp[0]) / denom;
        }

        // i = 2
        {
            double denom = b[2] - a[2] * cp[1];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 2");
            cp[2] = c[2] / denom;
            dp[2] = (d[2] - a[2] * dp[1]) / denom;
        }

        // i = 3
        {
            double denom = b[3] - a[3] * cp[2];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 3");
            cp[3] = c[3] / denom;
            dp[3] = (d[3] - a[3] * dp[2]) / denom;
        }

        // i = 4
        {
            double denom = b[4] - a[4] * cp[3];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 4");
            cp[4] = c[4] / denom;
            dp[4] = (d[4] - a[4] * dp[3]) / denom;
        }

        // i = 5
        {
            double denom = b[5] - a[5] * cp[4];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 5");
            cp[5] = c[5] / denom;
            dp[5] = (d[5] - a[5] * dp[4]) / denom;
        }

        // i = 6
        {
            double denom = b[6] - a[6] * cp[5];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 6");
            cp[6] = c[6] / denom;
            dp[6] = (d[6] - a[6] * dp[5]) / denom;
        }

        // i = 7
        {
            double denom = b[7] - a[7] * cp[6];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 7");
            cp[7] = c[7] / denom;
            dp[7] = (d[7] - a[7] * dp[6]) / denom;
        }

        // i = 8
        {
            double denom = b[8] - a[8] * cp[7];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 8");
            cp[8] = c[8] / denom;
            dp[8] = (d[8] - a[8] * dp[7]) / denom;
        }

        // i = 9 (last row)
        {
            double denom = b[9] - a[9] * cp[8];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 9");
            cp[9] = 0.0;
            dp[9] = (d[9] - a[9] * dp[8]) / denom;
        }

        // ---------------- Back substitution ----------------
        x[9] = dp[9];
        x[8] = dp[8] - cp[8] * x[9];
        x[7] = dp[7] - cp[7] * x[8];
        x[6] = dp[6] - cp[6] * x[7];
        x[5] = dp[5] - cp[5] * x[6];
        x[4] = dp[4] - cp[4] * x[5];
        x[3] = dp[3] - cp[3] * x[4];
        x[2] = dp[2] - cp[2] * x[3];
        x[1] = dp[1] - cp[1] * x[2];
        x[0] = dp[0] - cp[0] * x[1];
    }

    //=============================================================================================//
    void BodyStatesRecordingToVtpIncludeNode::writeWithFileName(const std::string& sequence)
    {
        for (SPHBody* body : bodies_)
        {
            if (body->checkNewlyUpdated())
            {
                BaseParticles& base_particles = body->getBaseParticles();

                if (state_recording_)
                {
                    std::string filefullpath = io_environment_.OutputFolder() + "/" + body->getName() + "_" + sequence + ".vtp";
                    if (fs::exists(filefullpath))
                    {
                        fs::remove(filefullpath);
                    }
                    std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
                    // begin of the XML file
                    out_file << "<?xml version=\"1.0\"?>\n";
                    out_file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
                    out_file << " <PolyData>\n";

                    // physical time
                    if (sph_system_.isPhysical())
                    {
                        out_file << "<FieldData>\n";
                        out_file << "<DataArray type=\"Float64\"  Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">\n";
                        out_file << std::fixed << std::setprecision(9) << sv_physical_time_->getValue() << "\n";
                        out_file << " </DataArray>\n";
                        out_file << "</FieldData>\n";
                    }

                    size_t total_real_particles = base_particles.TotalRealParticles();
                    out_file << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles
                        << "\" NumberOfVerts=\"" << total_real_particles << "\">\n";

                    // write current/final particle positions first
                    out_file << "   <Points>\n";
                    out_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
                    out_file << "    ";
                    for (size_t i = 0; i != total_real_particles; ++i)
                    {
                        Vec3d particle_position = upgradeToVec3d(base_particles.ParticlePositions()[i]);
                        out_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
                    }
                    out_file << std::endl;
                    out_file << "    </DataArray>\n";
                    out_file << "   </Points>\n";

                    // write header of particles data
                    out_file << "   <PointData  Vectors=\"vector\">\n";
                    writeParticlesToVtk(out_file, base_particles);
                    out_file << "   </PointData>\n";

                    // write empty cells
                    out_file << "   <Verts>\n";
                    out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
                    out_file << "    ";
                    for (size_t i = 0; i != total_real_particles; ++i)
                    {
                        out_file << i << " ";
                    }
                    out_file << std::endl;
                    out_file << "    </DataArray>\n";
                    out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
                    out_file << "    ";
                    for (size_t i = 0; i != total_real_particles; ++i)
                    {
                        out_file << i + 1 << " ";
                    }
                    out_file << std::endl;
                    out_file << "    </DataArray>\n";
                    out_file << "   </Verts>\n";

                    out_file << "  </Piece>\n";
                    out_file << " </PolyData>\n";
                    out_file << "</VTKFile>\n";

                    out_file.close();
                }
            }
            body->setNotNewlyUpdated();
        }
    }
} // namespace udf
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//