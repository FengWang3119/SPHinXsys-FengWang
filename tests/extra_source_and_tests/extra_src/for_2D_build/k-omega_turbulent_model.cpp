//#pragma once
#include "k-omega_turbulent_model.hpp"
namespace SPH
{
//=================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
kOmega_BaseTurbuClosureCoeff::kOmega_BaseTurbuClosureCoeff()
    : std_kw_beta_star_(0.09), std_kw_sigma_star_(0.6),
      std_kw_alpha_(0.52), std_kw_sigma_(0.5), std_kw_f_beta_(1.0), std_kw_beta_0_(0.0708),
      std_kw_sigma_do_(0.125), std_kw_C_lim_(0.875), std_kw_beta_i_(0.075)
{
    std_kw_beta_ = std_kw_beta_0_ * std_kw_f_beta_;
    std_kw_beta_star_25_ = pow(std_kw_beta_star_, 0.25);
}
//=================================================================================================//
kOmegaTurbulentEddyViscosity::
    kOmegaTurbulentEddyViscosity(SPHBody &sph_body)
    : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
      rho_(*particles_->getVariableDataByName<Real>("Density")),
      turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(*particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
      wall_Y_plus_(*particles_->getVariableDataByName<Real>("WallYplus")),
      wall_Y_star_(*particles_->getVariableDataByName<Real>("WallYstar")),
      turbu_strain_rate_traceless_magnitude_(*particles_->getVariableDataByName<Real>("TurbulentStrainRateTracelessMagnitude")),
      mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()) {}
//=================================================================================================//
void kOmegaTurbulentEddyViscosity::update(size_t index_i, Real dt)
{
    Real limited_omega = std_kw_C_lim_ * turbu_strain_rate_traceless_magnitude_[index_i] / sqrt(std_kw_beta_star_);
    Real turbu_omega_tilde_ = SMAX(turbu_omega_[index_i], limited_omega);
    turbu_mu_[index_i] = rho_[index_i] * turbu_k_[index_i] / turbu_omega_tilde_;
}
//=================================================================================================//
kOmegaStdWallFuncCorrection::
    kOmegaStdWallFuncCorrection(BaseInnerRelation &inner_relation,
                                BaseContactRelation &contact_relation, Real y_p_constant)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateContact(contact_relation),
      y_p_(*particles_->registerSharedVariable<Real>("Y_P")),
      wall_Y_plus_(*particles_->registerSharedVariable<Real>("WallYplus")),
      wall_Y_star_(*particles_->registerSharedVariable<Real>("WallYstar")),
      velo_tan_(*particles_->registerSharedVariable<Real>("TangentialVelocity")),
      velo_friction_(*particles_->registerSharedVariable<Vecd>("FrictionVelocity")),
      vel_(*particles_->getVariableDataByName<Vecd>("Velocity")), rho_(*particles_->getVariableDataByName<Real>("Density")),
      molecular_viscosity_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()),
      turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(*particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
      turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      is_near_wall_P1_(*particles_->getVariableDataByName<int>("IsNearWallP1")),
      is_near_wall_P2_(*particles_->getVariableDataByName<int>("IsNearWallP2")),
      velocity_gradient_(*particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
      k_production_(*particles_->getVariableDataByName<Real>("K_Production")),
      distance_to_dummy_interface_(*particles_->getVariableDataByName<Real>("DistanceToDummyInterface")),
      distance_to_dummy_interface_up_average_(*particles_->getVariableDataByName<Real>("DistanceToDummyInterfaceUpAver")),
      index_nearest(*particles_->getVariableDataByName<int>("NearestIndex")),
      e_nearest_tau_(*particles_->getVariableDataByName<Vecd>("WallNearestTangentialUnitVector")),
      e_nearest_normal_(*particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_n_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("NormalDirection"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }

    //particles_->registerSharedVariable(y_p_, "Y_P");
    particles_->addVariableToSort<Real>("Y_P");
    particles_->addVariableToWrite<Real>("Y_P");

    //** Fixed y_p_ as a constant distance *
    std::fill(y_p_.begin(), y_p_.end(), y_p_constant);

    //particles_->registerSharedVariable(wall_Y_plus_, "WallYplus");
    particles_->addVariableToSort<Real>("WallYplus");
    particles_->addVariableToWrite<Real>("WallYplus");

    //** Initial value is important, especially when use log law *
    //particles_->registerSharedVariable(wall_Y_star_, "WallYstar", TinyReal);
    particles_->addVariableToSort<Real>("WallYstar");
    particles_->addVariableToWrite<Real>("WallYstar");

    //particles_->registerSharedVariable(velo_tan_, "TangentialVelocity");
    particles_->addVariableToSort<Real>("TangentialVelocity");
    particles_->addVariableToWrite<Real>("TangentialVelocity");

    //particles_->registerSharedVariable(velo_friction_, "FrictionVelocity");
    particles_->addVariableToSort<Vecd>("FrictionVelocity");
    particles_->addVariableToWrite<Vecd>("FrictionVelocity");
};
//=================================================================================================//
void kOmegaStdWallFuncCorrection::interaction(size_t index_i, Real dt)
{
    velo_tan_[index_i] = 0.0;
    velo_friction_[index_i] = Vecd::Zero();
    wall_Y_plus_[index_i] = 0.0;
    wall_Y_star_[index_i] = 0.0;

    //** If use level-set to get distance from P to wall, activate this *
    //y_p_[index_i]= distance_to_dummy_interface_levelset_[index_i];

    if (is_near_wall_P2_[index_i] == 10)
    {
        Real y_p_constant_i = y_p_[index_i];

        Real turbu_k_i_05 = pow(turbu_k_[index_i], 0.5);
        //Real turbu_k_i_15 = pow(turbu_k_[index_i], 1.5);

        //** Choose one kind of the distance to calculate the wall-nearest values *
        //Real r_dummy_normal = distance_to_dummy_interface_up_average_[index_i];
        //Real r_dummy_normal = distance_to_dummy_interface_[index_i];
        //Real r_dummy_normal = distance_to_dummy_interface_levelset_[index_i];

        Vecd e_i_nearest_tau = e_nearest_tau_[index_i];
        Vecd e_i_nearest_n = e_nearest_normal_[index_i];
        const Vecd &vel_i = vel_[index_i];
        Real rho_i = rho_[index_i];
        Real nu_i = molecular_viscosity_ / rho_i;

        //** Calculate Y_star, note the current code is based on Y_star *
        wall_Y_star_[index_i] = y_p_constant_i * C_mu_25_ * turbu_k_i_05 / nu_i;

        //** Calculate friction velocity, including P2 region. *
        Real velo_fric_mag = 0.0;
        Real velo_tan_mag = 0.0; //** tangential velo magnitude for fluid particle i *

        velo_tan_mag = abs(e_i_nearest_tau.dot(vel_i));
        velo_tan_[index_i] = velo_tan_mag;

        if (wall_Y_star_[index_i] != static_cast<Real>(wall_Y_star_[index_i]))
        {
            std::cout << "y* is not a real value, please check" << std::endl;
            std::cin.get();
        }

        Real u_star = get_dimensionless_velocity(wall_Y_star_[index_i]);
        velo_fric_mag = sqrt(C_mu_25_ * turbu_k_i_05 * velo_tan_mag / u_star);

        if (velo_fric_mag != static_cast<Real>(velo_fric_mag))
        {
            std::cout << "friction velocity is not a real, please check" << std::endl;
            std::cout << "velo_fric=" << velo_fric_mag << std::endl
                      << "velo_tan_mag=" << velo_tan_mag << std::endl;
            std::cout << "turbu_k_=" << pow(turbu_k_[index_i], 0.5) << std::endl;
            std::cout << "sum=" << (Karman_ * velo_tan_mag * C_mu_25_ * pow(turbu_k_[index_i], 0.5) / log(turbu_const_E_ * C_mu_25_ * pow(turbu_k_[index_i], 0.5) * y_p_constant_i * rho_i / molecular_viscosity_)) << std::endl;
            std::cout << "numerator=" << Karman_ * velo_tan_mag * C_mu_25_ * pow(turbu_k_[index_i], 0.5) << std::endl;
            std::cout << "denominator=" << log(turbu_const_E_ * C_mu_25_ * pow(turbu_k_[index_i], 0.5) * y_p_constant_i * rho_i / molecular_viscosity_) << std::endl;
            Real temp = C_mu_25_ * pow(turbu_k_[index_i], 0.5) * velo_tan_mag / u_star;

            std::cout << "temp =" << temp << std::endl;

            std::cout << "pow(turbu_k_[index_i], 0.5) =" << pow(turbu_k_[index_i], 0.5) << std::endl;
            std::cout << "velo_tan_mag / u_star =" << velo_tan_mag / u_star << std::endl;
            std::cout << "velo_tan_mag =" << velo_tan_mag << std::endl;
            std::cout << " u_star =" << u_star << std::endl;
            std::cin.get();
        }

        //** friction velocity have the same direction of vel_i, if not, change its direction *
        velo_friction_[index_i] = velo_fric_mag * e_i_nearest_tau;
        if (vel_i.dot(velo_friction_[index_i]) < 0.0)
            velo_friction_[index_i] = -1.0 * velo_friction_[index_i];

        //** Calculate Y_plus  *
        wall_Y_plus_[index_i] = y_p_constant_i * velo_fric_mag / nu_i;

        // ** Correct the near wall values, only for P1 region *
        if (is_near_wall_P1_[index_i] == 1)
        {
            Matd vel_grad_i_tn = Matd::Zero(); //** velocity gradient of wall-nearest fluid particle i on t-n plane *
            Matd Q = Matd::Zero();
            Real total_weight = 0.0;

            Real dudn_p_weighted_sum = 0.0;
            Real G_k_p_weighted_sum = 0.0;
            Real omega_p_weighted_sum = 0.0;

            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
                StdLargeVec<Vecd> &n_k = *(contact_n_[k]);
                Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    Real dudn_p_j = 0.0;
                    Real G_k_p_j = 0.0;
                    Real omega_p_j = 0.0;

                    Real y_p_j = 0.0;

                    Vecd e_j_tau = Vecd::Zero();

                    size_t index_j = contact_neighborhood.j_[n];
                    Vecd e_j_n = n_k[index_j];

                    //** Get tangential unit vector, temporarily only suitable for 2D*
                    e_j_tau[0] = e_j_n[1];
                    e_j_tau[1] = e_j_n[0] * (-1.0);

                    y_p_j = get_distance_from_P_to_wall(y_p_constant_i);

                    Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
                    total_weight += weight_j;

                    Real denominator_log_law_j = C_mu_25_ * turbu_k_i_05 * Karman_ * y_p_j;

                    Real vel_i_tau_mag = abs(vel_i.dot(e_j_tau));
                    Real y_star_j = C_mu_25_ * turbu_k_i_05 * y_p_j / nu_i;
                    Real u_star_j = get_dimensionless_velocity(y_star_j);
                    Real fric_vel_mag_j = sqrt(C_mu_25_ * turbu_k_i_05 * vel_i_tau_mag / u_star_j);

                    Real dudn_p_mag_j = get_near_wall_velocity_gradient_magnitude(y_star_j, fric_vel_mag_j, denominator_log_law_j, nu_i);
                    dudn_p_j = dudn_p_mag_j * boost::qvm::sign(vel_i.dot(e_j_tau));
                    dudn_p_weighted_sum += weight_j * dudn_p_j;

                    if (y_star_j < y_star_threshold_laminar_ && GlobalStaticVariables::physical_time_ > start_time_laminar_)
                    {
                        G_k_p_j = 0.0;
                        omega_p_j = 6.0 * nu_i / (std_kw_beta_i_ * y_p_j * y_p_j);
                    }
                    else
                    {
                        G_k_p_j = rho_i * fric_vel_mag_j * fric_vel_mag_j * dudn_p_mag_j;
                        omega_p_j = turbu_k_i_05 / (std_kw_beta_star_25_ * Karman_ * y_p_j);
                    }
                    G_k_p_weighted_sum += weight_j * G_k_p_j;
                    omega_p_weighted_sum += weight_j * omega_p_j;
                }
            }
            turbu_omega_[index_i] = omega_p_weighted_sum / total_weight;

            vel_grad_i_tn(0, 0) = 0.0;
            vel_grad_i_tn(0, 1) = dudn_p_weighted_sum / total_weight;
            vel_grad_i_tn(1, 0) = 0.0;
            vel_grad_i_tn(1, 1) = 0.0;

            Q = getTransformationMatrix(e_i_nearest_n);

            velocity_gradient_[index_i] = Q.transpose() * vel_grad_i_tn * Q;

            k_production_[index_i] = G_k_p_weighted_sum / total_weight;
        }
    }
}
//=================================================================================================//
kOmega_kTransportEquationInner::kOmega_kTransportEquationInner(BaseInnerRelation &inner_relation, const StdVec<Real> &initial_values, int is_extr_visc_dissipa)
    : kOmega_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      dk_dt_(*particles_->registerSharedVariable<Real>("ChangeRateOfTKE")),
      dk_dt_without_dissipation_(*particles_->registerSharedVariable<Real>("ChangeRateOfTKEWithoutDissipation")),
      k_production_(*particles_->registerSharedVariable<Real>("K_Production")),
      is_near_wall_P1_(*particles_->getVariableDataByName<int>("IsNearWallP1")),
      velocity_gradient_(*particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
      turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(*particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
      turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_strain_rate_(*particles_->getVariableDataByName<Matd>("TurbulentStrainRate")),
      turbu_strain_rate_magnitude_(*particles_->getVariableDataByName<Real>("TurbulentStrainRateMagnitude")),
      is_extra_viscous_dissipation_(*particles_->registerSharedVariable<int>("TurbulentExtraViscousDissipation")),
      turbu_indicator_(*particles_->registerSharedVariable<int>("TurbulentIndicator")),
      k_diffusion_(*particles_->registerSharedVariable<Real>("K_Diffusion")),
      vel_x_(*particles_->registerSharedVariable<Real>("Velocity_X"))
{
    particles_->addVariableToSort<Real>("ChangeRateOfTKE");
    particles_->addVariableToSort<Real>("ChangeRateOfTKEWithoutDissipation");

    particles_->addVariableToSort<Real>("K_Production");
    particles_->addVariableToWrite<Real>("K_Production");

    particles_->addVariableToSort<Real>("TurbulenceKineticEnergy");
    particles_->addVariableToWrite<Real>("TurbulenceKineticEnergy");

    particles_->addVariableToSort<Real>("TurbulentSpecificDissipation");
    particles_->addVariableToWrite<Real>("TurbulentSpecificDissipation");

    particles_->addVariableToSort<Real>("TurbulentViscosity");
    particles_->addVariableToWrite<Real>("TurbulentViscosity");

    particles_->addVariableToSort<Matd>("TurbulentStrainRate");
    particles_->addVariableToWrite<Matd>("TurbulentStrainRate");

    //** Obtain Initial values for transport equations *
    std::fill(turbu_k_.begin(), turbu_k_.end(), initial_values[0]);
    std::fill(turbu_omega_.begin(), turbu_omega_.end(), initial_values[1]);
    std::fill(turbu_mu_.begin(), turbu_mu_.end(), initial_values[2]);

    //** for test */
    particles_->addVariableToSort<Real>("K_Diffusion");
    particles_->addVariableToWrite<Real>("K_Diffusion");

    particles_->addVariableToWrite<Real>("ChangeRateOfTKE");

    particles_->addVariableToSort<Real>("Velocity_X");

    particles_->addVariableToSort<int>("TurbulentIndicator");
    particles_->addVariableToWrite<int>("TurbulentIndicator");

    std::fill(is_extra_viscous_dissipation_.begin(), is_extra_viscous_dissipation_.end(), is_extr_visc_dissipa);
}
//=================================================================================================//
void kOmega_kTransportEquationInner::interaction(size_t index_i, Real dt)
{
    //Vecd vel_i = vel_[index_i];
    Real rho_i = rho_[index_i];
    Real turbu_mu_i = turbu_mu_[index_i];
    Real turbu_k_i = turbu_k_[index_i];
    Real turbu_omega_i = turbu_omega_[index_i];
    Real mu_eff_i = mu_ + std_kw_sigma_star_ * turbu_k_i / turbu_omega_i;

    dk_dt_[index_i] = 0.0;
    dk_dt_without_dissipation_[index_i] = 0.0;
    Real k_derivative(0.0);
    Real k_lap(0.0);
    Matd strain_rate = Matd::Zero();
    Matd strain_rate_traceless = Matd::Zero();
    Matd Re_stress = Matd::Zero();

    Real k_production(0.0);
    Real k_dissipation(0.0);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mu_eff_j = mu_ + std_kw_sigma_star_ * turbu_k_[index_j] / turbu_omega_[index_j];
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        k_derivative = (turbu_k_i - turbu_k_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        k_lap += 2.0 * mu_harmo * k_derivative * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
    }
    k_lap /= rho_i;

    strain_rate = 0.5 * (velocity_gradient_[index_i].transpose() + velocity_gradient_[index_i]);
    Real strain_rate_trace = strain_rate.trace();
    strain_rate_traceless = strain_rate - (1.0 / Dimensions) * strain_rate_trace * Matd::Identity(); //** For [2008 wilcox AIAA] *

    Real strain_rate_squire = (strain_rate.array() * strain_rate.array()).sum();
    turbu_strain_rate_magnitude_[index_i] = sqrt(2.0 * strain_rate_squire);
    Real strain_rate_traceless_squire = (strain_rate_traceless.array() * strain_rate_traceless.array()).sum();
    turbu_strain_rate_traceless_magnitude_[index_i] = sqrt(2.0 * strain_rate_traceless_squire);

    Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();
    //Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i;

    Matd k_production_matrix = Re_stress.array() * velocity_gradient_[index_i].array();
    //** The near wall k production is updated in wall function part *
    //** So this k_production variable should not be initialised *
    if (is_near_wall_P1_[index_i] != 1)
        k_production_[index_i] = k_production_matrix.sum();

    k_production = k_production_[index_i];
    k_dissipation = std_kw_beta_star_ * turbu_k_i * turbu_omega_i;

    dk_dt_[index_i] = k_production - k_dissipation + k_lap;
    dk_dt_without_dissipation_[index_i] = k_production + k_lap;

    //** for record */
    k_diffusion_[index_i] = k_lap;
    vel_x_[index_i] = vel_[index_i][0];
    turbu_strain_rate_[index_i] = strain_rate;
}
//=================================================================================================//
void kOmega_kTransportEquationInner::update(size_t index_i, Real dt)
{
    turbu_k_[index_i] += dk_dt_[index_i] * dt;
}
//=================================================================================================//
kOmega_omegaTransportEquationInner::kOmega_omegaTransportEquationInner(BaseInnerRelation &inner_relation)
    : kOmega_BaseTurbulentModel<Base, DataDelegateInner>(inner_relation),
      domega_dt_(*particles_->registerSharedVariable<Real>("ChangeRateOfTDR")),
      domega_dt_without_dissipation_(*particles_->registerSharedVariable<Real>("ChangeRateOfTDRWithoutDissp")),
      omega_production_(*particles_->registerSharedVariable<Real>("omega_Production")),
      omega_dissipation_(*particles_->registerSharedVariable<Real>("omega_Dissipation")),
      omega_diffusion_(*particles_->registerSharedVariable<Real>("omega_Diffusion")),
      omega_cross_diffusion_(*particles_->registerSharedVariable<Real>("omega_Cross_Diffusion")),
      turbu_mu_(*particles_->getVariableDataByName<Real>("TurbulentViscosity")),
      turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(*particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
      k_production_(*particles_->getVariableDataByName<Real>("K_Production")),
      is_near_wall_P1_(*particles_->getVariableDataByName<int>("IsNearWallP1"))
{
    particles_->addVariableToSort<Real>("ChangeRateOfTDR");
    particles_->addVariableToWrite<Real>("ChangeRateOfTDR");

    particles_->addVariableToSort<Real>("omega_Production");
    particles_->addVariableToWrite<Real>("omega_Production");

    particles_->addVariableToSort<Real>("omega_Dissipation");
    particles_->addVariableToWrite<Real>("omega_Dissipation");

    particles_->addVariableToSort<Real>("omega_Diffusion");
    particles_->addVariableToWrite<Real>("omega_Diffusion");

    particles_->addVariableToSort<Real>("omega_Cross_Diffusion");
    particles_->addVariableToWrite<Real>("omega_Cross_Diffusion");
}
//=================================================================================================//
void kOmega_omegaTransportEquationInner::
    interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    //Real turbu_mu_i = turbu_mu_[index_i];
    Real turbu_k_i = turbu_k_[index_i];
    Real turbu_omega_i = turbu_omega_[index_i];

    Real mu_eff_i = mu_ + std_kw_sigma_ * turbu_k_i / turbu_omega_i;

    domega_dt_[index_i] = 0.0;
    domega_dt_without_dissipation_[index_i] = 0.0;
    Real omega_production(0.0);
    Real omega_derivative(0.0);
    Real omega_lap(0.0);
    Real omega_dissipation(0.0);
    Real omega_cross_diffusion(0.0);
    //std_kw_alpha_[index_i] = get_alpha_standard_kw();
    Vecd k_gradient = Vecd::Zero();
    Vecd omega_gradient = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mu_eff_j = mu_ + std_kw_sigma_ * turbu_k_[index_j] / turbu_omega_[index_j];
        Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        omega_derivative = (turbu_omega_i - turbu_omega_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        omega_lap += 2.0 * mu_harmo * omega_derivative * inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] / rho_i;

        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
        //** non-conservative form *
        k_gradient += -1.0 * (turbu_k_i - turbu_k_[index_j]) * nablaW_ijV_j;
        omega_gradient += -1.0 * (turbu_omega_i - turbu_omega_[index_j]) * nablaW_ijV_j;
    }

    omega_production = std_kw_alpha_ * turbu_omega_i * k_production_[index_i] / turbu_k_i;
    omega_dissipation = std_kw_beta_ * turbu_omega_i * turbu_omega_i;

    Real grad_dot_k_omega = k_gradient.dot(omega_gradient);
    std_kw_sigma_d_ = grad_dot_k_omega > 0.0 ? std_kw_sigma_do_ : 0.0;

    omega_cross_diffusion = std_kw_sigma_d_ / turbu_omega_i * grad_dot_k_omega;

    domega_dt_[index_i] = omega_production - omega_dissipation + omega_lap + omega_cross_diffusion;
    domega_dt_without_dissipation_[index_i] = omega_production + omega_lap + omega_cross_diffusion;

    //** for test */
    omega_production_[index_i] = omega_production;
    omega_dissipation_[index_i] = omega_dissipation;
    omega_diffusion_[index_i] = omega_lap;
    omega_cross_diffusion_[index_i] = omega_cross_diffusion;
}
//=================================================================================================//
void kOmega_omegaTransportEquationInner::update(size_t index_i, Real dt)
{
    //** The near wall omega value is updated in wall function part *
    if (is_near_wall_P1_[index_i] != 1)
    {
        turbu_omega_[index_i] += domega_dt_[index_i] * dt;
    }
}
//=================================================================================================//
kOmegaInflowTurbulentCondition::kOmegaInflowTurbulentCondition(BodyPartByCell &body_part, Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet)
    : BaseFlowBoundaryCondition(body_part), type_turbu_inlet_(type_turbu_inlet),
      relaxation_rate_(relaxation_rate),
      CharacteristicLength_(CharacteristicLength),
      turbu_k_(*particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_omega_(*particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation"))
{
    TurbulentLength_ = turbulent_length_ratio_for_epsilon_inlet_ * CharacteristicLength_;
}
//=================================================================================================//
void kOmegaInflowTurbulentCondition::update(size_t index_i, Real dt)
{
    Real target_inflow_turbu_k = getTurbulentInflowK(pos_[index_i], vel_[index_i], turbu_k_[index_i]);
    turbu_k_[index_i] += relaxation_rate_ * (target_inflow_turbu_k - turbu_k_[index_i]);
    Real target_inflow_temp_turbu_E = getTurbulentInflowTemporaryEpsilon(pos_[index_i], turbu_k_[index_i], 0.0);

    //** Calculate inlet omega from k and epsilon */
    Real inflow_turbu_omega_temp = target_inflow_temp_turbu_E / std_kw_beta_star_ / target_inflow_turbu_k;
    //** Also has a temporary treatment because I don't know how to set frame transfer */
    Real target_inflow_turbu_omega = (pos_[index_i][0] < 0.0) ? inflow_turbu_omega_temp : turbu_omega_[index_i];
    turbu_omega_[index_i] += relaxation_rate_ * (target_inflow_turbu_omega - turbu_omega_[index_i]);
}
//=================================================================================================//
Real kOmegaInflowTurbulentCondition::getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k)
{
    Real u = velocity[0];
    Real temp_in_turbu_k = 1.5 * pow((turbulent_intensity_ * u), 2);
    Real turbu_k_original = turbu_k;
    if (type_turbu_inlet_ == 1)
    {
        Real channel_height = CharacteristicLength_; //** Temporarily treatment *

        //** Impose fully-developed K from PYTHON result */
        //** Calculate the distance to wall, Y. position here is the actual postion in x-y coordinate, no transformation*/
        Real Y = (position[1] < channel_height / 2.0) ? position[1] : channel_height - position[1];

        int polynomial_order = 8;
        int num_coefficient = polynomial_order + 1;
        //** Coefficient of the polynomial, 8th-order, from py21 dp=0.1 */
        Real coeff[] = {
            1.159981e-02, -4.662944e-02, 2.837400e-01,
            -1.193955e+00, 3.034851e+00, -4.766077e+00,
            4.529136e+00, -2.380854e+00, 5.307586e-01};
        Real polynomial_value = 0.0;
        for (int i = 0; i < num_coefficient; ++i)
        {
            polynomial_value += coeff[i] * std::pow(Y, i);
        }

        if (Y > channel_height / 2.0 || Y < 0.0)
        {
            std::cout << "position[1]=" << position[1] << std::endl;
            std::cout << "Y=" << Y << std::endl;
            std::cout << "polynomial_value=" << polynomial_value << std::endl;
            std::cout << "Stop" << std::endl;
            std::cout << "=================" << std::endl;
            std::cin.get();
        }

        temp_in_turbu_k = polynomial_value;
    }
    return (position[0] < 0.0) ? temp_in_turbu_k : turbu_k_original; //** Temporarily treatment *
}
//=================================================================================================//
Real kOmegaInflowTurbulentCondition::getTurbulentInflowTemporaryEpsilon(Vecd &position, Real &turbu_k, Real turbu_E)
{
    Real temp_in_turbu_E = C_mu_75_ * pow(turbu_k, 1.5) / TurbulentLength_;
    Real turbu_E_original = turbu_E;
    if (type_turbu_inlet_ == 1)
    {
        Real channel_height = CharacteristicLength_; //** Temporarily treatment *

        //** Impose fully-developed K from PYTHON result */
        //** Calculate the distance to wall, Y. position here is the actual postion in x-y coordinate, no transformation*/
        Real Y = (position[1] < channel_height / 2.0) ? position[1] : channel_height - position[1];

        int polynomial_order = 8;
        int num_coefficient = polynomial_order + 1;
        //** Coefficient of the polynomial, 8th-order, from py21 dp=0.1 */
        Real coeff[] = {
            1.428191e-02, -1.766636e-01, 1.153107e+00,
            -4.515606e+00, 1.103752e+01, -1.694146e+01,
            1.584534e+01, -8.241577e+00, 1.825421e+00};

        Real polynomial_value = 0.0;
        for (int i = 0; i < num_coefficient; ++i)
        {
            polynomial_value += coeff[i] * std::pow(Y, i);
        }
        temp_in_turbu_E = polynomial_value;
    }

    return (position[0] < 0.0) ? temp_in_turbu_E : turbu_E_original; //** Temporarily treatment *
}
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//