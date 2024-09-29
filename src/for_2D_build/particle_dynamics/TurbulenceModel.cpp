#ifndef TURBULENCEMODEL_CPP
#define TURBULENCEMODEL_CPP
#include "TurbulenceModel.h"
#include "sphinxsys.h"
namespace SPH
{  
    namespace fluid_dynamics
    { 
    //=================================================================================================//
        BaseTurbulence::BaseTurbulence(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator)
            : BaseIntegration<DataDelegateInner>(inner_relation), ghost_creator_(ghost_creator),
            mom_(this->particles_->template registerStateVariable<Vecd>("Momentum")),
          dmom_dt_(this->particles_->template registerStateVariable<Vecd>("MomentumChangeRate")),
          dmass_dt_(this->particles_->template registerStateVariable<Real>("MassChangeRate")),
          K_prod_p_(this->particles_->template registerStateVariable<Real>("TKEProductionInWallAdjCell")),
          K_prod_(this->particles_->template registerStateVariable<Real>("TKEProduction")),
          Eps_p_(this->particles_->template registerStateVariable<Real>("DissipationRateInWallAdjCell")),
          Eps_sum_(this->particles_->template registerStateVariable<Real>("DissipationSum")),
          K_adv_(this->particles_->template registerStateVariable<Real>("TKEAdvection")),
          K_lap_(this->particles_->template registerStateVariable<Real>("TKELaplacian")),
          Eps_adv_(this->particles_->template registerStateVariable<Real>("DissipationAdvection")),
          Tau_wall_(this->particles_->template registerStateVariable<Real>("WallShearStress")),
          Eps_lap_(this->particles_->template registerStateVariable<Real>("DissipationLaplacian")),
          Eps_prodscalar_(this->particles_->template registerStateVariable<Real>("DissipationProdscalar")),
          Eps_scalar_(this->particles_->template registerStateVariable<Real>("DissipationScalar")),
          Cmu_(0.09), sigmak_(1.0), sigmaeps_(1.3), C1eps_(1.44), C2eps_(1.92),
          K_(this->particles_->template registerStateVariable<Real>("TKE")),
          Eps_(this->particles_->template registerStateVariable<Real>("Dissipation")),
          mu_t_(this->particles_->template registerStateVariable<Real>("TurblunetViscosity"))
            {}
        //=================================================================================================//
           
            WallAdjacentCells::WallAdjacentCells(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
                : BaseTurbulence(inner_relation, ghost_creator), yp_(this->particles_->template registerStateVariable<Real>("WallNormalDistance")),
                  walladjacentcellflag_(this->particles_->template registerStateVariable<Real>("FlagForWallAdjacentCells")),
                  wallnormal_(this->particles_->template registerStateVariable<Vecd>("WallNormal")), ymax_(0.0), 
                  bounds_(inner_relation.getSPHBody())
            {
            walladjacentcellyp();
            }
        //=================================================================================================//
        void WallAdjacentCells::walladjacentcellyp()
        {
            walladjacentindex_.resize(particles_->ParticlesBound());
            wallghostindex_.resize(particles_->ParticlesBound());
            walleij_.resize(particles_->ParticlesBound());
            Real boundary_type = 3;
                
            if (!ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
            {
                for (size_t ghost_number = 0; ghost_number != ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
                {

                    size_t ghost_index = ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                    size_t index_real = ghost_creator_.each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                    walleij_[index_real] = ghost_creator_.each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];
                    walladjacentindex_[index_real] = index_real;
                    wallghostindex_[index_real] = ghost_index;
                    if ((pos_[index_real] - pos_[ghost_index]).dot(walleij_[index_real]) > ymax_)
                    {
                        ymax_ = (pos_[index_real] - pos_[ghost_index]).dot(walleij_[index_real]);

                    }
                    
                }
            }
        }
        //=================================================================================================//
        
        void WallAdjacentCells::update(size_t index_i, Real dt)
        {
         Vecd wallnormal = Vecd::Zero();
            Vecd lower_wall = {pos_[index_i][0], 0.0};
            Vecd upper_wall = {pos_[index_i][0], 2.0};
            Vecd lower_wall_normal = {0.0, 1.0};
            Vecd upper_wall_normal = {0.0, -1.0};
            /*
            BoundingBox bounds = bounds_.getSPHSystemBounds();
            Real channelheight = bounds.second_[1];
            Real halfwidth = 0.5 * channelheight;
            */ 

            bool lower_wall_condition = ((pos_[index_i] - lower_wall).dot(lower_wall_normal) <= 1.0 * ymax_);
            bool upper_wall_condition = ((pos_[index_i] - upper_wall).dot(upper_wall_normal) <= 1.0 * ymax_);

            if (lower_wall_condition)
            {
                yp_[index_i] = (pos_[index_i] - lower_wall).dot(lower_wall_normal);
                walladjacentcellflag_[index_i] = 1.0;
                wallnormal_[index_i] = lower_wall_normal;
            }
            else if (upper_wall_condition)
            {
                yp_[index_i] = (pos_[index_i] - upper_wall).dot(upper_wall_normal);
                walladjacentcellflag_[index_i] = 1.0;    
                wallnormal_[index_i] = upper_wall_normal;
            }
        }
        
        
        //=================================================================================================//

        KEpsilonStd1stHalf::KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator) 
            : StdWallFunctionFVM(inner_relation, ghost_creator),
              dK_dt_(this->particles_->template registerStateVariable<Real>("TKEChangeRate")),
              walladjacentcellflag_(this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells")),
              strain_rate_(this->particles_->template registerStateVariable<Real>("StrainRate")),
              dudx_(this->particles_->template registerStateVariable<Real>("dudx")),
              dudy_(this->particles_->template registerStateVariable<Real>("dudy")),
              dvdx_(this->particles_->template registerStateVariable<Real>("dvdx")),
              dvdy_(this->particles_->template registerStateVariable<Real>("dvdy")),
              vel_gradient_mat_(this->particles_->template getVariableDataByName<Matd>("VelocityGradient"))
        {}
        //=================================================================================================//
        void KEpsilonStd1stHalf::interaction(size_t index_i, Real dt)
        {
            
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            Matd K_prod = Matd::Zero(), K_prod_iso = Matd::Zero(), K_prod_total = Matd::Zero();
            Matd vel_matrix = Matd::Zero(), vel_gradient_mat = Matd::Zero();
            Matd strain_tensor = Matd::Zero(), strain_rate_modulus = Matd::Zero();
            Real Kprodtot = 0.0;
            K_prod_[index_i] = 0.0, K_adv_[index_i] = 0.0, K_lap_[index_i] = 0.0, strain_rate_[index_i] = 0.0;
            Eps_sum_[index_i] = 0.0, dudx_[index_i] = 0.0, dudy_[index_i] = 0.0, dvdx_[index_i] = 0.0, dvdy_[index_i] = 0.0;
            vel_gradient_mat_[index_i] = Matd::Zero();
            Real mu_t_upperlimit = 1e4 * fluid_.ReferenceViscosity();
            Real mu_t_lowerlimit = 1e-3 * fluid_.ReferenceViscosity();
            Real mu_t = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
            mu_t_[index_i] = std::max(std::min(mu_t_upperlimit, mu_t), mu_t_lowerlimit);

            if (walladjacentcellflag_[index_i] == 1.0)
            {
                nearwallquantities(index_i);
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (K_[index_i] - K_[index_j])) * ((vel_[index_i]).dot(e_ij));
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigmak_) * (K_[index_i] - K_[index_j]) / (r_ij));    
                }
                strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
                strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
                strain_rate_[index_i] = std::sqrt(strain_rate_modulus.sum());
                
                K_prod_[index_i] = K_prod_p_[index_i];
                Eps_[index_i] = Eps_p_[index_i];
                
                dudx_[index_i] = vel_gradient_mat_[index_i](0, 0);
                dudy_[index_i] = vel_gradient_mat_[index_i](0, 1);
                dvdx_[index_i] = vel_gradient_mat_[index_i](1, 0);
                dvdy_[index_i] = vel_gradient_mat_[index_i](1, 1);

                dK_dt_[index_i] = K_adv_[index_i] + K_prod_[index_i] - rho_[index_i] * Eps_[index_i] + K_lap_[index_i]; 
               
            }
            else
            {
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (K_[index_i] - K_[index_j])) * ((vel_[index_i]).dot(e_ij));
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigmak_) * (K_[index_i] - K_[index_j]) / (r_ij));
                    
                    vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                    vel_gradient_mat_[index_i] += dW_ij * Vol_[index_j] * vel_matrix;

                }
                strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
                strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
                strain_rate_[index_i] = std::sqrt(strain_rate_modulus.sum());
                
                K_prod = (mu_t_[index_i] * strain_rate_modulus);
                K_prod_[index_i] = K_prod.sum();

                dudx_[index_i] = vel_gradient_mat_[index_i](0, 0);
                dudy_[index_i] = vel_gradient_mat_[index_i](0, 1);
                dvdx_[index_i] = vel_gradient_mat_[index_i](1, 0);
                dvdy_[index_i] = vel_gradient_mat_[index_i](1, 1);

                dK_dt_[index_i] = K_adv_[index_i] + K_prod_[index_i] - rho_[index_i] * Eps_[index_i] + K_lap_[index_i];
            }
        }
        //=================================================================================================//
        void KEpsilonStd1stHalf::update(size_t index_i, Real dt)
        {
            K_[index_i] += (dK_dt_[index_i] / rho_[index_i]) * dt;
            if (K_[index_i] < 0.0)
            {
                K_[index_i] = 1e-7;
            }
        }
        //=================================================================================================//
        KEpsilonStd2ndHalf::KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator) 
            : BaseTurbulence(inner_relation, ghost_creator),
            dEps_dt_(this->particles_->template registerStateVariable<Real>("DissipationChangeRate")),
            walladjacentcellflag_(this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells"))
        {}
        //=================================================================================================//
        
        void KEpsilonStd2ndHalf::interaction(size_t index_i, Real dt)
        {
            Real Eps_changerate = 0.0;
            Eps_adv_[index_i] = 0.0, Eps_lap_[index_i] = 0.0, Eps_prodscalar_[index_i] = 0.0, Eps_scalar_[index_i] = 0.0;
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            if (walladjacentcellflag_[index_i] != 1)
            {
                Real mu_t_upperlimit = 1e4 * fluid_.ReferenceViscosity();
                Real mu_t_lowerlimit = 1e-3 * fluid_.ReferenceViscosity();
                Real mu_t = rho_[index_i] * Cmu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
                mu_t_[index_i] = std::max(std::min(mu_t_upperlimit, mu_t), mu_t_lowerlimit);
            }   
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ij = inner_neighborhood.dW_ij_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];
                Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

                Eps_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (Eps_[index_i] - Eps_[index_j])) * ((vel_[index_i]).dot(e_ij));
                Eps_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * (fluid_.ReferenceViscosity() + mu_t_avg / sigmaeps_) * ((Eps_[index_i] - Eps_[index_j]) / (r_ij));

                Eps_changerate = Eps_adv_[index_i] + Eps_lap_[index_i];
            }
                Eps_prodscalar_[index_i] = C1eps_ * (Eps_[index_i] / (K_[index_i])) * K_prod_[index_i];
                Eps_scalar_[index_i] = -C2eps_ * rho_[index_i] * (Eps_[index_i] * Eps_[index_i]) / (K_[index_i]);
                dEps_dt_[index_i] = Eps_changerate + Eps_prodscalar_[index_i] + Eps_scalar_[index_i];   
        }
        //=================================================================================================//
        void KEpsilonStd2ndHalf::update(size_t index_i, Real dt)
        {
            if (walladjacentcellflag_[index_i] != 1)
            {
                Eps_[index_i] += (dEps_dt_[index_i] / rho_[index_i]) * dt;
            }
            if (Eps_[index_i] < 0.0)
            {
                Eps_[index_i] = 1e-7;
            }
        }
        //=================================================================================================//
         StdWallFunctionFVM::StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
            : WallAdjacentCells(inner_relation, ghost_creator), ystar_(this->particles_->template registerStateVariable<Real>("Ystar")),
              vel_gradient_mat_(this->particles_->template registerStateVariable< Matd > ("VelocityGradient")),
              vonkar_(0.4187), E_(9.793)
           {}
        //=================================================================================================//
        void StdWallFunctionFVM::nearwallquantities(size_t index_i)
        {
            ystar_[index_i] = (rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * yp_[index_i]) / (fluid_.ReferenceViscosity());
            Real u_star, mu_eff_wall, friction_velocity;
            Vecd veltangential = (vel_[index_i] - wallnormal_[index_i].dot(vel_[index_i]) * (wallnormal_[index_i]));

            if (ystar_[index_i] >= 11.225)
            {
                u_star = (1.0 / vonkar_) * std::log(E_ * ystar_[index_i]);
                mu_t_[index_i] = fluid_.ReferenceViscosity() * ((ystar_[index_i]) / (1 / vonkar_ * std::log(E_ * ystar_[index_i])) - 1.0);
                
                Tau_wall_[index_i] = (veltangential.norm() * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                vel_gradient_mat_[index_i] = Matd::Zero();
                vel_gradient_mat_[index_i](0, 1) = Tau_wall_[index_i] / (rho_[index_i] * pow(Cmu_, 0.25) * pow(K_[index_i], 0.5) * vonkar_ * yp_[index_i]);
                //vel_gradient_mat_[index_i](0, 1) = veltangential.norm() / (yp_[index_i] * log(E_ * ystar_[index_i]));

                K_prod_p_[index_i] = std::pow(Tau_wall_[index_i], 2.0) / (vonkar_ * rho_[index_i] * std::pow(Cmu_, 0.25) * std::pow(K_[index_i], 0.5) * yp_[index_i]);
                Eps_p_[index_i] = (std::pow(Cmu_, 3.0 / 4.0) * std::pow(K_[index_i], 1.5)) / (vonkar_ * yp_[index_i]);
               
            }
            else if (ystar_[index_i] < 11.225)
            {
                u_star = ystar_[index_i];
                Tau_wall_[index_i] = fluid_.ReferenceViscosity() * veltangential.norm() / yp_[index_i];
                vel_gradient_mat_[index_i] = Matd::Zero();
                vel_gradient_mat_[index_i](0, 1) = Tau_wall_[index_i] / fluid_.ReferenceViscosity();
                K_prod_p_[index_i] = 0.0;
                Eps_p_[index_i] = (K_[index_i] * 2.0 * fluid_.ReferenceViscosity()) / (rho_[index_i] * yp_[index_i] * yp_[index_i]);
            }  
        }
        //=================================================================================================// 
    }// namespace fluid_dynamics

}// namespace SPH
#endif // TURBULENCEMODEL_CPP