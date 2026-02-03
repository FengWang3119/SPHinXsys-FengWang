/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	k-epsilon_turbulent_model.h
 * @brief 	
 * @details     
 * @author Xiangyu Hu
 */

#ifndef UDF_K_EPSILON_TURBULENT_MODEL_H
#define UDF_K_EPSILON_TURBULENT_MODEL_H

#include "udf_common_turbulence_model.h"
#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
namespace udf
{
class kEpsilon_TurbulentClosureCoefficient
{
  public:
    explicit kEpsilon_TurbulentClosureCoefficient();
    virtual ~kEpsilon_TurbulentClosureCoefficient(){};

  protected:
    Real Karman_;
    Real turbu_const_E_;
    Real inv_turbu_E_;
    Real C_mu_, C_mu_25_, C_mu_75_;
    Real turbulent_intensity_;

    //** Closure coefficients for K *
    Real sigma_k_;

    //** Closure coefficients for Epsilon *
    Real C_l_, C_2_;
    Real sigma_E_;
    Real turbulent_length_ratio_for_epsilon_inlet_;
};
//=================================================================================================//
template <typename... T>
class kEpsilon_BaseTurbulentModel;

template <class DataDelegationType>
class kEpsilon_BaseTurbulentModel<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType, public kEpsilon_TurbulentClosureCoefficient
{
  public:
    template <class BaseRelationType>
    explicit kEpsilon_BaseTurbulentModel(BaseRelationType &base_relation);
    virtual ~kEpsilon_BaseTurbulentModel(){};

  protected:
    Matd *turbu_strain_rate_; //** temporary naming to distinguish the regular strain rate *
    Viscosity &viscosity_;
    Real mu_, smoothing_length_, particle_spacing_min_;
    Real *rho_, *Vol_;
    Vecd *vel_;
    int dimension_;
};
//=================================================================================================//
class kEpsilon_kTransportEquationInner : public kEpsilon_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kEpsilon_kTransportEquationInner(BaseInnerRelation &inner_relation, const StdVec<Real> &initial_values, int is_extr_visc_dissipa, bool is_STL);
    virtual ~kEpsilon_kTransportEquationInner(){};

    //inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *dk_dt_;
    Real *dk_dt_without_dissipation_;
    Real *k_production_;

    int *is_near_wall_P1_;
    Matd *velocity_gradient_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real *turbu_mu_;
    Matd *turbu_strain_rate_;
    int *is_extra_viscous_dissipation_;
    bool is_STL_;

    int *turbu_indicator_;
    Real *k_diffusion_;
};
//=================================================================================================//
class kEpsilon_TKE_Diffusion : public kEpsilon_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kEpsilon_TKE_Diffusion(BaseInnerRelation &inner_relation);
    virtual ~kEpsilon_TKE_Diffusion(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *turbu_k_;
    Real *turbu_mu_;
    Real *k_diffusion_;
};
//=================================================================================================//
class kEpsilon_epsilonTransportEquationInner : public kEpsilon_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kEpsilon_epsilonTransportEquationInner(BaseInnerRelation &inner_relation, bool is_STL);
    virtual ~kEpsilon_epsilonTransportEquationInner(){};

    //inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *depsilon_dt_;
    Real *depsilon_dt_without_dissipation_;
    Real *ep_production;
    Real *ep_dissipation_;
    Real *ep_diffusion_;

    Real *turbu_mu_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real *k_production_;
    int *is_near_wall_P1_;
    bool is_STL_;
};
//=================================================================================================//
class kEpsilon_TDR_Diffusion : public kEpsilon_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kEpsilon_TDR_Diffusion(BaseInnerRelation &inner_relation);
    virtual ~kEpsilon_TDR_Diffusion(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *ep_diffusion_;
    Real *turbu_mu_;
    Real *turbu_epsilon_;
};
//=================================================================================================//
class kEpsilon_TurbulentEddyViscosity : public LocalDynamics, public kEpsilon_TurbulentClosureCoefficient
{
  public:
    explicit kEpsilon_TurbulentEddyViscosity(SPHBody &sph_body);
    virtual ~kEpsilon_TurbulentEddyViscosity(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *rho_;
    Real *turbu_mu_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real *wall_Y_plus_, *wall_Y_star_;
    Viscosity &viscosity_;
    Real mu_;
};
//=================================================================================================//
/**
	* @class   kEpsilon_InflowTurbulentCondition
	* @brief   Inflow boundary condition which imposes directly to a given velocity profile.
	*          TargetVelocity gives the velocity profile along the inflow direction,
	*          i.e. x direction in local frame.
	*/
class kEpsilon_InflowTurbulentCondition : public BaseFlowBoundaryCondition, public kEpsilon_TurbulentClosureCoefficient
{
  public:
    explicit kEpsilon_InflowTurbulentCondition(BodyPartByCell &body_part,
                                               Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet);
    virtual ~kEpsilon_InflowTurbulentCondition(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int type_turbu_inlet_;
    Real relaxation_rate_;
    Real CharacteristicLength_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real TurbulentLength_;

    virtual Real getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k);
    virtual Real getTurbulentInflowE(Vecd &position, Real &turbu_k, Real &turbu_E);
};
//=================================================================================================//
class kEpsilon_StandardWallFunctionCorrection : public LocalDynamics, public DataDelegateContact, public WallFunction
{
  public:
    kEpsilon_StandardWallFunctionCorrection(BaseInnerRelation &inner_relation,
                                            BaseContactRelation &contact_relation);
    virtual ~kEpsilon_StandardWallFunctionCorrection(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *wall_Y_plus_, *wall_Y_star_;
    Real *velo_tan_;
    Vecd *velo_friction_;

    Real *y_p_;
    Vecd *vel_;
    Real *rho_;
    Viscosity &viscosity_;
    Real molecular_viscosity_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real *turbu_mu_;
    int *is_near_wall_P1_;
    int *is_near_wall_P2_;
    Matd *velocity_gradient_;
    Real *k_production_;
    Real *distance_to_dummy_interface_;
    Real *distance_to_dummy_interface_up_average_;
    int *index_nearest;
    Vecd *e_nearest_tau_;
    Vecd *e_nearest_normal_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_n_;
    Real *physical_time_;
};
//=================================================================================================//
}
} // namespace fluid_dynamics
} // namespace SPH
#endif // K_EPSILON_TURBULENT_MODEL_H