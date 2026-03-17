#pragma once

#include "udf_P_refinement.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
namespace udf
{
//=================================================================================================//
    template <class DataDelegationType>
    template <class BaseRelationType>
    P_refinement_GetVelocityGradientInnerOnlyP<DataDelegationType>::
        P_refinement_GetVelocityGradientInnerOnlyP(BaseRelationType& base_relation)
        : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
        velocity_gradient_inner_only_P_(this->particles_->template registerStateVariableData<Matd>("VelocityGradientInnerOnlyP")),
        Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
        vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
        is_near_wall_P1_(this->particles_->template getVariableDataByName<int>("IsNearWallP1")),
        is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")) {}
//=================================================================================================//
} // namespace udf
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//