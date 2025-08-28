#pragma once

#include "mph_method.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
CalculateVolumeStrain<Base, DataDelegationType>::
    CalculateVolumeStrain(BaseRelationType& base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
    volume_strain_(this->particles_->template registerStateVariableData<Real>("VolumeStrain")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
CalculateVelocityDivergence<Base, DataDelegationType>::
CalculateVelocityDivergence(BaseRelationType& base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
    vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
    velocity_divergence_(this->particles_->template registerStateVariableData<Real>("VelocityDivergence")) {}
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//