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
    volume_strain_(this->particles_->template registerStateVariableData<Real>("VolumeStrain")),
    sum_weight_(this->particles_->template registerStateVariableData<Real>("WeightSummation")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
CalculateVelocityDivergence<Base, DataDelegationType>::
CalculateVelocityDivergence(BaseRelationType& base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
    vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
    velocity_divergence_(this->particles_->template registerStateVariableData<Real>("VelocityDivergence")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
CalculatePressureGradientForce<Base, DataDelegationType>::
CalculatePressureGradientForce(BaseRelationType& base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
    force_(this->particles_->template getVariableDataByName<Vecd>("Force")),
    p_(this->particles_->template getVariableDataByName<Real>("Pressure")),
    Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//