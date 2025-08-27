#pragma once

#include "mph_method.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
CalculateVolumeStrain<Base, DataDelegationType>::
    CalculateVolumeStrain(BaseRelationType& base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
    volume_strain_(this->particles_->template registerStateVariableData<Real>("VolumeStrain")) {}
//=================================================================================================//
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//