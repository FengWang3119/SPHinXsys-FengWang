#pragma once

#include "k-omega_turbulent_model.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
kOmega_BaseTurbulentModel<Base, DataDelegationType>::kOmega_BaseTurbulentModel(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      turbu_strain_rate_(this->particles_->template registerStateVariable<Matd>("TurbulentStrainRate")),
      turbu_strain_rate_magnitude_(this->particles_->template registerStateVariable<Real>("TurbulentStrainRateMagnitude")),
      turbu_strain_rate_traceless_magnitude_(this->particles_->template registerStateVariable<Real>("TurbulentStrainRateTracelessMagnitude")),
      mu_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      particle_spacing_min_(base_relation.getSPHBody().sph_adaptation_->MinimumSpacing()),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      dimension_(2) {}
//A temporarily treatment for dimension
//=================================================================================================//

} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//