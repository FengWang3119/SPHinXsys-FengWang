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
      turbu_k_(*(this->particles_->template registerSharedVariable<Real>("TurbulenceKineticEnergy"))),
      turbu_epsilon_(*(this->particles_->template registerSharedVariable<Real>("TurbulentDissipation"))),
      turbu_mu_(*(this->particles_->template registerSharedVariable<Real>("TurbulentViscosity"))),
      turbu_strain_rate_(*(this->particles_->template registerSharedVariable<Matd>("TurbulentStrainRate"))),
      turbu_strain_rate_magnitude_(*(this->particles_->template registerSharedVariable<Real>("TurbulentStrainRateMagnitude"))),
      turbu_omega_(*(this->particles_->template registerSharedVariable<Real>("TurbulentSpecificDissipation"))),
      mu_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      particle_spacing_min_(base_relation.getSPHBody().sph_adaptation_->MinimumSpacing()),
      rho_(*this->particles_->template getVariableDataByName<Real>("Density")),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(*this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      dimension_(2)
{
}
//A temporarily treatment for dimension
//=================================================================================================//

} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//