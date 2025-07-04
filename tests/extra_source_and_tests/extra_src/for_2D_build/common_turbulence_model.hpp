#pragma once

#include "common_turbulence_model.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TKEnergyForce<Base, DataDelegationType>::
    TKEnergyForce(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      force_(this->particles_->template registerStateVariable<Vecd>("Force")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      indicator_(this->particles_->template getVariableDataByName<int>("Indicator")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      turbu_k_(this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      test_k_grad_rslt_(this->particles_->template registerStateVariable<Vecd>("TkeGradResult")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
GetVelocityGradient<DataDelegationType>::
    GetVelocityGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      is_near_wall_P1_(this->particles_->template getVariableDataByName<int>("IsNearWallP1")),
      is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      velocity_gradient_(this->particles_->template registerStateVariable<Matd>("TurbulentVelocityGradient")),
      velocity_gradient_wall(this->particles_->template registerStateVariable<Matd>("Velocity_Gradient_Wall")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TurbuViscousForce<DataDelegationType>::TurbuViscousForce(BaseRelationType &base_relation)
    : ViscousForce<DataDelegationType>(base_relation),
      turbu_k_(this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_mu_(this->particles_->template getVariableDataByName<Real>("TurbulentViscosity")),
      wall_Y_plus_(this->particles_->template getVariableDataByName<Real>("WallYplus")),
      wall_Y_star_(this->particles_->template getVariableDataByName<Real>("WallYstar")),
      velo_friction_(this->particles_->template getVariableDataByName<Vecd>("FrictionVelocity")),
      y_p_(this->particles_->template getVariableDataByName<Real>("Y_P")),
      is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      viscosity_(DynamicCast<Viscosity>(this, this->particles_->getBaseMaterial())),
      molecular_viscosity_(viscosity_.ReferenceViscosity()),
      c0_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceSoundSpeed()) {}

//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TurbulentLinearGradientCorrectionMatrix<DataDelegationType>::
    TurbulentLinearGradientCorrectionMatrix(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      turbu_B_(this->particles_->template registerStateVariable<Matd>(
          "TurbulentLinearGradientCorrectionMatrix", IdentityMatrix<Matd>::value)),
      B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
{
    this->particles_->template addVariableToWrite<Matd>("TurbulentLinearGradientCorrectionMatrix");
    this->particles_->template addVariableToSort<Matd>("TurbulentLinearGradientCorrectionMatrix");
    this->particles_->template addVariableToWrite<Matd>("LinearGradientCorrectionMatrix");
    this->particles_->template addVariableToSort<Matd>("LinearGradientCorrectionMatrix");
}

//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//