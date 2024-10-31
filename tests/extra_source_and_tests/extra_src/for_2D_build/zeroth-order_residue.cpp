//#pragma once
#include "zeroth-order_residue.h"
namespace SPH
{
//=================================================================================================//
GetLimiterOfTransportVelocityCorrection::
    GetLimiterOfTransportVelocityCorrection(SPHBody &sph_body, Real slope)
    : LocalDynamics(sph_body),
      h_ref_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      zero_gradient_residue_(particles_->getVariableDataByName<Vecd>("ZeroGradientResidue")),
      slope_(slope),
      limiter_tvc_(particles_->registerStateVariable<Real>("LimiterOfTVC"))
{
    particles_->addVariableToWrite<Real>("LimiterOfTVC");
}
//=================================================================================================//
void GetLimiterOfTransportVelocityCorrection::update(size_t index_i, Real dt)
{
    Real squared_norm = zero_gradient_residue_[index_i].squaredNorm();
    limiter_tvc_[index_i] = SMIN(slope_ * squared_norm * h_ref_ * h_ref_, Real(1));
}
//=================================================================================================//
GetPressureGradientResidue::
    GetPressureGradientResidue(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      zero_gradient_residue_(particles_->getVariableDataByName<Vecd>("ZeroGradientResidue")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      pressure_gradient_residue_(particles_->registerStateVariable<Vecd>("PressureGradientResidue"))
{
    particles_->addVariableToWrite<Vecd>("PressureGradientResidue");
}
//=================================================================================================//
void GetPressureGradientResidue::update(size_t index_i, Real dt)
{
    pressure_gradient_residue_[index_i] = zero_gradient_residue_[index_i] * p_[index_i];
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//