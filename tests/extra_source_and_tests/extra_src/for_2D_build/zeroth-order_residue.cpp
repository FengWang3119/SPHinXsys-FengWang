//#pragma once
#include "zeroth-order_residue.hpp"
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
GetPressureGradientResidue_RKGC<Inner<>>::GetPressureGradientResidue_RKGC(BaseInnerRelation &inner_relation)
    : GetPressureGradientResidue_RKGC<Base, DataDelegateInner>(inner_relation) {}
//=================================================================================================//
void GetPressureGradientResidue_RKGC<Inner<>>::interaction(size_t index_i, Real dt)
{
    pressure_gradient_residue_RKGC_[index_i] = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
        pressure_gradient_residue_RKGC_[index_i] -= (B_[index_j] * p_[index_i] + B_[index_i] * p_[index_j]) * nablaW_ijV_j;
    }
}
//=================================================================================================//
GetPressureGradientResidue_RKGC<Contact<>>::GetPressureGradientResidue_RKGC(BaseContactRelation &contact_relation)
    : GetPressureGradientResidue_RKGC<Base, DataDelegateContact>(contact_relation) {}
//=================================================================================================//
void GetPressureGradientResidue_RKGC<Contact<>>::interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * this->Vol_[index_j] * contact_neighborhood.e_ij_[n];
            pressure_gradient_residue_RKGC_[index_i] -= (p_[index_i] + p_[index_i]) * B_[index_i] * nablaW_ijV_j;
        }
    }
}
//=================================================================================================//
NonDimensionalisePressure::
    NonDimensionalisePressure(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      p_dimensionless_(particles_->registerStateVariable<Real>("PressureDimensionless"))
{
    particles_->addVariableToWrite<Real>("PressureDimensionless");
}
//=================================================================================================//
void NonDimensionalisePressure::update(size_t index_i, Real dt)
{
    p_dimensionless_[index_i] = p_[index_i] / rho_[index_i];
}
//=================================================================================================//
CalculateAverageKGS::
CalculateAverageKGS(SPHBody& sph_body, Real channel_length)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
    zero_gradient_residue_(particles_->getVariableDataByName<Vecd>("ZeroGradientResidue")),
    pos_(particles_->getVariableDataByName<Vecd>("Position")),
    channel_length_(channel_length),
    particle_spacing_(sph_body_.sph_adaptation_->ReferenceSpacing()),
    total_real_particles_(particles_->TotalRealParticles()) {}
//=================================================================================================//
Real CalculateAverageKGS::reduce(size_t index_i, Real dt)
{
    Real pos_i_x = pos_[index_i][xAxis];
    if (pos_i_x > 0.0 && pos_i_x < (channel_length_ - 10.0 * particle_spacing_))
    {
        return zero_gradient_residue_[index_i].norm();
    }
    return 0.0;
}
//=================================================================================================//
Real CalculateAverageKGS::outputResult(Real reduced_value)
{
    return reduced_value / (Real(num_particle_in_domain_));
}
//=================================================================================================//
CalculateParticleInDomain::
CalculateParticleInDomain(SPHBody& sph_body, Real channel_length)
    : LocalDynamicsReduce<ReduceSum<size_t>>(sph_body),
    pos_(particles_->getVariableDataByName<Vecd>("Position")),
    channel_length_(channel_length),
    particle_spacing_(sph_body_.sph_adaptation_->ReferenceSpacing()) {}
//=================================================================================================//
size_t CalculateParticleInDomain::reduce(size_t index_i, Real dt)
{
    Real pos_i_x = pos_[index_i][xAxis];
    if (pos_i_x > 0.0 && pos_i_x < (channel_length_ - 10.0 * particle_spacing_))
    {
        return 1;
    }
    return 0;
}
//=================================================================================================//
size_t CalculateParticleInDomain::outputResult(size_t reduced_value)
{
    return reduced_value;
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//