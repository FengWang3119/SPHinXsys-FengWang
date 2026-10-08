//#pragma once
#include "zeroth-order_residue.hpp"
namespace SPH
{
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
    //if (pos_i_x > 0.0 && pos_i_x < 15.0)
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
    //if (pos_i_x > 0.0 && pos_i_x < 15.0)
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
CalculateAverageRiemannDissipation::
CalculateAverageRiemannDissipation(SPHBody& sph_body, Real channel_length)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
    dissipation_riemann_(particles_->registerStateVariable<Real>("RiemannDissipation")),
    pos_(particles_->getVariableDataByName<Vecd>("Position")),
    channel_length_(channel_length),
    particle_spacing_(sph_body_.sph_adaptation_->ReferenceSpacing()) {}
//=================================================================================================//
Real CalculateAverageRiemannDissipation::reduce(size_t index_i, Real dt)
{
    Real pos_i_x = pos_[index_i][xAxis];
    if (pos_i_x > 0.0 && pos_i_x < (channel_length_ - 10.0 * particle_spacing_))
    //if (pos_i_x > 0.0 && pos_i_x < 15.0)
    {
        return std::abs(dissipation_riemann_[index_i]);
    }
    return 0.0;
}
//=================================================================================================//
Real CalculateAverageRiemannDissipation::outputResult(Real reduced_value)
{
    return reduced_value / (Real(num_particle_in_domain_));
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//