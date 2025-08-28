// #pragma once
#include "mph_method.hpp"

namespace SPH
{
//=================================================================================================//    
KernelOfMPH::KernelOfMPH(Real h, Real particle_spacing)
    : Kernel(h, 2.0, 2.0, "KernelOfMPH")
{
    Real Swg_2D = 1.0 / 2.0 * 1.0 / 3.0 * Pi / particle_spacing / particle_spacing;
    Real h_p = 2.0 * h;
    factor_W_2D_ = (1.0 / Swg_2D) * (1.0 / (h_p * h_p));

    factor_W_1D_ = 0.0;
    factor_W_3D_ = 0.0;

    setDerivativeParameters();
}
//=================================================================================================//
Real KernelOfMPH::W_1D(const Real q) const
{
    return 0.0;
}
//=================================================================================================//
Real KernelOfMPH::W_2D(const Real q) const
{
    return pow(1.0 - 0.5 * q, 2);
}
//=================================================================================================//
Real KernelOfMPH::W_3D(const Real q) const
{
    return 0.0;
}
//=================================================================================================//
Real KernelOfMPH::dW_1D(const Real q) const
{
    return 0.0;
}
//=================================================================================================//
Real KernelOfMPH::dW_2D(const Real q) const
{
    return (1.0 - 0.5 * q)*(-1.0);
}
//=================================================================================================//
Real KernelOfMPH::dW_3D(const Real q) const
{
    return 0.0;
}
//=================================================================================================//
Real KernelOfMPH::d2W_1D(const Real q) const
{
    return 0.0;
}
//=================================================================================================//
Real KernelOfMPH::d2W_2D(const Real q) const
{
    return 0.0;
}
//=================================================================================================//
Real KernelOfMPH::d2W_3D(const Real q) const
{
    return 0.0;
}
//=================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
CalculateVolumeStrain<Inner<>>::CalculateVolumeStrain(BaseInnerRelation& inner_relation)
    : CalculateVolumeStrain<Base, DataDelegateInner>(inner_relation), 
    number_density_Lattice_SPH_(sph_body_.getSPHAdaptation().LatticeNumberDensity()),
    smoothing_length_(sph_body_.getSPHAdaptation().ReferenceSmoothingLength()),
    particle_spacing_(sph_body_.getSPHAdaptation().ReferenceSpacing())
{
    Real Swg_2D = 1.0 / 2.0 * 1.0 / 3.0 * Pi / particle_spacing_ / particle_spacing_;
    Real h_p = 2.0 * smoothing_length_;
    Real W_ij_R0_2D = (1.0 / Swg_2D) * (1.0 / (h_p * h_p));
    number_density_Lattice_MPH_ = number_density_Lattice_SPH_ - W_ij_R0_2D;
    std::cout << "The particle number density is: " << number_density_Lattice_MPH_ << std::endl;
    std::cout << "Press anykey to continue " << std::endl;
    std::cin.get();
}
//=================================================================================================//
void CalculateVolumeStrain<Inner<>>::interaction(size_t index_i, Real dt)
{
    volume_strain_[index_i] = 0.0;

    Real sum_W_ij = 0.0;
    const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (index_i != index_j)
        {
            sum_W_ij += inner_neighborhood.W_ij_[n];
        }
    }
    volume_strain_[index_i] += sum_W_ij;
}
//=================================================================================================//
void CalculateVolumeStrain<Inner<>>::update(size_t index_i, Real dt)
{
    volume_strain_[index_i] -= number_density_Lattice_MPH_;
}
//=================================================================================================//
CalculateVolumeStrain<Contact<>>::CalculateVolumeStrain(BaseContactRelation& contact_relation)
    : CalculateVolumeStrain<Base, DataDelegateContact>(contact_relation) {}
//=================================================================================================//
void CalculateVolumeStrain<Contact<>>::interaction(size_t index_i, Real dt)
{
    Real sum_W_ij = 0.0;
    for (size_t k = 0; k < DataDelegateContact::contact_configuration_.size(); ++k)
    {
        Neighborhood& contact_neighborhood = (*DataDelegateContact::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (index_i != index_j)
            {
                sum_W_ij += contact_neighborhood.W_ij_[n];
            }
        }
    }
    volume_strain_[index_i] += sum_W_ij;
}
//=================================================================================================//
UpdateVelocity::UpdateVelocity(SPHBody& sph_body)
    : LocalDynamics(sph_body),
    vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
    pos_(particles_->getVariableDataByName<Vecd>("Position")) {}
//=================================================================================================//
void UpdateVelocity::update(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt;
}
//=================================================================================================//
ResetForce::ResetForce(SPHBody& sph_body)
    : LocalDynamics(sph_body),
    force_(particles_->registerStateVariableData<Vecd>("Force")),
    force_prior_(particles_->registerStateVariableData<Vecd>("ForcePrior")) {}
//=================================================================================================//
void ResetForce::update(size_t index_i, Real dt)
{
    force_[index_i] = Vecd::Zero();
    force_prior_[index_i] = Vecd::Zero();
}
//=================================================================================================//
}// namespace fluid_dynamics
//=================================================================================================//
}// namespace SPH