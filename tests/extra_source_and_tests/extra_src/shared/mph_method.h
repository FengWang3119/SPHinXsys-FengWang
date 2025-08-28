#ifndef MPH_METHOD_H
#define MPH_METHOD_H

#include "sphinxsys.h"

namespace SPH
{
int test_ccc = 1;

class KernelOfMPH : public Kernel
{
public:
    KernelOfMPH(Real h, Real particle_spacing);

    /**
     * Calculates the kernel value for
     * the given distance of two particles
     */
    virtual Real W_1D(const Real q) const override;
    virtual Real W_2D(const Real q) const override;
    virtual Real W_3D(const Real q) const override;

    virtual Real dW_1D(const Real q) const override;
    virtual Real dW_2D(const Real q) const override;
    virtual Real dW_3D(const Real q) const override;

    virtual Real d2W_1D(const Real q) const override;
    virtual Real d2W_2D(const Real q) const override;
    virtual Real d2W_3D(const Real q) const override;
};
namespace fluid_dynamics
{
//=================================================================================================//
template <typename... InteractionTypes>
class CalculateVolumeStrain;

template <class DataDelegationType>
class CalculateVolumeStrain<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
public:
    template <class BaseRelationType>
    explicit CalculateVolumeStrain(BaseRelationType& base_relation);
    virtual ~CalculateVolumeStrain() {};

protected:
    Real *volume_strain_;
};

template <>
class CalculateVolumeStrain<Inner<>> : public CalculateVolumeStrain<Base, DataDelegateInner>
{
public:
    explicit CalculateVolumeStrain(BaseInnerRelation& inner_relation);
    virtual ~CalculateVolumeStrain() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

protected:
    Real number_density_Lattice_SPH_;
    Real number_density_Lattice_MPH_;
    Real smoothing_length_;
    Real particle_spacing_;
};

template <>
class CalculateVolumeStrain<Contact<>> : public CalculateVolumeStrain<Base, DataDelegateContact>
{
public:
    explicit CalculateVolumeStrain(BaseContactRelation& contact_relation);
    virtual ~CalculateVolumeStrain() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:

};

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseCalculateVolumeStrainComplex = ComplexInteraction<CalculateVolumeStrain<InnerInteractionType, ContactInteractionTypes...>>;

using CalculateVolumeStrainComplex = BaseCalculateVolumeStrainComplex<Inner<>, Contact<>>;
//=================================================================================================//
class UpdateVelocity : public LocalDynamics
{
public:
    explicit UpdateVelocity(SPHBody& sph_body);
    virtual ~UpdateVelocity() {};

    void update(size_t index_i, Real dt = 0.0);

protected:
    Vecd* vel_;
    Vecd* pos_;
};
//=================================================================================================//
class ResetForce : public LocalDynamics
{
public:
    explicit ResetForce(SPHBody& sph_body);
    virtual ~ResetForce() {};

    void update(size_t index_i, Real dt = 0.0);

protected:
    Vecd *force_, *force_prior_;
};
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
#endif // MPH_METHOD_H