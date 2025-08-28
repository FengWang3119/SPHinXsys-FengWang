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
template <typename... InteractionTypes>
class CalculateVelocityDivergence;

template <class DataDelegationType>
class CalculateVelocityDivergence<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
public:
    template <class BaseRelationType>
    explicit CalculateVelocityDivergence(BaseRelationType& base_relation);
    virtual ~CalculateVelocityDivergence() {};

protected:
    Vecd *vel_;
    Real *velocity_divergence_;
};

template <>
class CalculateVelocityDivergence<Inner<>> : public CalculateVelocityDivergence<Base, DataDelegateInner>
{
public:
    explicit CalculateVelocityDivergence(BaseInnerRelation& inner_relation);
    virtual ~CalculateVelocityDivergence() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:

};

template <>
class CalculateVelocityDivergence<Contact<>> : public CalculateVelocityDivergence<Base, DataDelegateContact>
{
public:
    explicit CalculateVelocityDivergence(BaseContactRelation& contact_relation);
    virtual ~CalculateVelocityDivergence() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:

};

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseCalculateVelocityDivergenceComplex = ComplexInteraction<CalculateVelocityDivergence<InnerInteractionType, ContactInteractionTypes...>>;

using CalculateVelocityDivergenceComplex = BaseCalculateVelocityDivergenceComplex<Inner<>, Contact<>>;
//=================================================================================================//
class CalculatePhysicalCoefficients : public LocalDynamics
{
public:
    explicit CalculatePhysicalCoefficients(SPHBody& sph_body, Real bulk_modulus, Real bulk_viscosity);
    virtual ~CalculatePhysicalCoefficients() {};

    void update(size_t index_i, Real dt = 0.0);

protected:
    Real *bulk_modulus_, *bulk_viscosity_;
    Real *Vol_, *rho_, *mass_;
    Real *volume_strain_;
    Real bulk_modulus_ref_;
};
//=================================================================================================//
class CalculatePressure : public LocalDynamics
{
public:
    explicit CalculatePressure(SPHBody& sph_body);
    virtual ~CalculatePressure() {};

    void update(size_t index_i, Real dt = 0.0);

protected:
    Real *bulk_modulus_, *bulk_viscosity_;
    Real *volume_strain_;
    Real *p_;
    Real *velocity_divergence_;

};
//=================================================================================================//
template <typename... InteractionTypes>
class CalculatePressureGradientForce;

template <class DataDelegationType>
class CalculatePressureGradientForce<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
public:
    template <class BaseRelationType>
    explicit CalculatePressureGradientForce(BaseRelationType& base_relation);
    virtual ~CalculatePressureGradientForce() {};

protected:
    Vecd *force_;
    Real *p_;
    Real *Vol_;
};

template <>
class CalculatePressureGradientForce<Inner<>> : public CalculatePressureGradientForce<Base, DataDelegateInner>
{
public:
    explicit CalculatePressureGradientForce(BaseInnerRelation& inner_relation);
    virtual ~CalculatePressureGradientForce() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:

};

template <>
class CalculatePressureGradientForce<Contact<>> : public CalculatePressureGradientForce<Base, DataDelegateContact>
{
public:
    explicit CalculatePressureGradientForce(BaseContactRelation& contact_relation);
    virtual ~CalculatePressureGradientForce() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:

};

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseCalculatePressureGradientForceComplex = ComplexInteraction<CalculatePressureGradientForce<InnerInteractionType, ContactInteractionTypes...>>;

using CalculatePressureGradientForceComplex = BaseCalculatePressureGradientForceComplex<Inner<>, Contact<>>;
//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
#endif // MPH_METHOD_H