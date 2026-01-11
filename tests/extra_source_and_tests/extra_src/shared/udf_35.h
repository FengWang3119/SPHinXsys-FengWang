/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	zeroth-order_residue.h
 * @brief 	
 * @details     
 * @author Xiangyu Hu
 */

#ifndef UDF_35_H
#define UDF_35_H

#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
//=================================================================================================//
class GetLimiterOfTransportVelocityCorrection : public LocalDynamics
{
  public:
    explicit GetLimiterOfTransportVelocityCorrection(SPHBody &sph_body, Real slope = 100.0);
    virtual ~GetLimiterOfTransportVelocityCorrection(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    const Real h_ref_;
    Vecd *zero_gradient_residue_;
    Real slope_;
    Real *limiter_tvc_;
};
//=================================================================================================//
class GetPressureGradientResidue : public LocalDynamics
{
  public:
    explicit GetPressureGradientResidue(SPHBody &sph_body);
    virtual ~GetPressureGradientResidue(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *zero_gradient_residue_;
    Real *p_;
    Vecd *pressure_gradient_residue_;
};
//=================================================================================================//
template <typename... InteractionTypes>
class GetPressureGradientResidue_RKGC;

template <class DataDelegationType>
class GetPressureGradientResidue_RKGC<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit GetPressureGradientResidue_RKGC(BaseRelationType &base_relation);
    virtual ~GetPressureGradientResidue_RKGC(){};

  protected:
    Real *p_;
    Vecd *pressure_gradient_residue_RKGC_;
    Matd *B_;
    Real *Vol_;
};
//** Inner part *
template <>
class GetPressureGradientResidue_RKGC<Inner<>> : public GetPressureGradientResidue_RKGC<Base, DataDelegateInner>
{
  public:
    explicit GetPressureGradientResidue_RKGC(BaseInnerRelation &inner_relation);
    virtual ~GetPressureGradientResidue_RKGC(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
};
//** Wall part *
template <>
class GetPressureGradientResidue_RKGC<Contact<>> : public GetPressureGradientResidue_RKGC<Base, DataDelegateContact>
{
  public:
    explicit GetPressureGradientResidue_RKGC(BaseContactRelation &contact_relation);
    virtual ~GetPressureGradientResidue_RKGC(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
};

//** Interface part *
template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseGetPressureGradientResidueComplex_RKGC = ComplexInteraction<GetPressureGradientResidue_RKGC<InnerInteractionType, ContactInteractionTypes...>>;
using GetPressureGradientResidueComplex_RKGC = BaseGetPressureGradientResidueComplex_RKGC<Inner<>, Contact<>>;
//=================================================================================================//
class NonDimensionalisePressure : public LocalDynamics
{
  public:
    explicit NonDimensionalisePressure(SPHBody &sph_body);
    virtual ~NonDimensionalisePressure(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *rho_;
    Real *p_;
    Real *p_dimensionless_;
};
class DisposerInBufferDeletion : public BaseLocalDynamics<BodyPartByCell>
{
  public:
    DisposerInBufferDeletion(AlignedBoxByCell &aligned_box_part);
    virtual ~DisposerInBufferDeletion(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    std::mutex mutex_switch_to_buffer_; /**< mutex exclusion for memory conflict */
    Vecd *pos_;
    AlignedBox &aligned_box_;
};

class InitialiseColorIndicator : public LocalDynamics
{
  public:
    explicit InitialiseColorIndicator(SPHBody &sph_body);
    virtual ~InitialiseColorIndicator(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int *color_indicator_;
    Vecd *pos_;
};
class InitialiseColorIndicator2 : public LocalDynamics
{
  public:
    explicit InitialiseColorIndicator2(SPHBody &sph_body, const StdVec<Vecd> &box);
    virtual ~InitialiseColorIndicator2() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int *color_indicator_;
    Vecd *pos_;
    const StdVec<Vecd> &box_;
};

class ClearBufferParticleIndicator : public LocalDynamics
{
  public:
    explicit ClearBufferParticleIndicator(SPHBody &sph_body, int third_dimension, Real lower_bound, Real upper_bound);
    virtual ~ClearBufferParticleIndicator(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int *buffer_particle_indicator_;
    int third_dimension_;
    Real lower_bound_;
    Vecd *pos_;
    Real upper_bound_;
};

class DisposerForInitialParticleDeletion : public LocalDynamics
{
  public:
    DisposerForInitialParticleDeletion(SPHBody &sph_body);
    virtual ~DisposerForInitialParticleDeletion(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *pos_;
    std::mutex mutex_switch_to_buffer_; /**< mutex exclusion for memory conflict */
    int *buffer_particle_indicator_;
};
class DisposerForSplashParticleDeletion : public LocalDynamics
{
  public:
    DisposerForSplashParticleDeletion(SPHBody &sph_body);
    virtual ~DisposerForSplashParticleDeletion(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *pos_;
    std::mutex mutex_switch_to_buffer_; /**< mutex exclusion for memory conflict */
    Real *pos_div_;
};
//=================================================================================================//
class TagMixedParticle : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit TagMixedParticle(BaseInnerRelation &inner_relation, Real mixing_rate_interactive_radius);
    virtual ~TagMixedParticle() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    int *is_mixed_;
    int *color_indicator_;
    Real interactive_radius_;
};
//=================================================================================================//
class CalculateFluidParticleNumberInOutletChannel : public LocalDynamicsReduce<ReduceSum<int>>
{
  public:
    explicit CalculateFluidParticleNumberInOutletChannel(SPHBody &sph_body, Real radius_chamber, Real h_inlet);
    virtual ~CalculateFluidParticleNumberInOutletChannel() {};

    int reduce(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *pos_;
    Real radius_chamber_;
    Real height_inlet_channel_;
};
//=================================================================================================//
class CalculateMixedFluidParticleNumberInOutletChannel : public LocalDynamicsReduce<ReduceSum<int>>
{
  public:
    explicit CalculateMixedFluidParticleNumberInOutletChannel(SPHBody &sph_body, Real radius_chamber, Real h_inlet);
    virtual ~CalculateMixedFluidParticleNumberInOutletChannel() {};

    int reduce(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *pos_;
    int *is_mixed_;
    Real radius_chamber_;
    Real height_inlet_channel_;
};
//=================================================================================================//
class UpdateVolumeAndAddIndicatorAsEvolving : public LocalDynamics
{
public:
    explicit UpdateVolumeAndAddIndicatorAsEvolving(SPHBody& sph_body);
    virtual ~UpdateVolumeAndAddIndicatorAsEvolving() {};

    void update(size_t index_i, Real dt = 0.0);

protected:
    int *indicator_;
    Real *rho_, *mass_, *Vol_;
};
//=================================================================================================//
class ClearColorIndex : public LocalDynamics
{
public:
    explicit ClearColorIndex(SPHBody& sph_body);
    virtual ~ClearColorIndex() {};

    void update(size_t index_i, Real dt = 0.0);

protected:
    int* color_indicator_;
};
//=================================================================================================//
} // namespace SPH
#endif // K_EPSILON_TURBULENT_MODEL_H