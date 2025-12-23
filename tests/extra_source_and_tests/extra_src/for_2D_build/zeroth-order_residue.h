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

#ifndef ZEROTH_ORDER_RESIDUE_H
#define ZEROTH_ORDER_RESIDUE_H

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
} // namespace SPH
#endif // K_EPSILON_TURBULENT_MODEL_H