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
 * @file 	k-epsilon_turbulent_model.hpp
 * @brief 	k-epsilon_turbulent_model.hpp
 * @author	 Xiangyu Hu
 */
#pragma once

#include "udf_k-epsilon_turbulent_model.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
namespace udf
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
kEpsilon_BaseTurbulentModel<Base, DataDelegationType>::kEpsilon_BaseTurbulentModel(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      turbu_strain_rate_(this->particles_->template registerStateVariableData<Matd>("TurbulentStrainRate")),
      viscosity_(DynamicCast<Viscosity>(this, this->particles_->getBaseMaterial())),
      mu_(viscosity_.ReferenceViscosity()),
      smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      particle_spacing_min_(base_relation.getSPHBody().getSPHAdaptation().MinimumSpacing()),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      dimension_(2) {}
//** A temporarily treatment for dimension **
//=================================================================================================//
} // namespace udf
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//