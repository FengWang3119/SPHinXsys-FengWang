
#pragma once

#include "zeroth-order_residue.h"

namespace SPH
{
template <class DataDelegationType>
template <class BaseRelationType>
GetPressureGradientResidue_RKGC<Base, DataDelegationType>::
    GetPressureGradientResidue_RKGC(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      p_(this->particles_->template getVariableDataByName<Real>("Pressure")),
      pressure_gradient_residue_RKGC_(this->particles_->template registerStateVariable<Vecd>("PressureGradientResidueRKGC")),
      B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure"))
{
    this->particles_->template addVariableToWrite<Vecd>("PressureGradientResidueRKGC");
}

//=================================================================================================//
} // namespace SPH
  //=================================================================================================//