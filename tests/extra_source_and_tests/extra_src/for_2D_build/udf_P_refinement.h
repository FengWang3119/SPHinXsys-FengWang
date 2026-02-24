#ifndef UDF_P_REFINEMENT_H
#define UDF_P_REFINEMENT_H

#include "udf_k-omega_turbulent_model.h"
#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
namespace udf
{
//=================================================================================================//
    class P_refinement : public LocalDynamics, public kOmega_BaseTurbuClosureCoeff
    {
    public:
        explicit P_refinement(SPHBody& sph_body);
        virtual ~P_refinement(){};

        void update(size_t index_i, Real dt = 0.0);

    protected:
        int num_sub_node_;
        //
        int* is_near_wall_P1_;
        Real* y_p_;
        Vecd* vel_;
        Real* turbu_k_;
        Real* turbu_omega_;
        Real* rho_;
        Viscosity& viscosity_;
        Real mu_;
    };
//=================================================================================================//
} // namespace udf
} // namespace fluid_dynamics
} // namespace SPH
#endif // UDF_COMMON_TURBULENCE_MODEL_H