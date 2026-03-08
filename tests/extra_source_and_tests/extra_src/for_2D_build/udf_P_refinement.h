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
    class P_refinement_GetVelocityGradientInnerOnlyP : public LocalDynamics, public DataDelegateInner
    {
    public:
        explicit P_refinement_GetVelocityGradientInnerOnlyP(BaseInnerRelation& inner_relation);
        virtual ~P_refinement_GetVelocityGradientInnerOnlyP() {};
        void interaction(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);

    protected:
        Matd* velocity_gradient_inner_only_P_;
        //
        Matd* turbu_B_;
        Real* Vol_;
        Vecd* vel_;
        int* is_near_wall_P1_;
    };
//=================================================================================================//
    class P_refinement : public LocalDynamics, public kOmega_BaseTurbuClosureCoeff
    {
    public:
        explicit P_refinement(SPHBody& sph_body);
        virtual ~P_refinement(){};

        void update(size_t index_i, Real dt = 0.0);
        double solve_1D_sublayer(double kinematic_viscosity, double u_p_outer, double k_p_outer,
            double w_p_outer, double vel_grad_p_outer, double nut_p_outer, double h_sublayer, 
            double utau_outer, double Q_target);
        std::vector<double> tdma(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d);
        void tdma5(const double a[5], const double b[5], const double c[5], const double d[5], double x[5]);

        inline Real obtainTangentialComponent(const Vecd& vec, const Vecd& normal)
        {
            Real dot_un = vec.dot(normal);                  
            Real norm_sqr = vec.dot(vec) - dot_un * dot_un;   
            return std::sqrt(std::max(0.0, norm_sqr));          
        }

    protected:
        int num_sub_node_;
        Real* friction_velocity_from_sublayer_;
        //
        int* is_near_wall_P1_;
        Real* y_p_;
        Vecd* vel_;
        Real* turbu_k_;
        Real* turbu_omega_;
        Real* rho_;
        Viscosity& viscosity_;
        Real mu_;
        Matd* velocity_gradient_inner_only_P_;
        Real* turbu_mu_;
        Real* distance_to_dummy_interface_;
        Real* wall_shear_stress_;
        Vecd* e_nearest_normal_;
        Real fluid_particle_spacing_;
    };
//=================================================================================================//
//=================================================================================================//
} // namespace udf
} // namespace fluid_dynamics
} // namespace SPH
#endif // UDF_COMMON_TURBULENCE_MODEL_H