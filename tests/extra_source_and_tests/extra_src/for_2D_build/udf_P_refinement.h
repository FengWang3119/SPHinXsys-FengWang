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
    template <typename... InteractionTypes>
    class P_refinement_GetVelocityGradient;

    template <class DataDelegationType>
    class P_refinement_GetVelocityGradient<DataDelegationType>
        : public LocalDynamics, public DataDelegationType
    {
    public:
        template <class BaseRelationType>
        explicit P_refinement_GetVelocityGradient(BaseRelationType& base_relation);
        virtual ~P_refinement_GetVelocityGradient() {};

    protected:
        Matd* velocity_gradient_only_P_;
        //
        Real* Vol_;
        Vecd* vel_;
        int* is_near_wall_P1_;
        int* is_near_wall_P2_;
    };

    //** Inner part *
    template <>
    class P_refinement_GetVelocityGradient<Inner<>> : public P_refinement_GetVelocityGradient<DataDelegateInner>
    {
    public:
        explicit P_refinement_GetVelocityGradient(BaseInnerRelation& inner_relation);
        virtual ~P_refinement_GetVelocityGradient() {};
        void interaction(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);

    protected:
        Matd* turbu_B_;
        Matd* B_;
    };
    using P_refinement_GetVelocityGradientInner = P_refinement_GetVelocityGradient<Inner<>>;

    //** Wall part *
    template <>
    class P_refinement_GetVelocityGradient<Contact<Wall>> : public InteractionWithWall<P_refinement_GetVelocityGradient>
    {
    public:
        explicit P_refinement_GetVelocityGradient(BaseContactRelation& contact_relation);
        virtual ~P_refinement_GetVelocityGradient() {};
        void interaction(size_t index_i, Real dt = 0.0);

    protected:

    };

    //** Interface part *
    using P_refinement_GetVelocityGradientComplex = ComplexInteraction<P_refinement_GetVelocityGradient<Inner<>, Contact<Wall>>>;
//=================================================================================================//
    class P_refinement : public LocalDynamics, public kOmega_BaseTurbuClosureCoeff, public WallFunctionCoefficient
    {
    public:
        explicit P_refinement(SPHBody& sph_body);
        virtual ~P_refinement(){};

        void update(size_t index_i, Real dt = 0.0);
        Vec6d solve_1D_sublayer(double kinematic_viscosity, double u_p_outer, double k_p_outer,
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
        Real* target_flow_rate_in_sublayer_;
        Real* vel_ps_magnitude_;
        Real* dudn_;
        Real* utau_node_;
        Vec6d* node_value_; // ** Temporary treatment only valid for 5-node configuration, first is utau, then velocity *
        Vecd* node_vel_first_second_; // ** Temporary treatment
        Vecd* node_vel_third_fourth_; // ** Temporary treatment
        Real* node_vel_fifth_; // ** Temporary treatment
        Real* dUdn_P_sublayer_magnitude_;
        Matd* dUdn_P_sublayer_;
        //
        int* is_near_wall_P1_;
        Real* y_p_;
        Vecd* vel_;
        Real* turbu_k_;
        Real* turbu_omega_;
        Real* rho_;
        Viscosity& viscosity_;
        Real mu_;
        Matd* velocity_gradient_only_P_;
        Real* turbu_mu_;
        Real* distance_to_dummy_interface_;
        Real* wall_shear_stress_;
        Vecd* e_nearest_normal_;
        Real fluid_particle_spacing_;
        Real* physical_time_;
        Matd* velocity_gradient_;
    };
//=================================================================================================//
//=================================================================================================//
} // namespace udf
} // namespace fluid_dynamics
} // namespace SPH
#endif // UDF_COMMON_TURBULENCE_MODEL_H