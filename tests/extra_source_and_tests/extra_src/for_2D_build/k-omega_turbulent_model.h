#ifndef K_OMEGA_TURBULENT_MODEL_H
#define K_OMEGA_TURBULENT_MODEL_H

#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
class BaseTurbuClosureCoeff
{
  public:
    explicit BaseTurbuClosureCoeff();
    virtual ~BaseTurbuClosureCoeff(){};

  protected:
    Real Karman_;
    Real turbu_const_E_;
    Real C_mu_, C_mu_25_, C_mu_75_;
    Real turbulent_intensity_;

    //** Closure coefficients for K *
    Real sigma_k_;

    //** Closure coefficients for Epsilon *
    Real C_l_, C_2_;
    Real sigma_E_;
    Real turbulent_length_ratio_for_epsilon_inlet_;

    //** Start time for laminar law *
    Real start_time_laminar_;
    Real y_star_threshold_laminar_;

    //** Closure coefficients for Omega *
    Real std_kw_sigma_k_;
    Real std_kw_sigma_omega_;
    Real std_kw_beta_star_;
    Real std_kw_sigma_star_;
    Real std_kw_alpha_;
    Real std_kw_sigma_;
    Real std_kw_beta_;
    Real std_kw_f_beta_; //** Temporarily treat for 2d */
    Real std_kw_beta_0_;
    Real std_kw_sigma_do_;
    Real std_kw_sigma_d_;
    Real std_kw_C_lim_;
    Real std_kw_beta_i_;
};
//=================================================================================================//
class WallFunction : public BaseTurbuClosureCoeff
{
  public:
    explicit WallFunction(){};
    virtual ~WallFunction(){};

    Real get_dimensionless_velocity(Real y_star);
    Real get_near_wall_velocity_gradient_magnitude(Real y_star, Real vel_fric_mag, Real denominator_log_law, Real dynamic_viscosity);
    Real get_distance_from_P_to_wall(Real y_p_constant);

    Real log_law_wall_function(Real y_star);
    Real laminar_law_wall_function(Real y_star);
    Real log_law_velocity_gradient(Real vel_fric_mag, Real denominator_log_law);
    Real laminar_law_velocity_gradient(Real vel_fric_mag, Real dynamic_viscosity);
};
//=================================================================================================//
template <typename... InteractionTypes>
class GetVelocityGradient;

template <class DataDelegationType>
class GetVelocityGradient<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit GetVelocityGradient(BaseRelationType &base_relation);
    virtual ~GetVelocityGradient(){};

  protected:
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Vecd> &vel_, &pos_;
    StdLargeVec<int> &is_near_wall_P1_; //** This is used to specially treat near wall region  *
    StdLargeVec<int> &is_near_wall_P2_;

    StdLargeVec<Matd> &velocity_gradient_;
    //**For test*
    StdLargeVec<Matd> &velocity_gradient_wall;
};
//** Inner part *
template <>
class GetVelocityGradient<Inner<>> : public GetVelocityGradient<DataDelegateInner>
{
  public:
    explicit GetVelocityGradient(BaseInnerRelation &inner_relation, Real weight_sub);
    virtual ~GetVelocityGradient(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &velocity_gradient_;
    StdLargeVec<Matd> &B_;
    StdLargeVec<Matd> &turbu_B_;
    Real weight_sub_nearwall_;
};
using GetVelocityGradientInner = GetVelocityGradient<Inner<>>;

//** Updated Wall part *
template <>
class GetVelocityGradient<Contact<Wall>> : public InteractionWithWall<GetVelocityGradient>
{
  public:
    explicit GetVelocityGradient(BaseContactRelation &contact_relation);
    virtual ~GetVelocityGradient(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &velocity_gradient_;
};

//** Interface part *
using GetVelocityGradientComplex = ComplexInteraction<GetVelocityGradient<Inner<>, Contact<Wall>>>;

//using GetVelocityGradientComplex = BaseGetVelocityGradientComplex<Inner<>, Contact<>>;
//=================================================================================================//
template <typename... T>
class BaseTurbulentModel;

template <class DataDelegationType>
class BaseTurbulentModel<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType, public BaseTurbuClosureCoeff
{
  public:
    template <class BaseRelationType>
    explicit BaseTurbulentModel(BaseRelationType &base_relation);
    virtual ~BaseTurbulentModel(){};

  protected:
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_epsilon_;
    StdLargeVec<Real> &turbu_mu_;
    StdLargeVec<Matd> &turbu_strain_rate_; //** temporary naming to distinguish the regular strain rate *
    StdLargeVec<Real> &turbu_strain_rate_magnitude_;
    StdLargeVec<Real> &turbu_omega_;

    Real mu_, smoothing_length_, particle_spacing_min_;
    StdLargeVec<Real> &rho_, &Vol_;
    StdLargeVec<Vecd> &vel_;
    int dimension_;
};
//=================================================================================================//
/**
	 * @class kOmegaSST_TurbulentModelInner
	 * @brief  kOmegaSST_TurbulentModelInner
	 */
class kOmega_kTransportEquationInner : public BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kOmega_kTransportEquationInner(BaseInnerRelation &inner_relation, const StdVec<Real> &initial_values, int is_extr_visc_dissipa = 0);
    virtual ~kOmega_kTransportEquationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &dk_dt_;
    StdLargeVec<Real> &dk_dt_without_dissipation_;
    StdLargeVec<Real> &k_production_;

    StdLargeVec<int> &is_near_wall_P1_; //** This is used to specially treat near wall region  *
    StdLargeVec<Matd> &velocity_gradient_;
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_omega_;
    StdLargeVec<Real> &turbu_mu_;
    StdLargeVec<Matd> &turbu_strain_rate_;
    StdLargeVec<Real> &turbu_strain_rate_magnitude_;
    StdLargeVec<int> &is_extra_viscous_dissipation_;

    //** for test */
    StdLargeVec<int> &turbu_indicator_;
    StdLargeVec<Real> &k_diffusion_, &vel_x_;
};
//=================================================================================================//
/**
	 * @class kOmegaSST_TurbulentModelInner
	 * @brief  kOmegaSST_TurbulentModelInner
	 */
class kOmega_omegaTransportEquationInner : public BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kOmega_omegaTransportEquationInner(BaseInnerRelation &inner_relation);
    virtual ~kOmega_omegaTransportEquationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &domega_dt_, &domega_dt_without_dissipation_;
    StdLargeVec<Real> &omega_production_, &omega_dissipation_, &omega_diffusion_;
    StdLargeVec<Real> &omega_cross_diffusion_;

    StdLargeVec<Real> &turbu_mu_;
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_omega_;
    StdLargeVec<Real> &k_production_;
    StdLargeVec<int> &is_near_wall_P1_;
};
//=================================================================================================//

template <typename... InteractionTypes>
class TKEnergyForce;

template <class DataDelegationType>
class TKEnergyForce<Base, DataDelegationType>
    : public BaseTurbulentModel<Base, DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit TKEnergyForce(BaseRelationType &base_relation);
    virtual ~TKEnergyForce(){};

  protected:
    StdLargeVec<Vecd> &force_;
    StdLargeVec<Real> &mass_;
    StdLargeVec<int> &indicator_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Vecd> &test_k_grad_rslt_;

    //StdLargeVec<Vecd> &tke_acc_inner_, &tke_acc_wall_;
};
//** Inner part *
template <>
class TKEnergyForce<Inner<>> : public TKEnergyForce<Base, DataDelegateInner>
{
  public:
    explicit TKEnergyForce(BaseInnerRelation &inner_relation);
    virtual ~TKEnergyForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &test_k_grad_rslt_;
    StdLargeVec<Matd> &B_;
};
//** Wall part *
template <>
class TKEnergyForce<Contact<>> : public TKEnergyForce<Base, DataDelegateContact>
{
  public:
    explicit TKEnergyForce(BaseContactRelation &contact_relation);
    //: TKEnergyForce<Base, DataDelegateContact>(contact_relation) {};
    virtual ~TKEnergyForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &test_k_grad_rslt_;
    StdLargeVec<Matd> &B_;
};

//** Interface part *
template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseTKEnergyForceComplex = ComplexInteraction<TKEnergyForce<InnerInteractionType, ContactInteractionTypes...>>;

using TKEnergyForceComplex = BaseTKEnergyForceComplex<Inner<>, Contact<>>;
//=================================================================================================//
template <typename... InteractionTypes>
class TurbuViscousForce;

template <class DataDelegationType>
class TurbuViscousForce<DataDelegationType> : public ViscousForce<DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit TurbuViscousForce(BaseRelationType &base_relation);
    virtual ~TurbuViscousForce(){};

  protected:
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_mu_;
    StdLargeVec<Real> &wall_Y_plus_;
    StdLargeVec<Real> &wall_Y_star_;
    StdLargeVec<Vecd> &velo_friction_;
    StdLargeVec<Real> &y_p_;
    StdLargeVec<int> &is_near_wall_P2_;
    Real molecular_viscosity_;
    Real c0_;

    //** For test *
    //StdLargeVec<Matd> visc_direction_matrix_;
    //StdLargeVec<Vecd> &visc_acc_inner_, &visc_acc_wall_;
};

//** Inner part *
template <>
class TurbuViscousForce<Inner<>> : public TurbuViscousForce<DataDelegateInner>, public ForcePrior
{
  public:
    explicit TurbuViscousForce(BaseInnerRelation &inner_relation);
    virtual ~TurbuViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<int> &turbu_indicator_;
    StdLargeVec<int> &is_extra_viscous_dissipation_;
    StdLargeVec<Matd> &B_;
};

//** Wall part *
using BaseTurbuViscousForceWithWall = InteractionWithWall<TurbuViscousForce>;
template <>
class TurbuViscousForce<Contact<Wall>> : public BaseTurbuViscousForceWithWall, public WallFunction
{
  public:
    explicit TurbuViscousForce(BaseContactRelation &wall_contact_relation);
    virtual ~TurbuViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real wall_particle_spacing_;
    StdLargeVec<Matd> &B_;
};

//** Interface part *
using TurbulentViscousForceWithWall = ComplexInteraction<TurbuViscousForce<Inner<>, Contact<Wall>>>;
//=================================================================================================//
class kOmegaTurbulentEddyViscosity : public LocalDynamics,
                                     public DataDelegateSimple,
                                     public BaseTurbuClosureCoeff
{
  public:
    explicit kOmegaTurbulentEddyViscosity(SPHBody &sph_body);
    virtual ~kOmegaTurbulentEddyViscosity(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &rho_;
    StdLargeVec<Real> &turbu_mu_;
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_omega_;
    StdLargeVec<Real> &wall_Y_plus_, &wall_Y_star_;
    StdLargeVec<Real> &turbu_strain_rate_magnitude_;
    Real mu_;
};
//=================================================================================================//
/**
	 * @class TurbulentAdvectionTimeStepSize
	 * @brief Computing the turbulent advection time step size
	 */
class TurbulentAdvectionTimeStepSize : public LocalDynamicsReduce<ReduceMax>,
                                       public DataDelegateSimple
{
  public:
    explicit TurbulentAdvectionTimeStepSize(SPHBody &sph_body, Real U_max, Real advectionCFL = 0.25);
    virtual ~TurbulentAdvectionTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    StdLargeVec<Vecd> &vel_;
    Real smoothing_length_min_;
    Real speed_ref_turbu_, advectionCFL_;
    StdLargeVec<Real> &turbu_mu_;
    Fluid &fluid_;
};
//=================================================================================================//
/**
	* @class   InflowTurbulentCondition
	* @brief   Inflow boundary condition which imposes directly to a given velocity profile.
	*          TargetVelocity gives the velocity profile along the inflow direction,
	*          i.e. x direction in local frame.
	*/
class InflowTurbulentCondition : public BaseFlowBoundaryCondition, public BaseTurbuClosureCoeff
{
  public:
    explicit InflowTurbulentCondition(BodyPartByCell &body_part,
                                      Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet);
    virtual ~InflowTurbulentCondition(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int type_turbu_inlet_;
    Real relaxation_rate_;
    Real CharacteristicLength_;
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_epsilon_;
    Real TurbulentLength_;

    virtual Real getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k);
    virtual Real getTurbulentInflowE(Vecd &position, Real &turbu_k, Real &turbu_E);
};
//=================================================================================================//
class JudgeIsNearWall : public LocalDynamics, public DataDelegateContact, public BaseTurbuClosureCoeff
{
  public:
    JudgeIsNearWall(BaseInnerRelation &inner_relation,
                    BaseContactRelation &contact_relation);
    virtual ~JudgeIsNearWall(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &distance_to_dummy_interface_;
    StdLargeVec<Real> &distance_to_dummy_interface_up_average_;
    StdLargeVec<int> &is_near_wall_P1_;
    StdLargeVec<int> &is_near_wall_P2_;
    StdLargeVec<int> &index_nearest_;
    StdLargeVec<Vecd> &e_nearest_tau_, &e_nearest_normal_;

    StdLargeVec<Vecd> &pos_;
    int dimension_;
    Real fluid_particle_spacing_, wall_particle_spacing_;
    StdLargeVec<Vecd> &distance_from_wall_;
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Vecd> *> contact_n_;
};
//=================================================================================================//
class kOmegaStdWallFuncCorrection : public LocalDynamics, public DataDelegateContact, public WallFunction
{
  public:
    kOmegaStdWallFuncCorrection(BaseInnerRelation &inner_relation,
                                BaseContactRelation &contact_relation, Real y_p_constant);
    virtual ~kOmegaStdWallFuncCorrection(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    //Real offset_dist_;
    StdLargeVec<Real> &y_p_;
    StdLargeVec<Real> &wall_Y_plus_, &wall_Y_star_;
    StdLargeVec<Real> &velo_tan_;
    StdLargeVec<Vecd> &velo_friction_;

    StdLargeVec<Vecd> &vel_;
    StdLargeVec<Real> &rho_;
    Real molecular_viscosity_;
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_epsilon_;
    StdLargeVec<Real> &turbu_omega_;
    StdLargeVec<Real> &turbu_mu_;
    StdLargeVec<int> &is_near_wall_P1_;
    StdLargeVec<int> &is_near_wall_P2_;
    StdLargeVec<Matd> &velocity_gradient_;
    StdLargeVec<Real> &k_production_;
    StdLargeVec<Real> &distance_to_dummy_interface_;
    StdLargeVec<Real> &distance_to_dummy_interface_up_average_;
    StdLargeVec<int> &index_nearest;
    StdLargeVec<Vecd> &e_nearest_tau_;
    StdLargeVec<Vecd> &e_nearest_normal_;
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Vecd> *> contact_n_;
};
//=================================================================================================//
class ConstrainNormalVelocityInRegionP : public LocalDynamics,
                                         public DataDelegateSimple
{
  public:
    explicit ConstrainNormalVelocityInRegionP(SPHBody &sph_body);
    virtual ~ConstrainNormalVelocityInRegionP(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<int> &is_near_wall_P1_;
    StdLargeVec<Vecd> &e_nearest_normal_;
};
//=================================================================================================//
template <int INDICATOR>
class TurbulentIndicatedParticles : public WithinScope
{
    StdLargeVec<int> &indicator_;
    StdLargeVec<int> &turbu_plug_flow_indicator_;

  public:
    explicit TurbulentIndicatedParticles(BaseParticles *base_particles)
        : WithinScope(),
          indicator_(*base_particles->getVariableByName<int>("Indicator")),
          turbu_plug_flow_indicator_(*base_particles->getVariableByName<int>("TurbulentPlugFlowIndicator")){};
    bool operator()(size_t index_i)
    {
        return indicator_[index_i] == INDICATOR && turbu_plug_flow_indicator_[index_i] == INDICATOR;
    };
};

using TurbulentPlugFlowParticles = TurbulentIndicatedParticles<0>;
//=================================================================================================//
template <typename... InteractionTypes>
class TurbulentLinearGradientCorrectionMatrix;

template <class DataDelegationType>
class TurbulentLinearGradientCorrectionMatrix<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit TurbulentLinearGradientCorrectionMatrix(BaseRelationType &base_relation);
    virtual ~TurbulentLinearGradientCorrectionMatrix(){};

  protected:
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Matd> &turbu_B_;
    StdLargeVec<Matd> &B_;
};

template <>
class TurbulentLinearGradientCorrectionMatrix<Inner<>>
    : public TurbulentLinearGradientCorrectionMatrix<DataDelegateInner>
{
    Real turbu_alpha_;

  public:
    explicit TurbulentLinearGradientCorrectionMatrix(BaseInnerRelation &inner_relation, Real alpha = Real(0))
        : TurbulentLinearGradientCorrectionMatrix<DataDelegateInner>(inner_relation), turbu_alpha_(alpha){};
    template <typename BodyRelationType, typename FirstArg>
    explicit TurbulentLinearGradientCorrectionMatrix(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : TurbulentLinearGradientCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~TurbulentLinearGradientCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using TurbulentLinearGradientCorrectionMatrixInner = TurbulentLinearGradientCorrectionMatrix<Inner<>>;

//=================================================================================================//
class GetLimiterOfTransportVelocityCorrection : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit GetLimiterOfTransportVelocityCorrection(SPHBody &sph_body, Real slope);
    virtual ~GetLimiterOfTransportVelocityCorrection(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    const Real h_ref_;
    StdLargeVec<Vecd> &zero_gradient_residue_;
    Real slope_;
    StdLargeVec<Real> &limiter_tvc_;
};
//=================================================================================================//
template <class ParticleScope>
using TVC_Limited_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, TruncatedLinear, LinearGradientCorrection, ParticleScope>;
template <class ParticleScope>
using TVC_NoLimiter_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, NoLimiter, LinearGradientCorrection, ParticleScope>;
//=================================================================================================//
class ModifiedTruncatedLinear : public Limiter
{
    Real ref_, slope_;

  public:
    ModifiedTruncatedLinear(Real ref, Real slope = 1000.0)
        : Limiter(), ref_(ref), slope_(slope){};
    Real operator()(Real measure)
    {
        Real measure_scale = measure * ref_;
        return SMIN(slope_ * measure_scale, Real(1));
    };
};
template <class ParticleScope>
using TVC_ModifiedLimited_NoRKGC =
    BaseTransportVelocityCorrectionComplex<SingleResolution, ModifiedTruncatedLinear, NoKernelCorrection, ParticleScope>;

template <class ParticleScope>
using TVC_ModifiedLimited_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, ModifiedTruncatedLinear, LinearGradientCorrection, ParticleScope>;

template <class ParticleScope>
using TVC_ModifiedLimited_RKGC_OBFCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, ModifiedTruncatedLinear, LinearGradientCorrectionWithBulkScope, ParticleScope>;

template <class ParticleScope>
using TVC_NotLimited_RKGC_OBFCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, NoLimiter, LinearGradientCorrectionWithBulkScope, ParticleScope>;

template <class ParticleScope>
using TVC_ModifiedLimited_withoutLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, ModifiedTruncatedLinear, NoKernelCorrection, ParticleScope>;
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // K_EPSILON_TURBULENT_MODEL_H