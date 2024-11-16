#ifndef K_OMEGA_TURBULENT_MODEL_H
#define K_OMEGA_TURBULENT_MODEL_H

#include "k-epsilon_turbulent_model.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
class kOmega_BaseTurbuClosureCoeff
{
  public:
    explicit kOmega_BaseTurbuClosureCoeff();
    virtual ~kOmega_BaseTurbuClosureCoeff(){};

  protected:
    //** Closure coefficients for Omega *
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
    Real std_kw_beta_star_25_;
};
//=================================================================================================//
template <typename... T>
class kOmega_BaseTurbulentModel;

template <class DataDelegationType>
class kOmega_BaseTurbulentModel<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType, public kOmega_BaseTurbuClosureCoeff
{
  public:
    template <class BaseRelationType>
    explicit kOmega_BaseTurbulentModel(BaseRelationType &base_relation);
    virtual ~kOmega_BaseTurbulentModel(){};

  protected:
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_mu_;
    StdLargeVec<Matd> &turbu_strain_rate_; //** temporary naming to distinguish the regular strain rate *
    StdLargeVec<Real> &turbu_strain_rate_magnitude_;
    StdLargeVec<Real> &turbu_strain_rate_traceless_magnitude_;
    StdLargeVec<Real> &turbu_omega_;

    Real mu_, smoothing_length_, particle_spacing_min_;
    StdLargeVec<Real> &rho_;
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Vecd> &vel_;
    int dimension_;
};
//=================================================================================================//
/**
	 * @class kOmegaSST_TurbulentModelInner
	 * @brief  kOmegaSST_TurbulentModelInner
	 */
class kOmega_kTransportEquationInner : public kOmega_BaseTurbulentModel<Base, DataDelegateInner>
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
    StdLargeVec<Real> &k_diffusion_;
    StdLargeVec<Real> &vel_x_;
};
//=================================================================================================//
/**
	 * @class kOmegaSST_TurbulentModelInner
	 * @brief  kOmegaSST_TurbulentModelInner
	 */
class kOmega_omegaTransportEquationInner : public kOmega_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kOmega_omegaTransportEquationInner(BaseInnerRelation &inner_relation);
    virtual ~kOmega_omegaTransportEquationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &domega_dt_;
    StdLargeVec<Real> &domega_dt_without_dissipation_;
    StdLargeVec<Real> &omega_production_;
    StdLargeVec<Real> &omega_dissipation_;
    StdLargeVec<Real> &omega_diffusion_;
    StdLargeVec<Real> &omega_cross_diffusion_;

    StdLargeVec<Real> &turbu_mu_;
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_omega_;
    StdLargeVec<Real> &k_production_;
    StdLargeVec<int> &is_near_wall_P1_;
};
//=================================================================================================//
class kOmegaTurbulentEddyViscosity : public LocalDynamics,
                                     public DataDelegateSimple,
                                     public kOmega_BaseTurbuClosureCoeff
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
    StdLargeVec<Real> &wall_Y_plus_;
    StdLargeVec<Real> &wall_Y_star_;
    StdLargeVec<Real> &turbu_strain_rate_traceless_magnitude_;
    Real mu_;
};
//=================================================================================================//
class kOmegaStdWallFuncCorrection : public LocalDynamics,
                                    public DataDelegateContact,
                                    public WallFunction,
                                    public kOmega_BaseTurbuClosureCoeff
{
  public:
    kOmegaStdWallFuncCorrection(BaseInnerRelation &inner_relation,
                                BaseContactRelation &contact_relation, Real y_p_constant);
    virtual ~kOmegaStdWallFuncCorrection(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    //Real offset_dist_;
    StdLargeVec<Real> &y_p_;
    StdLargeVec<Real> &wall_Y_plus_;
    StdLargeVec<Real> &wall_Y_star_;
    StdLargeVec<Real> &velo_tan_;
    StdLargeVec<Vecd> &velo_friction_;

    StdLargeVec<Vecd> &vel_;
    StdLargeVec<Real> &rho_;
    Real molecular_viscosity_;
    StdLargeVec<Real> &turbu_k_;
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
/**
	* @class   InflowTurbulentCondition
	* @brief   Inflow boundary condition which imposes directly to a given velocity profile.
	*          TargetVelocity gives the velocity profile along the inflow direction,
	*          i.e. x direction in local frame.
	*/
class kOmegaInflowTurbulentCondition : public BaseFlowBoundaryCondition,
                                       public BaseTurbuClosureCoeff,
                                       public kOmega_BaseTurbuClosureCoeff
{
  public:
    explicit kOmegaInflowTurbulentCondition(BodyPartByCell &body_part,
                                            Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet);
    virtual ~kOmegaInflowTurbulentCondition(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int type_turbu_inlet_;
    Real relaxation_rate_;
    Real CharacteristicLength_;
    StdLargeVec<Real> &turbu_k_;
    StdLargeVec<Real> &turbu_omega_;
    Real TurbulentLength_;

    virtual Real getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k);
    virtual Real getTurbulentInflowTemporaryEpsilon(Vecd &position, Real &turbu_k, Real turbu_E);
};
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // K_EPSILON_TURBULENT_MODEL_H