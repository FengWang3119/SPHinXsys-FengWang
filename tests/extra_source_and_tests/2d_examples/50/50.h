#include "bidirectional_buffer.h"
#include "udf_common_turbulence_model.cpp"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "udf_k-omega_turbulent_model.cpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "udf_P_refinement.h"
#include "sphinxsys.h"
#include "py2_20_profile.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 2.0;  /**< Channel height. */
Real DL = 30.0; /**< Channel length. */
Real num_fluid_cross_section = 160.0;

Real time_gradually_increase_vel = 2.0;
//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------
Real characteristic_length = DH; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/

//** For K and Epsilon/Omega, type of the turbulent inlet, 2 is by table *
int type_turbulent_inlet_omega = 2;
std::string turbulent_inlet_omega_profile_source = "PY2-20";
//** For K and Epsilon/Omega, type of the turbulent inlet, 2 is by table *
int type_turbulent_inlet_k = 2;
std::string turbulent_inlet_k_profile_source = "PY2-20";
// ** 3 is using the table method *
int type_velocity_inlet = 3;
std::string inlet_vel_profile_source = "PY2-20";

Real relaxation_rate_turbulent_inlet = 0.8;
//** Tag for wall treatment *
int is_blended = 1;
//** Tag for AMRD *
int is_AMRD = 0;
bool is_constrain_normal_velocity_in_P_region = false;
//** Weight for correcting the velocity  gradient in the sub near wall region  *
//Real weight_vel_grad_sub_nearwall = 0.1;
//** Tag for Source Term Linearisation *
bool is_source_term_linearisation = false;

//** Tag for Sublayer Model *
static constexpr int num_node_sublayer_model = 5;
static constexpr int type_tdma_sublayer_model = 5;

//** Empirical parameter for initial stability*
Real turbulent_module_activate_time = 2.0;
//** Initial values for K, Omega and Mu_t *
StdVec<Real> initial_turbu_values = {0.01, 2.056, 0.02};

//Real y_p_constant = 0.05;
//Real y_p_constant = DH / 2.0 / num_fluid_cross_section; //** For the first try or Not use BOT *
//Real resolution_ref = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
//Real offset_distance = y_p_constant - resolution_ref / 2.0;                        //** Basically offset distance is large than or equal to 0 *

// ** If not use BOT *
Real y_p_constant = DH / 2.0 / num_fluid_cross_section;            
Real resolution_ref = DH / num_fluid_cross_section;
Real offset_distance = 0.0;                        

Real BW = resolution_ref * 4; /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;
Real half_channel_height = DH / 2.0;
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + 2.0 * BW));

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real U_inlet = 1.0;
Real U_f = U_inlet;         //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */

Real Re = 40000.0;

Real Outlet_pressure = 0.0;

Real mu_f = rho0_f * U_f * DH / Re;

Real Re_calculated = U_f * DH * rho0_f / mu_f;

Real DH_C = DH - 2.0 * offset_distance;
//----------------------------------------------------------------------
//	The emitter block with offset model.
//----------------------------------------------------------------------
Vec2d left_buffer_halfsize = Vec2d(0.5 * BW, 0.5 * DH_C + BW);
Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize + Vecd(0.0, offset_distance - BW);

Vec2d right_buffer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d right_buffer_translation = Vec2d(DL, DH + 0.25 * DH) - right_buffer_halfsize;
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
Real x_observe_start = 0.99 * DL;
int num_observer_points = std::round(DH_C / resolution_ref); //**Every particle is regarded as a cell monitor*
Real observe_spacing = DH_C / num_observer_points;

// By kernel weight.
StdVec<Vecd> observation_location;
StdVec<Vecd> observation_theoretical_locations;
Vecd pos_observe_start = Vecd(x_observe_start, resolution_ref / 2.0 + offset_distance);
Vecd unit_direction_observe = Vecd(0.0, 1.0);
Real observer_offset_distance = 2.0 * resolution_ref;

void get_observation_locations()
{
    for (int i = 0; i < num_observer_points; ++i)
    {
        Vecd pos_observer_i = pos_observe_start + i * observe_spacing * unit_direction_observe;
        Vecd pos_observer_i_no_offset = pos_observe_start + i * observe_spacing * unit_direction_observe;
        if (i == 0)
        {
            pos_observer_i -= observer_offset_distance * unit_direction_observe;
        }
        if (i == num_observer_points - 1)
        {
            pos_observer_i += observer_offset_distance * unit_direction_observe;
        }
        observation_location.push_back(pos_observer_i);
        observation_theoretical_locations.push_back(pos_observer_i_no_offset);
    }
}
void output_observer_theoretical_y()
{
    std::cout << "offset_distance=" << offset_distance << std::endl;
    std::string filename = "../bin/output/observer_theoretical_y.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (const Vecd &position : observation_theoretical_locations)
    {
        outfile << position[1] << "\n";
    }
    outfile.close();
}
//** For regression test *
StdVec<Vecd> observer_location_center_point = {Vecd(0.5 * DL, 0.5 * DH)};

//** For getting near wall friction velocity or wall shear stress *
namespace observe_nearwall
{
constexpr const char *namespace_prefix = "nearwall";
const int number_observe_line = 1;
Real sparse_ratio = 2.0;
Real observer_offset_distance = 0.0 * resolution_ref; //** Offset the first and last observing point *
Real observer_offset_distance_whole_line = 2.0 * resolution_ref;
Vec2d unit_direction_observe(1.0, 0.0);
// ** Determine the observing start point of the each line. *
Real observe_start_x[number_observe_line] = {0.0};
Real observe_start_y[number_observe_line] = {0.5 * resolution_ref - observer_offset_distance_whole_line};
// ** Determine the length of the observing line and other information. *
Real observe_line_length[number_observe_line] = {0.0};
int num_observer_points[number_observe_line] = {0};
void getObservingLineLengthAndEndPoints()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        observe_line_length[i] = DL;
        num_observer_points[i] = std::round(observe_line_length[i] / resolution_ref / sparse_ratio);
    }
}

StdVec<Vecd> observation_locations;
StdVec<Vecd> observation_theoretical_locations;
void getPositionsOfMultipleObserveLines()
{
    getObservingLineLengthAndEndPoints();
    for (int k = 0; k < number_observe_line; ++k)
    {
        Vecd pos_observe_start(observe_start_x[k], observe_start_y[k]);
        int num_observer_point = num_observer_points[k];
        Real observe_spacing = observe_line_length[k] / num_observer_point;
        for (int i = 0; i < num_observer_point; ++i)
        {
            Real offset = 0.0;
            offset = (i == 0 ? -observer_offset_distance : (i == num_observer_point - 1 ? observer_offset_distance : 0.0));
            Vecd pos_observer_i = pos_observe_start + (i * observe_spacing + offset) * unit_direction_observe;
            Vecd pos_observer_i_no_offset = pos_observe_start + i * observe_spacing * unit_direction_observe;
            observation_locations.push_back(pos_observer_i);
            observation_theoretical_locations.push_back(pos_observer_i_no_offset);
        }
    }
}
void output_observe_positions()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_positions.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (const Vecd &position : observation_locations)
    {
        outfile << position[0] << " " << position[1] << "\n";
    }
    outfile.close();
}
void output_observe_theoretical_x()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_theoretical_x.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (const Vecd &position : observation_theoretical_locations)
    {
        outfile << position[0] << "\n";
    }
    outfile.close();
}
void output_number_observe_points_on_lines()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_num_points_on_lines.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (const int &number : num_observer_points)
    {
        outfile << number << "\n";
    }
    outfile.close();
}
} // namespace observe_nearwall

namespace observe_node_cross_section
{
    Real x_observe_node = x_observe_start;
    // By kernel weight.
    StdVec<Vecd> observation_location;
    StdVec<Vecd> observation_theoretical_locations;
    void get_observation_locations()
    {
        observation_location.push_back(Vec2d(x_observe_node, resolution_ref / 2.0 - observer_offset_distance));
        observation_location.push_back(Vec2d(x_observe_node, (DH - resolution_ref / 2.0) + observer_offset_distance));
        observation_theoretical_locations.push_back(Vec2d(x_observe_node, resolution_ref / 2.0));
        observation_theoretical_locations.push_back(Vec2d(x_observe_node, DH - resolution_ref / 2.0));
    }
    constexpr const char* namespace_prefix = "Node_";
    void output_observe_positions()
    {
        std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_positions.dat";
        std::ofstream outfile(filename);
        if (!outfile.is_open())
        {
            std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
            return;
        }
        for (const Vecd& position : observation_location)
        {
            outfile << position[0] << " " << position[1] << "\n";
        }
        outfile.close();
    }
    void output_observer_theoretical_y()
    {
        std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_theoretical_y.dat";
        std::ofstream outfile(filename);
        if (!outfile.is_open())
        {
            std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
            return;
        }
        for (const Vecd& position : observation_theoretical_locations)
        {
            outfile << position[1] << "\n";
        }
        outfile.close();
    }
} // namespace observe_node_cross_section
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, DH));
    water_block_shape.push_back(Vecd(DL + offset_distance, DH));
    water_block_shape.push_back(Vecd(DL + offset_distance, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, 0.0));
    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon computational_domain(createWaterBlockShape());
        add<ExtrudeShape<MultiPolygonShape>>(-offset_distance, computational_domain, "ComputationalDomain");
    }
};

std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge - BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - BW, DH));
    water_block_shape.push_back(Vecd(DL + BW, DH));
    water_block_shape.push_back(Vecd(DL + BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - BW, 0.0));

    return water_block_shape;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
    water_block_shape.push_back(Vecd(DL + 2.0 * BW, DH));
    water_block_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

    return water_block_shape;
}

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_dummy_boundary(createOuterWallShape());
        add<ExtrudeShape<MultiPolygonShape>>(-offset_distance + BW, outer_dummy_boundary, "OuterDummyBoundary");

        MultiPolygon inner_dummy_boundary(createInnerWallShape());
        subtract<ExtrudeShape<MultiPolygonShape>>(-offset_distance, inner_dummy_boundary, "InnerDummyBoundary");
    }
};

inline Real polyEval(const std::vector<Real>& a, Real x)
{
    const int n = static_cast<int>(a.size());
    if (n == 0)
    {
        std::cout << "size of coefficient = 0 " << '\n'
            << "=================\n";
        std::cin.get();
    }
    Real result = a[n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        result = result * x + a[i];
    }
    return result;
}
inline Real imposeByCosineRamp(Real target_value, Real current_time, Real ramp_time)
{
    if (ramp_time <= 0.0)
    {
        return target_value;
    }

    return current_time < ramp_time
        ? 0.5 * target_value * (1.0 - cos(Pi * current_time / ramp_time))
        : target_value;
}
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBox &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_inlet), t_ref_(time_gradually_increase_vel),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Real getVel_table(Vecd& position)
    {
        //** Calculate the distance to wall, Y. position[1] is the distance to the centerline */
        Real Y = half_channel_height - std::abs(position[1]);
        return inletProfileTable(Y);
    }
    
    Real inletProfileTable(Real Y)
    {
        const Real Y_min = PY2_20_VelocityProfile::Y_min;
        const Real Y_max = PY2_20_VelocityProfile::Y_max;
        const Real inv_dY = PY2_20_VelocityProfile::inv_dY;

        constexpr int N =
            static_cast<int>(PY2_20_VelocityProfile::N);

        Y = SMAX(Y_min, SMIN(Y, Y_max));

        Real s = (Y - Y_min) * inv_dY;
        int i = static_cast<int>(s);

        if (i >= N - 1)
            i = N - 2;
        if (i < 0)
            i = 0;

        Real w = s - Real(i);

        return (1.0 - w) * PY2_20_VelocityProfile::U[i]
            + w * PY2_20_VelocityProfile::U[i + 1];
    }

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;

        if (type_velocity_inlet == 3)
        {
            if (inlet_vel_profile_source != "PY2-20")
            {
                std::cout << "Error: inlet velocity profile for table method" << std::endl;
                std::cin.get();
                exit(1);
            }
            target_velocity[0] = imposeByCosineRamp(getVel_table(position), current_time, t_ref_);
        }
        else
        {
            std::cout << "type_inlet_vel: Type wrongly defined! Stop here." << std::endl;
            std::cin.get();
        }

        if (position[1] > half_channel_height)
        {
            std::cout << "Particles out of domain, wrong inlet velocity." << std::endl;
            std::cout << position[1] << std::endl;
            std::cin.get();
        }

        target_velocity[1] = 0.0;
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//	Inflow tke
//----------------------------------------------------------------------
struct InflowTurbulentKineticEnergy
{
    Real u_ref_, t_ref_;
    AlignedBox& aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowTurbulentKineticEnergy(BoundaryConditionType& boundary_condition)
        : u_ref_(U_inlet), t_ref_(time_gradually_increase_vel),
        aligned_box_(boundary_condition.getAlignedBox()),
        halfsize_(aligned_box_.HalfSize()) {
    }

    Real getTKE_table(Vecd& position)
    {
        Real half_channel_height = DH / 2.0;

        //** Calculate the distance to wall, Y. */
        Real Y = half_channel_height - std::abs(position[yAxis]);

        const Real Y_min = PY2_20_TKEProfile::Y_min;
        const Real Y_max = PY2_20_TKEProfile::Y_max;
        const Real inv_dY = PY2_20_TKEProfile::inv_dY;

        constexpr int N =
            static_cast<int>(PY2_20_TKEProfile::N);

        Y = SMAX(Y_min, SMIN(Y, Y_max));

        Real s = (Y - Y_min) * inv_dY;
        int i = static_cast<int>(s);

        if (i >= N - 1)
            i = N - 2;
        if (i < 0)
            i = 0;

        Real w = s - Real(i);

        return (1.0 - w) * PY2_20_TKEProfile::K[i]
            + w * PY2_20_TKEProfile::K[i + 1];
    }

    Real operator()(Vecd& position, Vecd& velocity, Real current_tke, Real current_time)
    {
        Real target_inflow_turbu_k = 0.0;
        if (type_turbulent_inlet_k == 2)
        {
            if (turbulent_inlet_k_profile_source != "PY2-20")
            {
                std::cout << "Error: polynomial turbulent inlet k profile, " << "but the provided data source is: " << turbulent_inlet_k_profile_source << std::endl;
                std::cin.get();
                exit(1);
            }
            target_inflow_turbu_k = getTKE_table(position);
        }
        else
        {
            std::cout << "type_turbulent_inlet_k: Type wrongly defined! Stop here." << std::endl;
            std::cin.get();
        }
        return target_inflow_turbu_k;
    }
};


//----------------------------------------------------------------------
//	Inflow tsdr
//----------------------------------------------------------------------
struct InflowTurbulentSpecificDissipationRate
{
    Real u_ref_, t_ref_;
    AlignedBox& aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowTurbulentSpecificDissipationRate(BoundaryConditionType& boundary_condition)
        : u_ref_(U_inlet), t_ref_(time_gradually_increase_vel),
        aligned_box_(boundary_condition.getAlignedBox()),
        halfsize_(aligned_box_.HalfSize()) {
    }

    Real getTSDR_table(Vecd& position, Real& turbu_omega)
    {
        Real half_channel_height = DH / 2.0;

        //** Calculate the distance to wall, Y. */
        Real Y = half_channel_height - std::abs(position[yAxis]);

        const Real Y_min = PY2_20_OmegaProfile::Y_min;
        const Real Y_max = PY2_20_OmegaProfile::Y_max;
        const Real inv_dY = PY2_20_OmegaProfile::inv_dY;

        constexpr int N =
            static_cast<int>(PY2_20_OmegaProfile::N);

        Y = SMAX(Y_min, SMIN(Y, Y_max));

        Real s = (Y - Y_min) * inv_dY;
        int i = static_cast<int>(s);

        if (i >= N - 1)
            i = N - 2;
        if (i < 0)
            i = 0;

        Real w = s - Real(i);

        return (1.0 - w) * PY2_20_OmegaProfile::Omega[i]
            + w * PY2_20_OmegaProfile::Omega[i + 1];
    }

    Real operator()(Vecd& position, Vecd& velocity, Real current_tsdr, Real current_time, Real current_tke)
    {
        Real target_inflow_turbu_omega = 0.0;
        if (type_turbulent_inlet_omega == 2)
        {
            if (turbulent_inlet_omega_profile_source != "PY2-20")
            {
                std::cout << "Error: polynomial turbulent inlet omega profile, " << "but the provided data source is: " << turbulent_inlet_omega_profile_source << std::endl;
                std::cin.get();
                exit(1);
            }
            target_inflow_turbu_omega = getTSDR_table(position, current_tsdr);
        }
        else
        {
            std::cout << "type_turbulent_inlet_omega: Type wrongly defined! Stop here." << std::endl;
            std::cin.get();
        }
        return target_inflow_turbu_omega;
    }
};

struct RightOutflowPressure
{
    template <class BoundaryConditionType>
    RightOutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

namespace SPH 
{
//=================================================================================================//
    class UpdateVolume : public LocalDynamics
    {
    public:
        explicit UpdateVolume(SPHBody& sph_body)
            : LocalDynamics(sph_body),
            rho_(particles_->getVariableDataByName<Real>("Density")),
            mass_(particles_->getVariableDataByName<Real>("Mass")),
            Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
            indicator_(particles_->getVariableDataByName<int>("Indicator")) {}
        virtual ~UpdateVolume() {};

        void update(size_t index_i, Real dt = 0.0) 
        {
            Vol_[index_i] = mass_[index_i] / rho_[index_i];
        };

    protected:
        Real* rho_, * mass_, * Vol_;
        int* indicator_;
    };
//=================================================================================================//
}
