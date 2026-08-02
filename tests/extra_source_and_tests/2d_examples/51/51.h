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
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 1.0;  /**< Channel height. */
Real DL = 1.0; /**< Channel length. */
Real num_fluid_cross_section = 40.0;

constexpr Real wave_amplitude = 0.1;
constexpr Real wave_length = 1.0;
constexpr Real pi = 3.14159265358979323846;

Vecd external_acc = Vecd(0.00258, 0.0);
Real external_acc_gradually_impose_t = 2.0;
//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------
//** Tag for wall treatment *
int is_blended = 0;
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

Real y_p_constant_sublayer = 0.0005;

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
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBoxd system_domain_bounds(Vecd(-2.0 * BW, -wave_amplitude - BW), Vecd(DL + 2.0 * BW, DH + BW));

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

//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
Real x_observe_start = 0.99 * DL;
int num_observer_points = std::round(DH / resolution_ref); //**Every particle is regarded as a cell monitor*
Real observe_spacing = DH / num_observer_points;

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
//  Case-dependent geometry: periodic wavy channel
//----------------------------------------------------------------------
/**
 * @brief Lower sinusoidal wall.
 *
 * x = 0.00 lambda: y =  0
 * x = 0.25 lambda: y = -A, wave trough
 * x = 0.75 lambda: y =  A, wave crest
 * x = 1.00 lambda: y =  0
 */
Real lowerWallHeight(Real x)
{
    return -wave_amplitude *
        std::sin(2.0 * pi * x / wave_length);
}

/**
 * @brief Fluid domain bounded by a flat upper wall and a wavy lower wall.
 */
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;

    // Number of straight segments used to approximate the sinusoidal wall.
    constexpr size_t number_of_segments = 200;

    // Left and upper boundaries.
    water_block_shape.push_back(Vecd(0.0, lowerWallHeight(0.0)));
    water_block_shape.push_back(Vecd(0.0, DH));
    water_block_shape.push_back(Vecd(DL, DH));

    // Lower wavy boundary, traversed from right to left.
    for (size_t i = 0; i <= number_of_segments; ++i)
    {
        Real x = DL *
            static_cast<Real>(number_of_segments - i) /
            static_cast<Real>(number_of_segments);

        water_block_shape.push_back(Vecd(x, lowerWallHeight(x)));
    }

    return water_block_shape;
}

class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name)
        : ComplexShape(shape_name)
    {
        MultiPolygon computational_domain(createWaterBlockShape());
        add<MultiPolygonShape>(
            computational_domain, "ComputationalDomain");
    }
};

/**
 * @brief Flat upper wall.
 *
 * The wall is extended in the streamwise direction to provide
 * sufficient wall particles near the periodic boundaries.
 */
std::vector<Vecd> createUpperWallShape()
{
    Real wall_x_min = -2.0 * BW;
    Real wall_x_max = DL + 2.0 * BW;

    std::vector<Vecd> upper_wall_shape;

    upper_wall_shape.push_back(Vecd(wall_x_min, DH));
    upper_wall_shape.push_back(Vecd(wall_x_min, DH + BW));
    upper_wall_shape.push_back(Vecd(wall_x_max, DH + BW));
    upper_wall_shape.push_back(Vecd(wall_x_max, DH));
    upper_wall_shape.push_back(Vecd(wall_x_min, DH));

    return upper_wall_shape;
}

/**
 * @brief Sinusoidal lower wall.
 */
std::vector<Vecd> createLowerWallShape()
{
    std::vector<Vecd> lower_wall_shape;

    Real wall_x_min = -2.0 * BW;
    Real wall_x_max = DL + 2.0 * BW;

    constexpr size_t number_of_segments = 240;

    // Fluid-facing surface, from left to right.
    for (size_t i = 0; i <= number_of_segments; ++i)
    {
        Real x = wall_x_min +
            (wall_x_max - wall_x_min) *
            static_cast<Real>(i) /
            static_cast<Real>(number_of_segments);

        lower_wall_shape.push_back(
            Vecd(x, lowerWallHeight(x)));
    }

    // Outer surface, from right to left.
    for (size_t i = 0; i <= number_of_segments; ++i)
    {
        Real x = wall_x_max -
            (wall_x_max - wall_x_min) *
            static_cast<Real>(i) /
            static_cast<Real>(number_of_segments);

        lower_wall_shape.push_back(
            Vecd(x, lowerWallHeight(x) - BW));
    }

    // Explicitly close the polygon.
    lower_wall_shape.push_back(
        Vecd(wall_x_min, lowerWallHeight(wall_x_min)));

    return lower_wall_shape;
}

/**
 * @brief Wall boundary body definition.
 */
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string& shape_name)
        : ComplexShape(shape_name)
    {
        MultiPolygon upper_wall(createUpperWallShape());
        add<MultiPolygonShape>(upper_wall, "UpperWall");

        MultiPolygon lower_wall(createLowerWallShape());
        add<MultiPolygonShape>(lower_wall, "LowerWavyWall");
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
