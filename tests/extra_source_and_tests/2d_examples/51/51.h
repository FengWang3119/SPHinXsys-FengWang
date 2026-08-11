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
Real num_fluid_cross_section = 20.0;

constexpr Real wave_amplitude = 0.1;
constexpr Real wave_length = 1.0;
constexpr Real pi = 3.14159265358979323846;
/**
 * @brief Lower sinusoidal wall
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

Vecd external_acc = Vecd(0.01687141, 0.0);
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

//Real y_p_constant_sublayer = 0.0005;

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

Real BW = resolution_ref * 4; /**< Reference size of the emitter. */
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBoxd system_domain_bounds(Vecd(-2.0 * BW, -wave_amplitude - BW), Vecd(DL + 2.0 * BW, DH + BW));

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real U_inlet = 0.816;
Real U_f = U_inlet;         //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */

Real mu_f = 1.0e-4;

Real Re_calculated = U_f * DH * rho0_f / mu_f;

//----------------------------------------------------------------------
//  Center-point observer
//----------------------------------------------------------------------

Real x_observe_center = 0.5 * DL;
Real y_observe_center =
0.5 * (DH + lowerWallHeight(x_observe_center));

StdVec<Vecd> observer_location_center_point = {
    Vecd(x_observe_center, y_observe_center) };

//----------------------------------------------------------------------
//  Cross-section observers at the wave trough and crest
//----------------------------------------------------------------------
namespace observe_cross_sections
{
    constexpr const char* namespace_prefix = "CrossSection_";
    constexpr int number_of_cross_sections = 2;

    // Based on:
    // y_w(x) = -A sin(2 pi x / lambda)
    //
    // x = 0.25 lambda: trough
    // x = 0.75 lambda: crest
    Real observe_x[number_of_cross_sections] = {
        0.25 * wave_length,
        0.75 * wave_length };

    const char* section_name[number_of_cross_sections] = {
        "trough",
        "crest" };

    // Preserve the original offset-observer treatment.
    // Set this value to zero if the offset is no longer needed.
    Real observer_offset_distance = 2.0 * resolution_ref;

    int num_observer_points[number_of_cross_sections] = { 0, 0 };
    Real observe_spacing[number_of_cross_sections] = { 0.0, 0.0 };

    StdVec<Vecd> observation_locations;
    StdVec<Vecd> observation_theoretical_locations;

    /**
     * @brief Generate observers at the trough and crest cross-sections.
     */
    void getObservationLocations()
    {
        observation_locations.clear();
        observation_theoretical_locations.clear();

        for (int k = 0; k < number_of_cross_sections; ++k)
        {
            Real x = observe_x[k];
            Real y_bottom = lowerWallHeight(x);
            Real local_channel_height = DH - y_bottom;

            num_observer_points[k] = std::max(1, static_cast<int>(std::round(local_channel_height / resolution_ref)));

            observe_spacing[k] = local_channel_height / static_cast<Real>(num_observer_points[k]);

            for (int i = 0; i < num_observer_points[k]; ++i)
            {
                // Cell-centered theoretical observation position.
                Real y_theoretical =
                    y_bottom +
                    (static_cast<Real>(i) + 0.5) *
                    observe_spacing[k];

                Vecd theoretical_position(x, y_theoretical);
                Vecd actual_position = theoretical_position;

                // Preserve the original offset treatment:
                // move the first point into the lower wall and
                // the last point into the upper wall.
                if (i == 0)
                {
                    actual_position[1] -= observer_offset_distance;
                }
                else if (i == num_observer_points[k] - 1)
                {
                    actual_position[1] += observer_offset_distance;
                }

                observation_locations.push_back(actual_position);
                observation_theoretical_locations.push_back(
                    theoretical_position);
            }
        }
    }

    /**
     * @brief Output actual observer positions.
     *
     * Columns:
     * section index, x, y
     */
    void outputObservePositions()
    {
        std::string filename =
            "../bin/output/" +
            std::string(namespace_prefix) +
            "observer_positions.dat";

        std::ofstream outfile(filename);

        if (!outfile.is_open())
        {
            std::cerr
                << "Error: Unable to open file "
                << filename
                << " for writing."
                << std::endl;
            return;
        }

        size_t global_index = 0;

        for (int k = 0; k < number_of_cross_sections; ++k)
        {
            for (int i = 0; i < num_observer_points[k]; ++i)
            {
                const Vecd& position =
                    observation_locations[global_index];

                outfile
                    << k << " "
                    << position[0] << " "
                    << position[1] << "\n";

                ++global_index;
            }
        }
    }

    /**
     * @brief Output theoretical y-coordinates separately
     *        for the trough and crest.
     */
    void outputTheoreticalY()
    {
        size_t global_index = 0;

        for (int k = 0; k < number_of_cross_sections; ++k)
        {
            std::string filename =
                "../bin/output/" +
                std::string(namespace_prefix) +
                section_name[k] +
                "_theoretical_y.dat";

            std::ofstream outfile(filename);

            if (!outfile.is_open())
            {
                std::cerr
                    << "Error: Unable to open file "
                    << filename
                    << " for writing."
                    << std::endl;
                return;
            }

            for (int i = 0; i < num_observer_points[k]; ++i)
            {
                outfile
                    << observation_theoretical_locations
                    [global_index][1]
                    << "\n";

                ++global_index;
            }
        }
    }

    /**
     * @brief Output the number of observers
     *        on each cross-section.
     */
    void outputNumberOfObserverPoints()
    {
        std::string filename =
            "../bin/output/" +
            std::string(namespace_prefix) +
            "num_points.dat";

        std::ofstream outfile(filename);

        if (!outfile.is_open())
        {
            std::cerr
                << "Error: Unable to open file "
                << filename
                << " for writing."
                << std::endl;
            return;
        }

        for (int k = 0; k < number_of_cross_sections; ++k)
        {
            outfile
                << section_name[k] << " "
                << num_observer_points[k] << "\n";
        }
    }

} // namespace observe_cross_sections

//----------------------------------------------------------------------
//  Closest SPH-node observers at the trough and crest
//----------------------------------------------------------------------
namespace observe_node_cross_sections
{
    constexpr const char* namespace_prefix = "Node_";
    constexpr int number_of_cross_sections = 2;

    // For y_w(x) = -A sin(2 pi x / lambda):
    //
    // x = 0.25 lambda: wave trough
    // x = 0.75 lambda: wave crest
    Real observe_x[number_of_cross_sections] = {
        0.25 * wave_length,
        0.75 * wave_length };

    const char* section_name[number_of_cross_sections] = {
        "trough",
        "crest" };

    // Same offset treatment as in the original straight-channel case.
    Real observer_offset_distance =
        2.0 * resolution_ref;

    StdVec<Vecd> observation_locations;
    StdVec<Vecd> observation_theoretical_locations;

    void getObservationLocations()
    {
        observation_locations.clear();
        observation_theoretical_locations.clear();

        for (int k = 0; k < number_of_cross_sections; ++k)
        {
            Real x = observe_x[k];
            Real y_bottom = lowerWallHeight(x);

            // Theoretical locations of the closest SPH points
            // adjacent to the lower and upper walls.
            Real lower_node_y =
                y_bottom + 0.5 * resolution_ref;

            Real upper_node_y =
                DH - 0.5 * resolution_ref;

            Vecd lower_node_theoretical(
                x,
                lower_node_y);

            Vecd upper_node_theoretical(
                x,
                upper_node_y);

            // Actual observer positions are shifted into the wall,
            // following exactly the original offset treatment.
            Vecd lower_observer_position(
                x,
                lower_node_y -
                observer_offset_distance);

            Vecd upper_observer_position(
                x,
                upper_node_y +
                observer_offset_distance);

            observation_locations.push_back(
                lower_observer_position);

            observation_locations.push_back(
                upper_observer_position);

            observation_theoretical_locations.push_back(
                lower_node_theoretical);

            observation_theoretical_locations.push_back(
                upper_node_theoretical);
        }
    }

    void outputObservePositions()
    {
        std::string filename =
            "../bin/output/" +
            std::string(namespace_prefix) +
            "observer_positions.dat";

        std::ofstream outfile(filename);

        if (!outfile.is_open())
        {
            std::cerr
                << "Error: Unable to open file "
                << filename
                << " for writing."
                << std::endl;
            return;
        }

        for (int k = 0; k < number_of_cross_sections; ++k)
        {
            size_t lower_index =
                static_cast<size_t>(2 * k);

            size_t upper_index =
                lower_index + 1;

            outfile
                << section_name[k] << " lower "
                << observation_locations[lower_index][0] << " "
                << observation_locations[lower_index][1]
                << "\n";

            outfile
                << section_name[k] << " upper "
                << observation_locations[upper_index][0] << " "
                << observation_locations[upper_index][1]
                << "\n";
        }
    }

    void outputTheoreticalPositions()
    {
        std::string filename =
            "../bin/output/" +
            std::string(namespace_prefix) +
            "observer_theoretical_positions.dat";

        std::ofstream outfile(filename);

        if (!outfile.is_open())
        {
            std::cerr
                << "Error: Unable to open file "
                << filename
                << " for writing."
                << std::endl;
            return;
        }

        for (int k = 0; k < number_of_cross_sections; ++k)
        {
            size_t lower_index =
                static_cast<size_t>(2 * k);

            size_t upper_index =
                lower_index + 1;

            outfile
                << section_name[k] << " lower "
                << observation_theoretical_locations
                [lower_index][0]
                << " "
                << observation_theoretical_locations
                [lower_index][1]
                << "\n";

            outfile
                << section_name[k] << " upper "
                << observation_theoretical_locations
                [upper_index][0]
                << " "
                << observation_theoretical_locations
                [upper_index][1]
                << "\n";
        }
    }

} // namespace observe_node_cross_sections

//----------------------------------------------------------------------
//  Case-dependent geometry: periodic wavy channel
//----------------------------------------------------------------------
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
