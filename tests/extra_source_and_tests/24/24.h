#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 2.0; /**< Channel height. */
Real half_channel_height = DH / 2.0;
Real characteristic_length = DH;
Real num_fluid_cross_section = 20.0;
Real resolution_ref = DH / num_fluid_cross_section; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                       /**< Extending width for BCs. */
Real buffer_thickness = 5.0 * resolution_ref;

Real extend_in = 0.0;
Real extend_out = 0.0;
Real extend_compensate_relaxation = 0.0;

// Real DH1 = 5.0 * DH;
// Real DL2 = 10.0 * DH;
Real DH1 = 20.0 * DH;
Real DL2 = 50.0 * DH;

Vec2d point_O(0.0, 0.0);
Vec2d point_A = point_O + Vec2d(0.0, DH + 2.0 * DH1);
Vec2d point_B = point_A + Vec2d(DL2, 0.0);
Vec2d point_C = point_O + Vec2d(DL2, 0.0);

Vec2d point_OA_half = (point_O + point_A) / 2.0;
Vec2d point_BC_half = (point_B + point_C) / 2.0;

Vec2d point_D = point_O + Vec2d(0.0, DH1);
Vec2d point_E = point_D + Vec2d(0.0, DH);
Vec2d point_F = point_E + Vec2d(buffer_thickness, 0.0);
Vec2d point_G = point_D + Vec2d(buffer_thickness, 0.0);
Vec2d point_H = point_A + Vec2d(buffer_thickness, 0.0);
Vec2d point_I = point_O + Vec2d(buffer_thickness, 0.0);

Vec2d point_AO_half = (point_A + point_O) / 2.0;
Vec2d point_AE_half = (point_A + point_E) / 2.0;
Vec2d point_ED_half = (point_E + point_D) / 2.0;
Vec2d point_DO_half = (point_D + point_O) / 2.0;

Vec2d point_AB_half = (point_A + point_B) / 2.0;

Vec2d point_J = point_H + Vec2d(0.0, -buffer_thickness);
Vec2d point_K = point_B + Vec2d(-buffer_thickness, -buffer_thickness);

Vec2d point_OC_half = (point_O + point_C) / 2.0;

BoundingBox system_domain_bounds(Vec2d(point_O[0], point_O[1]) + Vec2d(-BW, -BW), point_B + Vec2d(BW, BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real U_inlet = 1.0;
Real U_f = U_inlet;         //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real T_ref = 2.0;
Real rho0_f = 1.0;
Real Re = 40.0;

Real Outlet_pressure = 0.0;
Real Freestream_pressure = 0.0;

Real mu_f = rho0_f * U_f * characteristic_length / Re;

Real Re_calculated = U_f * characteristic_length * rho0_f / mu_f;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
// Vec2d left_bidirectional_buffer_halfsize = 0.5 * Vec2d(buffer_thickness, (point_A[1] - point_O[1]));
// Vec2d left_bidirectional_translation = point_AO_half + Vec2d(0.5 * buffer_thickness, 0.0);
Vec2d left_bidirectional_buffer_halfsize = 0.5 * Vec2d(buffer_thickness, (point_E[1] - point_D[1]));
Vec2d left_bidirectional_translation = point_ED_half + Vec2d(0.5 * buffer_thickness, 0.0);

Vec2d right_bidirectional_buffer_halfsize = 0.5 * Vec2d(buffer_thickness, (point_B[1] - point_C[1]));
Vec2d right_bidirectional_translation = point_BC_half - Vec2d(0.5 * buffer_thickness, 0.0);
Vec2d normal = Vec2d(1.0, 0.0);

Vec2d static_buffer_halfsize_up = 0.5 * Vec2d(buffer_thickness, (point_A[1] - point_E[1]));
Vec2d static_translation_up = point_AE_half + Vec2d(0.5 * buffer_thickness, 0.0);
Vec2d static_buffer_halfsize_down = 0.5 * Vec2d(buffer_thickness, (point_D[1] - point_O[1]));
Vec2d static_translation_down = point_DO_half + Vec2d(0.5 * buffer_thickness, 0.0);

Vec2d up_buffer_halfsize = 0.5 * Vec2d(buffer_thickness, (point_K[0] - point_J[0]));
Vec2d up_buffer_translation = point_AB_half + Vec2d(0.0, -0.5 * buffer_thickness);
Real up_buffer_rotation_angle = -0.5 * Pi; //** Negative means clock-wise */

Vec2d down_buffer_halfsize = 0.5 * Vec2d(buffer_thickness, (point_K[0] - point_J[0]));
Vec2d down_buffer_translation = point_OC_half + Vec2d(0.0, 0.5 * buffer_thickness);
Real down_buffer_rotation_angle = 0.5 * Pi; //** Negative means clock-wise */

//----------------------------------------------------------------------
// Observation for regression test.
//----------------------------------------------------------------------
StdVec<Vecd> observer_location = {Vecd(0.5 * DL2, 0.5 * DH)}; /**< Displacement observation point. */
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
// ** By kernel weight. *
//** For getting centerline velocity *
namespace observe_centerline
{
Real offset_distance = 0.0;
const int number_observe_line = 1;
Real observer_offset_distance = 0.0 * resolution_ref;
Vec2d unit_direction_observe(1.0, 0.0);
// ** Determine the observing start point of the each line. *
Real observe_start_x[number_observe_line] = {point_G[0]};
Real observe_start_y[number_observe_line] = {point_ED_half[1]};
// ** Determine the length of the observing line and other information. *
Real observe_line_length[number_observe_line] = {0.0};
int num_observer_points[number_observe_line] = {0};
void getObservingLineLengthAndEndPoints()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        observe_line_length[i] = point_K[0] - point_G[0];
        num_observer_points[i] = std::round(observe_line_length[i] / resolution_ref);
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
    std::string filename = "../bin/output/observer_positions.dat";
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
    std::string filename = "../bin/output/observer_theoretical_x.dat";
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
    std::string filename = "../bin/output/observer_num_points_on_lines.dat";
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
} // namespace observe_centerline
//** For getting cross-section velocity *
namespace observe_cross_sections
{
constexpr const char *namespace_prefix = "cross_sections";
const int number_observe_line = 50;
Real observer_offset_distance = 0.0;
Vec2d unit_direction_observe(0.0, 1.0);
// ** Determine the observing start point. *
Real observe_start_x[number_observe_line] = {0.0};
Real observe_start_y[number_observe_line] = {0.0};
Real observe_start_x_spacing = (point_K[0] - point_J[0]) / Real(number_observe_line - 1); //** 3 lines means actually 2 spaceing */

void get_observe_start_coordinate()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        observe_start_x[i] = point_J[0] + i * observe_start_x_spacing;
        observe_start_y[i] = point_OA_half[1];
    }
}

// ** Determine the length of the observing line and other information. *
Real observe_line_length[number_observe_line] = {0.0};
int num_observer_points[number_observe_line] = {0};

void getObservingLineLengthAndEndPoints()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        observe_line_length[i] = point_J[1] - point_ED_half[1];
        num_observer_points[i] = std::round(observe_line_length[i] / resolution_ref);
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
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_observer_positions.dat";
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
void output_observe_theoretical_y()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_observer_theoretical_y.dat";
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
void output_observe_line_pos_x()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_observer_line_pos_x.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (const Real &position : observe_start_x)
    {
        outfile << position << "\n";
    }
    outfile.close();
}
void output_number_observe_points_on_lines()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_observer_num_points_on_lines.dat";
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
} // namespace observe_cross_sections
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

struct FreestreamPressure
{
    template <class BoundaryConditionType>
    FreestreamPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        /*constant pressure*/
        Real pressure = Freestream_pressure;
        return pressure;
    }
};
//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ave;
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ave(0.0), u_ref_(U_inlet), t_ref_(T_ref),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();

        u_ave = current_time < t_ref_ ? 0.5 * U_inlet * (1.0 - cos(Pi * current_time / t_ref_)) : U_inlet;

        //target_velocity[0] = u_ave;
        //target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / half_channel_height / half_channel_height);
        if (abs(position[1]) <= half_channel_height)
        {
            target_velocity[0] = u_ave;
            target_velocity[1] = 0.0;
        }

        return target_velocity;
    }
};
struct StaticInflowVelocity
{

    template <class BoundaryConditionType>
    StaticInflowVelocity(BoundaryConditionType &boundary_condition) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(point_O);
        water_block_shape.push_back(point_A);
        water_block_shape.push_back(point_B);
        water_block_shape.push_back(point_C);
        water_block_shape.push_back(point_O);
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
// class WallBoundary : public MultiPolygonShape
// {
//   public:
//     explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
//     {
// std::vector<Vecd> outer_wall_shape;
// outer_wall_shape.push_back(point_O + Vecd(0.0, -BW)); //** Keep the section neat */
// outer_wall_shape.push_back(point_A + Vecd(0.0, +BW));
// outer_wall_shape.push_back(point_B + Vecd(0.0, +BW));
// outer_wall_shape.push_back(point_C + Vecd(0.0, -BW));
// outer_wall_shape.push_back(point_O + Vecd(0.0, -BW));

// std::vector<Vecd> inner_wall_shape;
// inner_wall_shape.push_back(point_O + Vecd(-BW, 0.0));
// inner_wall_shape.push_back(point_A + Vecd(-BW, 0.0));
// inner_wall_shape.push_back(point_B + Vecd(+BW, 0.0));
// inner_wall_shape.push_back(point_C + Vecd(+BW, 0.0));
// inner_wall_shape.push_back(point_O + Vecd(-BW, 0.0));

// multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
// multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);

//         std::vector<Vecd> wall_shape;
//         wall_shape.push_back(point_O + Vecd(0.0, -BW));
//         wall_shape.push_back(point_C + Vecd(0.0, -BW));
//         wall_shape.push_back(point_C + Vecd(0.0, -BW) + Vecd(0.0, -BW));
//         wall_shape.push_back(point_O + Vecd(0.0, -BW) + Vecd(0.0, -BW));
//         wall_shape.push_back(point_O + Vecd(0.0, -BW));

//         multi_polygon_.addAPolygon(wall_shape, ShapeBooleanOps::add);
//     }
// };

template <int INDICATOR>
class IndicatedParticlesExcludeInletBuffer : public WithinScope
{
    int *indicator_;
    int *buffer_indicator_;
    Vecd *pos_;

  public:
    explicit IndicatedParticlesExcludeInletBuffer(BaseParticles *base_particles)
        : WithinScope(),
          indicator_(base_particles->getVariableDataByName<int>("Indicator")),
          buffer_indicator_(base_particles->getVariableDataByName<int>("BufferParticleIndicator")),
          pos_(base_particles->getVariableDataByName<Vecd>("Position")){};
    bool operator()(size_t index_i)
    {
        if (pos_[index_i][0] < point_F[0] && pos_[index_i][1] > point_F[1]) //** Exclude static inlet buffer up */
        {
            return false;
        }
        if (pos_[index_i][0] < point_G[0] && pos_[index_i][1] < point_G[1]) //** Exclude static inlet buffer down */
        {
            return false;
        }
        if (indicator_[index_i] == INDICATOR)
        {
            return true;
        }
        else
        {
            return false;
        }
    };
};

using BulkParticlesWithoutInlet = IndicatedParticlesExcludeInletBuffer<0>;