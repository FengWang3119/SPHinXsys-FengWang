/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D. before is 1
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
//#include "k-epsilon_turbulent_model.cpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
#include "zeroth-order_residue.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
//** Dimension: m s kg */
//** Z along the channel length */
Real scale = 1.0e-3;
Real D_thr = 4.0 * scale;
Real DH = 3.0 * D_thr; /**< Channel height. */

Real Radius_inlet = DH / 2.0;

Real DL = 50.0 * D_thr; /**< Total domain length, note that this is determined by the INPUT geometry. */

Vecd point_O(0.0, 0.0, 0.0); //** O A are the start and end point of the total domain */
Vecd point_A = point_O + Vecd(0.0, 0.0, DL);
Vecd point_OA_half = (point_O + point_A) / 2.0;

Real num_fluid_cross_section = 30.0;                //** On inlet */
Real resolution_ref = DH / num_fluid_cross_section; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                       /**< Reference size of the emitter. */
Real half_channel_height = DH / 2.0;
Real buffer_thickness = 5.0 * resolution_ref;
//Real DL_sponge = buffer_thickness; //** Note, no pre-left sponge since the geometry is already extended, and the quantitative data is obtained from the throat */

const int SimTK_resolution = 20;
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(point_O + Vecd(-Radius_inlet, -Radius_inlet, 0.0) + 2.0 * Vecd(-BW, -BW, -BW),
                                 point_A + Vecd(Radius_inlet, Radius_inlet, 0.0) + 2.0 * Vecd(BW, BW, BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real U_inlet = 0.046;
Real U_f = U_inlet; //*Characteristic velocity
Real U_thr = (DH / D_thr) * (DH / D_thr) * U_inlet;
Real U_max = 2.0 * U_thr; //** An estimated value, generally 2.0 U_inlet for 3D *

Real c_f = 10.0 * U_max;
Real rho0_f = 1056.0; /**< Density. */
Real Re = 500.0;

Real mu_f = rho0_f * U_thr * D_thr / Re;

Real start_up_time_ref = 0.05;

Real Outlet_pressure = 0.0;
//----------------------------------------------------------------------
//	The open boundary setting.
//----------------------------------------------------------------------
Vecd rotation_axis(1.0, 0.0, 0.0);
Vecd left_buffer_halfsize = Vecd(Radius_inlet, Radius_inlet, 0.5 * buffer_thickness);
Vecd left_buffer_translation = point_O + Vecd(0.0, 0.0, 0.5 * buffer_thickness);
Vecd right_buffer_halfsize = Vecd(Radius_inlet, Radius_inlet, 0.5 * buffer_thickness);
Vecd right_buffer_translation = point_A + Vecd(0.0, 0.0, -0.5 * buffer_thickness);
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
/** Set the file path to the stl file. */
std::string stl_fluid_path_volume = "./input/fda_nozzle.stl";
Real scale_factor = 1.0e-3;
Vecd translation_stl(0.0, 0.0, 0.0);
class WaterBlockFromSTL : public ComplexShape
{
  public:
    explicit WaterBlockFromSTL(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stl_fluid_path_volume, translation_stl, scale_factor);
    }
};

/** Set the file path to the stl file. */
std::string stl_wall_path_volume = "./input/fda_nozzle.stl";
Real scale_factor_wall = 1.0e-3;
Vecd translation_stl_wall(0.0, 0.0, 0.0);
class WallBoundaryFromSTL : public ComplexShape
{
  public:
    explicit WallBoundaryFromSTL(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(BW, stl_wall_path_volume, translation_stl_wall, scale_factor_wall);
        subtract<TriangleMeshShapeSTL>(stl_wall_path_volume, translation_stl_wall, scale_factor_wall);
        subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0.0, 0.0, 1.0), Radius_inlet,
                                            (1.8 * BW) * 0.5, SimTK_resolution,
                                            point_O);
        subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0.0, 0.0, 1.0), Radius_inlet,
                                            (1.8 * BW) * 0.5, SimTK_resolution,
                                            point_A);
    }
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_inlet), t_ref_(start_up_time_ref),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        //** 3D modification */
        Real local_radius_square = position[0] * position[0] + position[1] * position[1];
        Real Radius_inlet_square = Radius_inlet * Radius_inlet;
        target_velocity[2] = SMAX(2.0 * u_ave * (1.0 - local_radius_square / Radius_inlet_square), 0.0);
        //target_velocity[2] = u_ave;

        // if (local_radius_square > Radius_inlet_square)
        // {
        //     std::cout << "Particles out of domain, wrong inlet velocity." << std::endl;
        //     std::cout << "local_radius_square=" << local_radius_square << std::endl;
        //     std::cout << "Radius_inlet=" << Radius_inlet_square << std::endl;
        //     //std::this_thread::sleep_for(std::chrono::seconds(10));
        // }
        target_velocity[0] = 0.0;
        target_velocity[1] = 0.0;
        return target_velocity;
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
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
//** For getting centerline velocity *
namespace observe_centerline
{
constexpr const char *namespace_prefix = "centerline";
Vecd pos_observe_start = point_O;
Real sparsity_ratio = 5.0;
Real length_observing_line = DL;
Vecd unit_direction_observe = Vecd(0.0, 0.0, 1.0);
Real observer_offset_distance = 0.0;

int num_observer_points = std::round(length_observing_line / resolution_ref / sparsity_ratio); //**Every particle is regarded as a cell monitor*
Real observe_spacing = length_observing_line / num_observer_points;

StdVec<Vecd> observation_location;
void get_observation_locations()
{
    for (int i = 0; i < num_observer_points; ++i)
    {
        Vecd pos_observer_i = pos_observe_start + i * observe_spacing * unit_direction_observe;
        if (i == 0)
        {
            pos_observer_i -= observer_offset_distance * unit_direction_observe;
        }
        if (i == num_observer_points - 1)
        {
            pos_observer_i += observer_offset_distance * unit_direction_observe;
        }
        observation_location.push_back(pos_observer_i);
    }
}
void output_observer_theoretical_pos_on_line()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_theoretical_pos_on_line.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (int i = 0; i < num_observer_points; ++i)
    {
        outfile << observation_location[i].dot(unit_direction_observe) << "\n";
    }
    outfile.close();
}
} // namespace observe_centerline

//** For getting cross-section velocity *
namespace observe_cross_sections
{
Real observe_base_z = 25.0 * D_thr;
constexpr const char *namespace_prefix = "cross_sections";
const int number_observe_line = 5;
Real observer_offset_distance = 2.0 * resolution_ref;
Vecd unit_direction_observe(0.0, 1.0, 0.0);
// ** Determine the observing start point. *
Real observe_start_z[number_observe_line] = {
    observe_base_z - 0.088,
    observe_base_z - 0.008,
    observe_base_z + 0.008,
    observe_base_z + 0.024,
    observe_base_z + 0.032};

Real observe_start_x[number_observe_line] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0};

Real observe_start_y[number_observe_line] = {
    0.5 * resolution_ref - Radius_inlet,
    0.5 * resolution_ref - Radius_inlet,
    0.5 * resolution_ref - Radius_inlet,
    0.5 * resolution_ref - Radius_inlet,
    0.5 * resolution_ref - Radius_inlet};

// ** Determine the length of the observing line and other information. *
Real observe_line_length[number_observe_line] = {0.0};
int num_observer_points[number_observe_line] = {0};

void getObservingLineLengthAndEndPoints()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        observe_line_length[i] = DH;
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
        Vecd pos_observe_start(observe_start_x[k], observe_start_y[k], observe_start_z[k]);
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
        for (int i = 0; i < position.size(); ++i)
        {
            outfile << position[i] << " ";
        }
        outfile << "\n";
    }
    outfile.close();
}
void output_observer_theoretical_pos_on_line()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_theoretical_pos_on_line.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (int j = 0; j < number_observe_line; ++j)
    {
        for (int i = 0; i < num_observer_points[j]; ++i)
        {
            outfile << observation_theoretical_locations[i].dot(unit_direction_observe) << "\n";
        }
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
//** For getting cross-section velocity on Y=0 *
namespace observe_cross_sections_Y0
{
Real observe_base_z = 25.0 * D_thr;
constexpr const char *namespace_prefix = "cross_sections_Y0";
const int number_observe_line = 5;
Real observer_offset_distance = 2.0 * resolution_ref;
Vecd unit_direction_observe(1.0, 0.0, 0.0);
// ** Determine the observing start point. *
Real observe_start_z[number_observe_line] = {
    observe_base_z - 0.064,
    observe_base_z - 0.008,
    observe_base_z + 0.008,
    observe_base_z + 0.024,
    observe_base_z + 0.08};

Real observe_start_y[number_observe_line] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0};

Real observe_start_x[number_observe_line] = {
    0.5 * resolution_ref - Radius_inlet,
    0.5 * resolution_ref - Radius_inlet,
    0.5 * resolution_ref - Radius_inlet,
    0.5 * resolution_ref - Radius_inlet,
    0.5 * resolution_ref - Radius_inlet};

// ** Determine the length of the observing line and other information. *
Real observe_line_length[number_observe_line] = {0.0};
int num_observer_points[number_observe_line] = {0};

void getObservingLineLengthAndEndPoints()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        observe_line_length[i] = DH;
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
        Vecd pos_observe_start(observe_start_x[k], observe_start_y[k], observe_start_z[k]);
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
        for (int i = 0; i < position.size(); ++i)
        {
            outfile << position[i] << " ";
        }
        outfile << "\n";
    }
    outfile.close();
}
void output_observer_theoretical_pos_on_line()
{
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_theoretical_pos_on_line.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (int j = 0; j < number_observe_line; ++j)
    {
        for (int i = 0; i < num_observer_points[j]; ++i)
        {
            outfile << observation_theoretical_locations[i].dot(unit_direction_observe) << "\n";
        }
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
} // namespace observe_cross_sections_Y0
//** For regression test *
StdVec<Vecd> observer_location_center_point = {point_O + Vecd(0.0, 0.0, 0.5 * DL)};