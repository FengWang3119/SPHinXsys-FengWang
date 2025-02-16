/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "k-epsilon_turbulent_model.cpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
#include "zeroth-order_residue.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 2.0;   /**< Channel height. */
Real DL = 120.0; /**< Channel length. */
Real num_fluid_cross_section = 20.0;

Real characteristic_length = DH; /**<It needs characteristic Length to determine Re and calculate turbulent length and the inflow turbulent epsilon>*/
//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------
//** For K and Epsilon, type of the turbulent inlet, 0 is freestream, 1 is from interpolation from PY21 *
int type_turbulent_inlet = 1;
Real relaxation_rate_turbulent_inlet = 0.8;
//** Tag for AMRD *
int is_AMRD = 0;
bool is_constrain_normal_velocity_in_P_region = false;
//** Weight for correcting the velocity  gradient in the sub near wall region  *
Real weight_vel_grad_sub_nearwall = 0.1;
//** Tag for Source Term Linearisation *
bool is_source_term_linearisation = false;
//** Initial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = {0.000180001, 3.326679e-5, 1.0e-9};
//----------------------------------------------------------------------
//	Resolution for turbulence.
//----------------------------------------------------------------------
Real y_p_constant = 0.05;
//Real y_p_constant = DH/num_fluid_cross_section; //** For first try */

Real resolution_ref = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
Real offset_distance = y_p_constant - resolution_ref / 2.0;                        //** Basically offset distance is large than or equal to 0 *
Real DH_C = DH - 2.0 * offset_distance;
//----------------------------------------------------------------------
//	Other parameters.
//----------------------------------------------------------------------
Real BW = resolution_ref * 4;
Real DL_sponge = resolution_ref * 20;
Real half_channel_height = DH / 2.0;
Real buffer_thickness = 5.0 * resolution_ref;
Vec2d buffer_normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + 2.0 * BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real U_inlet = 1.0;
Real U_f = U_inlet;         //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */
Real Re = 20000.0;

Real Outlet_pressure = 0.0;

Real mu_f = rho0_f * U_f * characteristic_length / Re;

Real Re_calculated = U_f * characteristic_length * rho0_f / mu_f;
//----------------------------------------------------------------------
//	The buffer block with offset model.
//----------------------------------------------------------------------
//** Note that even if offset(reduce) the domain, the buffer region still cover the domain*/
//** the horizontal boundary offset loss will be compensated specifically*/
Vec2d left_buffer_halfsize = 0.5 * Vec2d(buffer_thickness, DH);
Vec2d left_buffer_translation = left_buffer_halfsize + Vec2d(-DL_sponge, 0.0); //** old version, the latest one not use Sponge */

Vec2d right_buffer_halfsize = 0.5 * Vec2d(buffer_thickness, DH);
Vec2d right_buffer_translation = Vec2d(DL, 0.5 * DH) + Vec2d(-0.5 * buffer_thickness, 0.0);
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
//** For getting centerline velocity *
namespace observe_centerline
{
Real x_observe_start = 0.0;
Real y_observe_start = DH / 2.0;

int num_observer_points = std::round(DL / resolution_ref / 5.0); //**Every particle is regarded as a cell monitor*
Real observe_spacing = DL / num_observer_points;

StdVec<Vecd> observation_location;
Vecd pos_observe_start = Vecd(x_observe_start, y_observe_start);
Vecd unit_direction_observe = Vecd(1.0, 0.0);
Real observer_offset_distance = 0.0;
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
void output_observer_theoretical_x()
{
    std::string filename = "../bin/output/observer_theoretical_x.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (const Vecd &position : observation_location)
    {
        outfile << position[0] << "\n";
    }
    outfile.close();
}
} // namespace observe_centerline

//** For getting cross-section velocity *
namespace observe_cross_sections
{
constexpr const char *namespace_prefix = "cross_sections";
const int number_observe_line = 5;
Real observer_offset_distance = 2.0 * resolution_ref;
Vec2d unit_direction_observe(0.0, 1.0);
// ** Determine the observing start point. *
Real observe_start_x[number_observe_line] = {
    0.0 * DL + observer_offset_distance,
    0.25 * DL,
    0.50 * DL,
    0.75 * DL,
    0.99 * DL - observer_offset_distance};

Real observe_start_y[number_observe_line] = {
    0.5 * resolution_ref,
    0.5 * resolution_ref,
    0.5 * resolution_ref,
    0.5 * resolution_ref,
    0.5 * resolution_ref};

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

//** For regression test *
StdVec<Vecd> observer_location_center_point = {Vecd(0.5 * DL, 0.5 * DH)};
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
        : u_ref_(U_inlet), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        //target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        //target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / half_channel_height / half_channel_height);
        target_velocity[0] = u_ave;
        if (1)
        {
            //** Impose fully-developed velocity from PYTHON result */
            //** Calculate the distance to wall, Y. position[1] is the distance to the centerline */
            Real Y = half_channel_height - std::abs(position[1]);
            int polynomial_order = 8;
            int num_coefficient = polynomial_order + 1;
            //** Coefficient of the polynomial, 8th-order, from py21 dp=0.024 */
            // Real coeff[] = {
            //     6.153336e-01, 3.095679e+00, -1.399783e+01,
            //     4.798221e+01, -1.100147e+02, 1.619762e+02,
            //     -1.464631e+02, 7.373006e+01, -1.577924e+01
            // };
            //** Coefficient of the polynomial, 8th-order, from py21 dp=0.1 */
            Real coeff[] = {
                6.492006e-01, 2.145673e+00, -7.442681e+00,
                2.148624e+01, -4.443593e+01, 6.171458e+01,
                -5.439313e+01, 2.726584e+01, -5.887918e+00};
            Real polynomial_value = 0.0;
            for (int i = 0; i < num_coefficient; ++i)
            {
                polynomial_value += coeff[i] * std::pow(Y, i);
            }

            if (Y > half_channel_height || Y < 0.0)
            {
                std::cout << "position[1]=" << position[1] << std::endl;
                std::cout << "Y=" << Y << std::endl;
                std::cout << "polynomial_value=" << polynomial_value << std::endl;
                std::cout << "Stop" << std::endl;
                std::cout << "=================" << std::endl;
                std::cin.get();
            }

            //** Impose inlet velocity gradually */
            target_velocity[0] = current_time < t_ref_ ? 0.5 * polynomial_value * (1.0 - cos(Pi * current_time / t_ref_)) : polynomial_value;
            //target_velocity[0] = polynomial_value;
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