#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "k-epsilon_turbulent_model.cpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real plate_length = 100.0;
Real num_fluid_cross_plate_length = 800.0;
Real resolution_ref = plate_length / num_fluid_cross_plate_length; /**< Initial reference particle spacing. */

Real y_p_constant = resolution_ref / 2.0; //** For the first try *
// Real y_p_constant = 0.05;
Real offset_distance = y_p_constant - resolution_ref / 2.0; //** Basically offset distance is large than or equal to 0 *
//----------------------------------------------------------------------
//	Unique control parameters for turbulence.
//----------------------------------------------------------------------
int inflow_velocity_type = 0;           // ** 0 uniform, 1 parabolic
int is_inflow_velocity_from_python = 0; // ** Overwrite the inflow_velocity_type by input from python

Real characteristic_length = plate_length; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
//** For K and Epsilon, type of the turbulent inlet, 0 is freestream, 1 is from interpolation from PY21, 2 is determined from inital values *
int type_turbulent_inlet = 0;
Real relaxation_rate_turbulent_inlet = 0.8;
//** Tag for AMRD *
int is_AMRD = 0;
bool is_constrain_normal_velocity_in_P_region = false;
//** Weight for correcting the velocity  gradient in the sub near wall region  *
Real weight_vel_grad_sub_nearwall = 0.1;
//** Tag for Source Term Linearisation *
bool is_source_term_linearisation = false;

//----------------------------------------------------------------------
//	Geometry settings.
//----------------------------------------------------------------------
Real plate_thickness = resolution_ref * 4;

Real BW = resolution_ref * 4; /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;

Real DH_half_flow_region = 0.03125 * plate_length;

Real DL_flow_region_temp = 1.5 * plate_length;
Real DL_flow_region = round(DL_flow_region_temp * 1.0e8) / 1.0e8; //* This size will affect define levelset, larger than 100 will cause failure

Real DL_upstream_distance = 0.25 * plate_length;

Vec2d point_O(0.0, 0.0);
Vec2d point_A = point_O + Vec2d(0.0, 2.0 * DH_half_flow_region + plate_thickness);
Vec2d point_B = point_A + Vec2d(DL_flow_region, 0.0);
Vec2d point_C = point_B - Vec2d(0.0, 2.0 * DH_half_flow_region + plate_thickness);

Vec2d point_D(DL_upstream_distance, DH_half_flow_region);
Vec2d point_E = point_D + Vec2d(0.0, plate_thickness);
Vec2d point_F = point_E + Vec2d(plate_length, 0.0);
Vec2d point_G = point_F - Vec2d(0.0, plate_thickness);

Real DH_total = point_A[yAxis] - point_O[yAxis];
Real DL_total = point_C[xAxis] - point_O[xAxis];
Real half_DH_total = DH_total / 2.0;
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
Vec2d point_left_down = point_O - 2.0 * Vec2d(BW, BW);
Vec2d point_right_up = point_B + 2.0 * Vec2d(BW, BW);
BoundingBox system_domain_bounds(point_left_down, point_right_up);

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real U_inlet = 1.0;
Real U_f = U_inlet;         //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */
Real Re = 4.2e6;

Real Outlet_pressure = 0.0;

Real mu_f = rho0_f * U_f * plate_length / Re;

Real Re_calculated = U_f * plate_length * rho0_f / mu_f;

//** Initial values for K, Epsilon and Mu_t *
Real initial_k = 0.013 * U_inlet * U_inlet;
Real initial_mut = 367.0 * mu_f;
Real initial_epsilon = 0.09 * rho0_f * initial_k * initial_k / initial_mut;
StdVec<Real> initial_turbu_values = {initial_k, initial_epsilon, initial_mut};
//----------------------------------------------------------------------
//	The emitter block with offset model.
//----------------------------------------------------------------------
Vec2d left_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH_total);
Vec2d left_buffer_translation = left_buffer_halfsize + Vec2d(-DL_sponge, 0.0);

Vec2d right_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH_total);
Vec2d right_buffer_translation = Vec2d(DL_total - 2.5 * resolution_ref, 0.5 * DH_total);
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
Real x_observe_start = 0.99 * DL_total;
int num_observer_points = std::round(DH_total / resolution_ref); //**Every particle is regarded as a cell monitor*
Real observe_spacing = DH_total / num_observer_points;

// By kernel weight.
StdVec<Vecd> observation_location;
Vecd pos_observe_start = Vecd(x_observe_start, resolution_ref / 2.0 + offset_distance);
Vecd unit_direction_observe = Vecd(0.0, 1.0);
Real observer_offset_distance = 2.0 * resolution_ref;

//** For regression test *
StdVec<Vecd> observer_location_center_point = {Vecd(0.5 * DL_total, 0.5 * DH_total)};
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape() //* Outflow case, no need to offset water block
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH_total));
    water_block_shape.push_back(Vecd(DL_total, DH_total));
    water_block_shape.push_back(Vecd(DL_total, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    return water_block_shape;
}

std::vector<Vecd> createFlatPlate()
{
    std::vector<Vecd> block_shape;
    block_shape.push_back(point_D);
    block_shape.push_back(point_E);
    block_shape.push_back(point_F);
    block_shape.push_back(point_G);
    block_shape.push_back(point_D);

    return block_shape;
}

// std::vector<Vecd> createOffsetFlatPlate()
// {
//     std::vector<Vecd> block_shape;
//     block_shape.push_back(point_D + Vec2d(-offset_distance, -offset_distance));
//     block_shape.push_back(point_E + Vec2d(-offset_distance, offset_distance));
//     block_shape.push_back(point_F + Vec2d(offset_distance, offset_distance));
//     block_shape.push_back(point_G + Vec2d(offset_distance, -offset_distance));
//     block_shape.push_back(point_D + Vec2d(-offset_distance, -offset_distance));

//     return block_shape;
// }

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_dummy_boundary(createFlatPlate());
        add<ExtrudeShape<MultiPolygonShape>>(offset_distance, outer_dummy_boundary, "OuterDummyBoundary");
    }
};

class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon computational_domain(createWaterBlockShape());
        add<MultiPolygonShape>(computational_domain, "ComputationalDomain");
        MultiPolygon solid_space(createFlatPlate());
        subtract<ExtrudeShape<MultiPolygonShape>>(offset_distance, solid_space, "SubSolid");
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
        switch (inflow_velocity_type)
        {
        case 0:
            target_velocity[0] = u_ave;
            break;
        case 1:
            target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / half_DH_total / half_DH_total);
            break;
        }

        if (is_inflow_velocity_from_python)
        {
            //** Impose fully-developed velocity from PYTHON result */
            //** Calculate the distance to wall, Y. position[1] is the distance to the centerline */
            Real Y = half_DH_total - std::abs(position[1]);
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

            if (Y > half_DH_total || Y < 0.0)
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
            // target_velocity[0] = polynomial_value;
        }

        if (position[1] > half_DH_total)
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