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
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 1.0e-3;
Real DH = 0.015; /**< Channel height. */
Real num_fluid_cross_section = 20.0;
Real extend_in = 5.0 * DH;
Real extend_out = 3.0 * DH;
Real DL1 = 3.0 * DH + extend_in;
Real DL2 = 3.0 * DH + extend_out;
Real DL3 = 21.0 * DH;
Real incline_angle = 0.0 * (2.0 * Pi / 360.0);
Vec2d point_A(0.0, DH);
Vec2d point_B(DL1, DH);
Vec2d point_C(DL1 + DL3, DH + DL3 * tan(incline_angle));
Vec2d point_D(DL1 + DL3 + DL2, DH + DL3 * tan(incline_angle));
Vec2d point_E(DL1 + DL3 + DL2, 0.0);
//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------
Real characteristic_length = DH; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
//** For K and Epsilon, type of the turbulent inlet, 0 is freestream, 1 is from interpolation from PY21 *
int type_turbulent_inlet = 0;
Real relaxation_rate_turbulent_inlet = 0.8;
//** Tag for AMRD *
int is_AMRD = 0;
bool is_constrain_normal_velocity_in_P_region = false;
//** Weight for correcting the velocity  gradient in the sub near wall region  *
Real weight_vel_grad_sub_nearwall = 0.1;
bool is_always_lattice_arrange_fluid = false;
//** Empirical parameter for initial stability*
Real turbulent_module_activate_time = 2.5;
//** Empirical parameter for initial stability, initial time-decaying ACC*
Real duration_initial_acceleration = 1.0;
//** Intial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = {0.000180001, 3.326679e-5, 1.0e-3};

Real y_p_constant = DH / 2.0 / num_fluid_cross_section; //** For the first try *
//Real y_p_constant = 0.05;
Real resolution_ref_temp = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
Real resolution_ref = round(resolution_ref_temp * 1.0e8) / 1.0e8;
Real offset_distance = y_p_constant - resolution_ref / 2.0; //** Basically offset distance is large than or equal to 0 *

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
//** Laminar *
//Real U_inlet = 0.5;
//Real U_max = 0.75;
//Real U_f = U_inlet; //*Characteristic velo is regarded as average velo here
//Real c_f = 10.0 * U_max;                                        /**< Speed of sound. */
//Real rho0_f = 1.0;                                            /**< Density. */
//Real mu_f = 0.01;

//** Same parameter as SPH_4 *
Real U_inlet = 1.0;
Real Outlet_pressure = 0.0;
Real U_f = U_inlet;         //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */
Real Re = 18000.0;
//Real Re = 100.0;
Real mu_f = rho0_f * U_f * DH / Re;

Real Re_calculated = U_f * DH * rho0_f / mu_f;

Real DH_C = DH - 2.0 * offset_distance;

//** Intial inlet pressure to drive flow *
//Real initial_inlet_pressure = 0.5 * rho0_f * U_inlet * U_inlet;
//Real initial_inlet_pressure = 5.0;
//----------------------------------------------------------------------
//	The emitter block with offset model.
//----------------------------------------------------------------------
Real BW = resolution_ref * 4; /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;
Real half_channel_height = DH / 2.0;

Vec2d left_buffer_halfsize = Vec2d(0.5 * BW, 0.5 * DH_C + BW);
Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize + Vecd(0.0, offset_distance - BW);

Real outlet_buffer_length = BW;
Real outlet_buffer_height = 5.0 * DH;

Real outlet_disposer_rotation_angel = 0.0; //** By default, counter-clockwise is positive *
Vec2d outlet_buffer_center_translation = (point_D + point_E) / 2.0 + Vecd(-1.0, 0.0) * outlet_buffer_length / 2.0;

Real outlet_emitter_rotation_angel = Pi + outlet_disposer_rotation_angel; //** By default, counter-clockwise is positive *

//** If return to the straight channel *
// Real outlet_disposer_rotation_angel = 0.0 * Pi ; //** By default, counter-clockwise is positive *
// Real outlet_emitter_rotation_angel =  Pi ; //** By default, counter-clockwise is positive *
// Vec2d outlet_buffer_center_translation = Vec2d(DL_domain - 0.5 * outlet_buffer_length , 0.5 * DH) ;

Vec2d right_buffer_halfsize = Vec2d(0.5 * outlet_buffer_length, 0.75 * outlet_buffer_height);
Vec2d right_buffer_translation = outlet_buffer_center_translation;

//Vec2d disposer_halfsize = Vec2d(0.75 * DH, 0.5 * BW);
//Vec2d disposer_translation = Vec2d(DL_domain + 0.25 * DH, DH_domain ) - disposer_halfsize;

//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
Real DL_domain = DL1 + DL3 + DL2 + offset_distance;
Real DH_domain = point_D[1];
Vec2d left_bottom_point(-DL_sponge - offset_distance, 0.0);
Vec2d right_up_point(DL_domain, DH_domain);
BoundingBox system_domain_bounds(left_bottom_point + Vec2d(-2.0 * BW, -2.0 * BW), right_up_point + Vec2d(2.0 * BW, 2.0 * BW));
//----------------------------------------------------------------------
// Output and time average control.
//----------------------------------------------------------------------
int screen_output_interval = 100;
Real end_time = 10.0;                /**< End time. */
Real Output_Time = end_time / 200.0; /**< Time stamps for output of body states. */
Real cutoff_time = end_time * 0.0;   //** cutoff_time should be a integral and the same as the PY script */
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
// ** By kernel weight. *
int number_observe_line = 7;
Real observer_offset_distance = 1.8 * resolution_ref;
Vec2d unit_direction_observe(0.0, -1.0); //** From end to start to ensure observe accuracy on the inclined wall */
// ** Determine the observing start point. *
Real observe_start_x[7] = {
    point_B[0] - 5.87 * DH,
    point_B[0] + 2.59 * DH,
    point_B[0] + 5.98 * DH,
    point_B[0] + 13.56 * DH,
    point_B[0] + 16.93 * DH,
    point_B[0] + 20.32 * DH,
    point_B[0] + 23.71 * DH};
Real observe_start_y = 0.0 + offset_distance;
// ** Determine the length of the observing line and other information. *
Real observe_line_length[7] = {0.0};
Real observe_end_y[7] = {0.0};
int num_observer_points[7] = {0};
void getObservingLineLengthAndEndPoints()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        if (observe_start_x[i] < point_B[0])
        {
            observe_end_y[i] = observe_start_y + DH_C;
            observe_line_length[i] = DH_C;
            num_observer_points[i] = std::round(observe_line_length[i] / resolution_ref);
        }
        else if (observe_start_x[i] >= point_B[0] && observe_start_x[i] < point_C[0])
        {
            Real x_diff = observe_start_x[i] - point_B[0];
            Real y_diff = x_diff * tan(incline_angle);
            observe_end_y[i] = y_diff + point_B[1] - offset_distance;
            observe_line_length[i] = observe_end_y[i] - observe_start_y;
            num_observer_points[i] = std::round(observe_line_length[i] / resolution_ref);
        }
        else if (observe_start_x[i] >= point_C[0])
        {
            observe_end_y[i] = point_D[1] - offset_distance;
            observe_line_length[i] = point_D[1] - 2.0 * offset_distance;
            num_observer_points[i] = std::round(observe_line_length[i] / resolution_ref);
        }
    }
}

StdVec<Vecd> observation_locations;
void getPositionsOfMultipleObserveLines()
{
    getObservingLineLengthAndEndPoints();
    for (int k = 0; k < number_observe_line; ++k)
    {
        Vecd pos_observe_start(observe_start_x[k], observe_end_y[k]);
        int num_observer_point = num_observer_points[k];
        Real observe_spacing = observe_line_length[k] / num_observer_point;
        for (int i = 0; i < num_observer_point; ++i)
        {
            Real offset = 0.0;
            offset = (i == 0 ? -observer_offset_distance : (i == num_observer_point - 1 ? observer_offset_distance : 0.0));
            Vecd pos_observer_i = pos_observe_start + (i * observe_spacing + offset) * unit_direction_observe;
            observation_locations.push_back(pos_observer_i);
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
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------

std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;

    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, DH));
    water_block_shape.push_back(point_A);
    water_block_shape.push_back(point_B);
    water_block_shape.push_back(point_C);
    water_block_shape.push_back(point_D);
    water_block_shape.push_back(Vecd(point_D[0] + offset_distance, point_D[1]));
    water_block_shape.push_back(Vecd(point_D[0] + offset_distance, point_E[1]));
    water_block_shape.push_back(point_E);
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

    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - BW, DH));
    water_block_shape.push_back(point_A);
    water_block_shape.push_back(point_B);
    water_block_shape.push_back(point_C);
    water_block_shape.push_back(point_D);
    water_block_shape.push_back(Vecd(point_D[0] + offset_distance + BW, point_D[1]));
    water_block_shape.push_back(Vecd(point_D[0] + offset_distance + BW, point_E[1]));
    water_block_shape.push_back(point_E);
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - BW, 0.0));

    return water_block_shape;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> water_block_shape;

    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - 2.0 * BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - 2.0 * BW, DH));
    water_block_shape.push_back(point_A);
    water_block_shape.push_back(point_B);
    water_block_shape.push_back(point_C);
    water_block_shape.push_back(point_D);
    water_block_shape.push_back(Vecd(point_D[0] + offset_distance + 2.0 * BW, point_D[1]));
    water_block_shape.push_back(Vecd(point_D[0] + offset_distance + 2.0 * BW, point_E[1]));
    water_block_shape.push_back(point_E);
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - 2.0 * BW, 0.0));

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
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_inlet), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        //target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        //target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / half_channel_height / half_channel_height);
        target_velocity[0] = u_ave;
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
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
    Real t_ref_, u_ref_, du_ave_dt_;

  public:
    explicit TimeDependentAcceleration(Vecd gravity_vector)
        : Gravity(gravity_vector), t_ref_(duration_initial_acceleration), u_ref_(U_inlet), du_ave_dt_(0) {}

    virtual Vecd InducedAcceleration(const Vecd &position) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
    }
};

struct RightOutflowPressure
{
    template <class BoundaryConditionType>
    RightOutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
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

    Real operator()(Real &p_)
    {
        return p_;
    }
};