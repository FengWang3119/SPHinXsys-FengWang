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
Real DH = 2.0; /**< Channel height. */
Real num_fluid_cross_section = 100.0;
Real DL = 3.0;

Real y_p_constant = DH / 2.0 / num_fluid_cross_section; //%% For the first try *
//Real y_p_constant = 0.025;

Real resolution_ref_temp = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
Real resolution_ref = round(resolution_ref_temp * 1.0e8) / 1.0e8;
Real offset_distance = y_p_constant - resolution_ref / 2.0; //%% Basically offset distance is large than or equal to 0 *
Real BW = resolution_ref * 4;

Vec2d point_O(0.0, 0.0);

Vec2d point_1(2.42139, 0.95617);
Vec2d point_2(1.82129, 0.95617);
Vec2d point_3(1.82129, 0.844755);
Vec2d point_4(2.41774, 0.612);
Vec2d point_5(2.498, 0.692264);
Vec2d point_6(2.63, 0.0);
Vec2d point_7(2.63, 0.4);
Vec2d point_8(2.705, 0.475);
Vec2d point_9(2.78, 0.4);
Vec2d point_10(2.78, 0.0);

Vec2d point_11 = point_O + Vec2d(-BW, DH);
Vec2d point_12(point_1[xAxis], point_11[yAxis]);
Vec2d point_13(point_5[xAxis], point_11[yAxis]);
Vec2d point_14 = point_11 + Vec2d(2.0 * BW + DL, 0.0);
Vec2d point_15 = point_14 + Vec2d(0.0, BW);
Vec2d point_16 = point_15 + Vec2d(-2.0 * BW - DL, 0.0);

Vec2d point_17 = point_O + Vec2d(-BW, -BW);
Vec2d point_18 = point_17 + Vec2d(2.0 * BW + DL, 0.0);
Vec2d point_19 = point_18 + Vec2d(0.0, BW);
Vec2d point_20 = point_O + Vec2d(-BW, 0.0);
//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------
//Real characteristic_length = DH; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
//** For K and Epsilon, type of the turbulent inlet, 0 is freestream, 1 is from interpolation from PY21 *
//int type_turbulent_inlet = 0;
//Real relaxation_rate_turbulent_inlet = 0.8;
//** Tag for AMRD *
int is_AMRD = 1;
bool is_constrain_normal_velocity_in_P_region = false;
//** Weight for correcting the velocity  gradient in the sub near wall region  *
Real weight_vel_grad_sub_nearwall = 0.1;
bool is_always_lattice_arrange_fluid = false;
//** Tag for Source Term Linearisation *
bool is_source_term_linearisation = true;
//** Empirical parameter for initial stability*
Real turbulent_module_activate_time = 0.0;
//** Initial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = {0.000180001, 3.326679e-5, 1.0e-3};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
//Real U_inlet = 1.0;
//Real Outlet_pressure = 0.0;
//Real U_f = U_inlet;         //*Characteristic velocity
//Real U_max = 3.0 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real DL_Sponge = 5.0 * resolution_ref;
Vec2d buffer_halfsize = 0.5 * Vec2d(DL_Sponge, DH);
Vec2d buffer_translation = Vec2d(DL_Sponge / 2.0, DH / 2.0) + Vec2d(BW, 0.0);

/*
//Vecd external_acc = Vecd(1.885, 0.0);
//Real external_acc_gradually_impose_t = 2.0;
*/

Real axis_vel_ref_ = 0.33;
Real external_acc_initial = 0.01;

/*
Vecd target_velocity = Vecd(38.0, 0.0);
Vecd external_acceleration_initial = Vecd(0.1, 0.0);
*/
Real U_max = 3.0; //** An estimated value, Periodic BC *
Real U_f = U_max; // * Characteristic velocity, Periodic BC

Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */
//Real Re = 40000.0;
//Real Re = 100.0;
//Real mu_f = rho0_f * U_f * DH / Re;

Real mu_f = 1.0e-6; //** Periodic BC */

Real Re_calculated = U_f * DH * rho0_f / mu_f;

Real DH_C = DH - 2.0 * offset_distance;

//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
// Real DL_domain = DL;
// Real DH_domain = DH;
Vec2d left_bottom_point = point_17;
Vec2d right_up_point = point_15;
BoundingBox system_domain_bounds(left_bottom_point + Vec2d(-2.0 * BW, -2.0 * BW), right_up_point + Vec2d(2.0 * BW, 2.0 * BW));
//----------------------------------------------------------------------
// Output and time average control.
//----------------------------------------------------------------------
int screen_output_interval = 100;
Real end_time = 30.0;                /**< End time. */
Real Output_Time = end_time / 200.0; /**< Time stamps for output of body states. */
Real cutoff_time = 25.0;             //%% cutoff_time should be a integral and the same as the PY script */
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
// ** By kernel weight. *
const int number_observe_line = 3;
Real observer_offset_distance = 2.0 * resolution_ref;
Vec2d unit_direction_observe(0.0, 1.0);
// ** Determine the observing start point. *
Real observe_start_x[number_observe_line] = {
    0.5,
    1.0,
    1.5};
Real observe_start_y[number_observe_line] = {
    point_O[yAxis] + 0.5 * resolution_ref + offset_distance,
    point_O[yAxis] + 0.5 * resolution_ref + offset_distance,
    point_O[yAxis] + 0.5 * resolution_ref + offset_distance};
// ** Determine the length of the observing line and other information. *
Real observe_line_length[number_observe_line] = {0.0};
int num_observer_points[number_observe_line] = {0};
void getObservingLineLengthAndEndPoints()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        observe_line_length[i] = DH_C;
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
void output_observe_theoretical_y()
{
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

//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;

    water_block_shape.push_back(point_O);
    water_block_shape.push_back(point_O + Vec2d(0.0, DH));
    water_block_shape.push_back(point_O + Vec2d(DL, DH));
    water_block_shape.push_back(point_O + Vec2d(DL, 0.0));
    water_block_shape.push_back(point_O);

    return water_block_shape;
}

//% For fish pass geometry
std::vector<Vecd> createUpperWallShape()
{
    std::vector<Vecd> shape;

    shape.push_back(point_11);
    shape.push_back(point_16);
    shape.push_back(point_15);
    shape.push_back(point_14);
    shape.push_back(point_13);
    shape.push_back(point_5);
    shape.push_back(point_4);
    shape.push_back(point_3);
    shape.push_back(point_2);
    shape.push_back(point_1);
    shape.push_back(point_12);
    shape.push_back(point_11);

    return shape;
}
std::vector<Vecd> createBottomWallShape()
{
    std::vector<Vecd> shape;

    shape.push_back(point_17);
    shape.push_back(point_20);
    shape.push_back(point_6);
    shape.push_back(point_7);
    shape.push_back(point_8);
    shape.push_back(point_9);
    shape.push_back(point_10);
    shape.push_back(point_19);
    shape.push_back(point_18);
    shape.push_back(point_17);

    return shape;
}
std::vector<Vecd> createHookWallShape()
{
    std::vector<Vecd> shape;

    shape.push_back(point_12);
    shape.push_back(point_13);
    shape.push_back(point_5);
    shape.push_back(point_4);
    shape.push_back(point_3);
    shape.push_back(point_2);
    shape.push_back(point_1);
    shape.push_back(point_12);

    return shape;
}
std::vector<Vecd> createPencilWallShape()
{
    std::vector<Vecd> shape;

    shape.push_back(point_6);
    shape.push_back(point_7);
    shape.push_back(point_8);
    shape.push_back(point_9);
    shape.push_back(point_10);
    shape.push_back(point_6);

    return shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon computational_domain(createWaterBlockShape());
        add<MultiPolygonShape>(computational_domain, "ComputationalDomain");

        MultiPolygon sub_upper_dummy_boundary(createHookWallShape());
        subtract<MultiPolygonShape>(sub_upper_dummy_boundary, "SubUpperDummyBoundary");
        MultiPolygon sub_bottom_dummy_boundary(createPencilWallShape());
        subtract<MultiPolygonShape>(sub_bottom_dummy_boundary, "SubBottomDummyBoundary");
    }
};

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon upper_dummy_boundary(createUpperWallShape());
        add<MultiPolygonShape>(upper_dummy_boundary, "UpperDummyBoundary");
        MultiPolygon bottom_dummy_boundary(createBottomWallShape());
        add<MultiPolygonShape>(bottom_dummy_boundary, "BottomDummyBoundary");
    }
};

/*
//% For straight channel geometry
std::vector<Vecd> createUpperWallShape()
{
    std::vector<Vecd> shape;

    shape.push_back(point_11);
    shape.push_back(point_16);
    shape.push_back(point_15);
    shape.push_back(point_14);
    shape.push_back(point_11);

    return shape;
}
std::vector<Vecd> createBottomWallShape()
{
    std::vector<Vecd> shape;

    shape.push_back(point_17);
    shape.push_back(point_20);
    shape.push_back(point_19);
    shape.push_back(point_18);
    shape.push_back(point_17);

    return shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon computational_domain(createWaterBlockShape());
        add<MultiPolygonShape>(computational_domain, "ComputationalDomain");
    }
};

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon upper_dummy_boundary(createUpperWallShape());
        add<MultiPolygonShape>(upper_dummy_boundary, "UpperDummyBoundary");
        MultiPolygon bottom_dummy_boundary(createBottomWallShape());
        add<MultiPolygonShape>(bottom_dummy_boundary, "BottomDummyBoundary");
    }
};
*/