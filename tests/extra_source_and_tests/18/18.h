/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
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
Real DH = 2.0;  /**< Channel height. */
Real DL = 30.0; /**< Channel length. */
Real num_fluid_cross_section = 20.0;

//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------

Real resolution_ref = DH / num_fluid_cross_section; /**< Initial reference particle spacing. */

Real BW = resolution_ref * 4; /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;
Real half_channel_height = DH / 2.0;
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
Real rho0_f = 1000.0; /**< Density. */
Real Re = 200.0;

Real Outlet_pressure = 100.0;

Real mu_f = rho0_f * U_f * DH / Re;

Real Re_calculated = U_f * DH * rho0_f / mu_f;

//----------------------------------------------------------------------
//	The emitter block with offset model.
//----------------------------------------------------------------------
// Vec2d left_buffer_halfsize = Vec2d(0.5 * BW, 0.5 * DH_C + BW);
// Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize + Vecd(0.0, offset_distance - BW);
Vec2d left_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH);
Vec2d left_buffer_translation = left_buffer_halfsize + Vec2d(-DL_sponge, 0.0);

// Vec2d right_buffer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
// Vec2d right_buffer_translation = Vec2d(DL, DH + 0.25 * DH) - right_buffer_halfsize;
Vec2d right_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH);
Vec2d right_buffer_translation = Vec2d(DL - 2.5 * resolution_ref, 0.5 * DH);
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
Real x_observe_start = 0.0;
Real y_observe_start = DH / 2.0;

int num_observer_points = std::round(DL / resolution_ref / 5.0); //**Every particle is regarded as a cell monitor*
Real observe_spacing = DH / num_observer_points;

StdVec<Vecd> observation_location;
Vecd pos_observe_start = Vecd(x_observe_start, y_observe_start);
Vecd unit_direction_observe = Vecd(1.0, 0.0);
Real observer_offset_distance = 0.0;

//** For regression test *
StdVec<Vecd> observer_location_center_point = {Vecd(0.5 * DL, 0.5 * DH)};
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon computational_domain(createWaterBlockShape());
        add<ExtrudeShape<MultiPolygonShape>>(0.0, computational_domain, "ComputationalDomain");
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
        add<ExtrudeShape<MultiPolygonShape>>(BW, outer_dummy_boundary, "OuterDummyBoundary");

        MultiPolygon inner_dummy_boundary(createInnerWallShape());
        subtract<ExtrudeShape<MultiPolygonShape>>(0.0, inner_dummy_boundary, "InnerDummyBoundary");
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

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / half_channel_height / half_channel_height);
        //target_velocity[0] = u_ave;

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