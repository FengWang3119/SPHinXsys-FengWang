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
Real num_fluid_cross_section = 10.0;
Real resolution_ref = DH / num_fluid_cross_section; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                       /**< Extending width for BCs. */
Real buffer_thickness = 5.0 * resolution_ref;

Real extend_in = 0.0;
Real extend_out = 0.0;
Real extend_compensate_relaxation = 0.0;
Real DH1 = 2.0 * DH;
Real DL1 = buffer_thickness; //** Inlet velocity as uniform U_inlet */
Real DL2 = 5.0 * DH;

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

StdVec<Vecd> observer_location = {Vecd(0.5 * DL2, 0.5 * DH)}; /**< Displacement observation point. */
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
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(point_O + Vecd(0.0, -BW)); //** Keep the section neat */
        outer_wall_shape.push_back(point_A + Vecd(0.0, +BW));
        outer_wall_shape.push_back(point_B + Vecd(0.0, +BW));
        outer_wall_shape.push_back(point_C + Vecd(0.0, -BW));
        outer_wall_shape.push_back(point_O + Vecd(0.0, -BW));

        std::vector<Vecd> inner_wall_shape;

        inner_wall_shape.push_back(point_O + Vecd(-BW, 0.0));
        inner_wall_shape.push_back(point_A + Vecd(-BW, 0.0));
        inner_wall_shape.push_back(point_B + Vecd(+BW, 0.0));
        inner_wall_shape.push_back(point_C + Vecd(+BW, 0.0));
        inner_wall_shape.push_back(point_O + Vecd(-BW, 0.0));

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};

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