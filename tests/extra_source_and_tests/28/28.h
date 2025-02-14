/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "k-epsilon_turbulent_model.cpp"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366; /**< Water tank length. */
Real DH = 5.366; /**< Water tank height. */
Real LL = 2.0;   /**< Water column length. */
Real LH = 1.0;   /**< Water column height. */
//Real particle_spacing_ref = 0.02;   /**< Initial reference particle spacing. */
//Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Resolution for turbulence.
//----------------------------------------------------------------------
Real num_fluid_cross_section = 50.0;
Real particle_spacing_estimated = LH / num_fluid_cross_section;

//Real y_p_constant = 0.05;
Real y_p_constant = particle_spacing_estimated / 2.0; //** For first try */

Real particle_spacing_ref = (LH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
Real offset_distance = y_p_constant - particle_spacing_ref / 2.0;                        //** Basically offset distance is large than or equal to 0 *
//----------------------------------------------------------------------
//	Other parameters.
//----------------------------------------------------------------------
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------
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
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Reference density of fluid. */
Real gravity_g = 1.0;                    /**< Gravity. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
// Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
// Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
// Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
// Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
// Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
// Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-BW, -BW));
    outer_wall_shape.push_back(Vecd(-BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, -BW));
    outer_wall_shape.push_back(Vecd(-BW, -BW));
    return outer_wall_shape;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(0.0, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, DH));
    inner_wall_shape.push_back(Vecd(DL, DH));
    inner_wall_shape.push_back(Vecd(DL, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, 0.0));
    return inner_wall_shape;
}
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
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, LH));
    water_block_shape.push_back(Vecd(LL, LH));
    water_block_shape.push_back(Vecd(LL, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));
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
//----------------------------------------------------------------------
//	Observation
//----------------------------------------------------------------------
MultiPolygon createProbeShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(LL, 0.0));
    pnts.push_back(Vecd(LL, BW));
    pnts.push_back(Vecd(DL, BW));
    pnts.push_back(Vecd(DL, 0.0));
    pnts.push_back(Vecd(LL, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}