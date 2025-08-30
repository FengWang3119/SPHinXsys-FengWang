#include "mph_method.cpp"
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup,
//	Dimension: m s kg
//	Reproduce [2021 Kondo CPM]
//----------------------------------------------------------------------
Real DL = 3.22;                    /**< Water tank length. */
Real LL = 1.2;                      /**< Water column length. */
Real LH = 0.6;                      /**< Water column height. */
Real DH = 4.0 * LH;                    /**< Water tank height. */
Real particle_spacing_ref = 0.02;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                       /**< Reference density of fluid. */
Real gravity_g = 9.8;                    /**< Gravity. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	MPH Parameters.
//----------------------------------------------------------------------
Real bulk_modulus_ref = 1.0e+6;
Real bulk_viscosity_ref = 2.0e+3;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};