#include "bidirectional_buffer.h"
#include "udf_common_turbulence_model.cpp"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "udf_k-omega_turbulent_model.cpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "udf_P_refinement.h"
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 2.0;  /**< Channel height. */
Real DL = 30.0; /**< Channel length. */
Real num_fluid_cross_section = 320.0;

Real time_gradually_increase_vel = 2.0;
//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------
Real characteristic_length = DH; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
//** For K and Epsilon/Omega, type of the turbulent inlet, 0 is freestream, 1 is from interpolation from PY21, 2 is from PY2-11 *
int type_turbulent_inlet_omega = 2;
//** For K and Epsilon/Omega, type of the turbulent inlet, 0 is freestream, 1 is from interpolation from PY21, 2 is from PY2-11*
int type_turbulent_inlet_k = 2;
// ** 0 is freestream, 1 is from interpolation from PY21, 2 is from OF6-28, 3 is from PY2-11 *
int type_velocity_inlet = 3;

Real relaxation_rate_turbulent_inlet = 0.8;
//** Tag for wall treatment *
int is_blended = 0;
//** Tag for AMRD *
int is_AMRD = 0;
bool is_constrain_normal_velocity_in_P_region = false;
//** Weight for correcting the velocity  gradient in the sub near wall region  *
//Real weight_vel_grad_sub_nearwall = 0.1;
//** Tag for Source Term Linearisation *
bool is_source_term_linearisation = false;
//** Empirical parameter for initial stability*
Real turbulent_module_activate_time = 2.0;
//** Initial values for K, Omega and Mu_t *
StdVec<Real> initial_turbu_values = {0.01, 2.056, 0.02};

//Real y_p_constant = 0.05;
//Real y_p_constant = DH / 2.0 / num_fluid_cross_section; //** For the first try or Not use BOT *
//Real resolution_ref = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
//Real offset_distance = y_p_constant - resolution_ref / 2.0;                        //** Basically offset distance is large than or equal to 0 *

// ** If not use BOT *
Real y_p_constant = DH / 2.0 / num_fluid_cross_section;            
Real resolution_ref = DH / num_fluid_cross_section;
Real offset_distance = 0.0;                        

Real BW = resolution_ref * 4; /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;
Real half_channel_height = DH / 2.0;
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + 2.0 * BW));

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real U_inlet = 1.0;
Real U_f = U_inlet;         //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */

//Real Re = 5714.0;
Real Re = 8.0e7;

Real Outlet_pressure = 0.0;

Real mu_f = rho0_f * U_f * DH / Re;

Real Re_calculated = U_f * DH * rho0_f / mu_f;

Real DH_C = DH - 2.0 * offset_distance;
//----------------------------------------------------------------------
//	The emitter block with offset model.
//----------------------------------------------------------------------
Vec2d left_buffer_halfsize = Vec2d(0.5 * BW, 0.5 * DH_C + BW);
Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize + Vecd(0.0, offset_distance - BW);

Vec2d right_buffer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d right_buffer_translation = Vec2d(DL, DH + 0.25 * DH) - right_buffer_halfsize;
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
Real x_observe_start = 0.99 * DL;
int num_observer_points = std::round(DH_C / resolution_ref); //**Every particle is regarded as a cell monitor*
Real observe_spacing = DH_C / num_observer_points;

// By kernel weight.
StdVec<Vecd> observation_location;
StdVec<Vecd> observation_theoretical_locations;
Vecd pos_observe_start = Vecd(x_observe_start, resolution_ref / 2.0 + offset_distance);
Vecd unit_direction_observe = Vecd(0.0, 1.0);
Real observer_offset_distance = 2.0 * resolution_ref;

void get_observation_locations()
{
    for (int i = 0; i < num_observer_points; ++i)
    {
        Vecd pos_observer_i = pos_observe_start + i * observe_spacing * unit_direction_observe;
        Vecd pos_observer_i_no_offset = pos_observe_start + i * observe_spacing * unit_direction_observe;
        if (i == 0)
        {
            pos_observer_i -= observer_offset_distance * unit_direction_observe;
        }
        if (i == num_observer_points - 1)
        {
            pos_observer_i += observer_offset_distance * unit_direction_observe;
        }
        observation_location.push_back(pos_observer_i);
        observation_theoretical_locations.push_back(pos_observer_i_no_offset);
    }
}
void output_observer_theoretical_y()
{
    std::cout << "offset_distance=" << offset_distance << std::endl;
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
//** For regression test *
StdVec<Vecd> observer_location_center_point = {Vecd(0.5 * DL, 0.5 * DH)};

//** For getting near wall friction velocity or wall shear stress *
namespace observe_nearwall
{
constexpr const char *namespace_prefix = "nearwall";
const int number_observe_line = 1;
Real sparse_ratio = 2.0;
Real observer_offset_distance = 0.0 * resolution_ref; //** Offset the first and last observing point *
Real observer_offset_distance_whole_line = 2.0 * resolution_ref;
Vec2d unit_direction_observe(1.0, 0.0);
// ** Determine the observing start point of the each line. *
Real observe_start_x[number_observe_line] = {0.0};
Real observe_start_y[number_observe_line] = {0.5 * resolution_ref - observer_offset_distance_whole_line};
// ** Determine the length of the observing line and other information. *
Real observe_line_length[number_observe_line] = {0.0};
int num_observer_points[number_observe_line] = {0};
void getObservingLineLengthAndEndPoints()
{
    for (int i = 0; i < number_observe_line; ++i)
    {
        observe_line_length[i] = DL;
        num_observer_points[i] = std::round(observe_line_length[i] / resolution_ref / sparse_ratio);
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
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_positions.dat";
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
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_theoretical_x.dat";
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
    std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_num_points_on_lines.dat";
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
} // namespace observe_nearwall

namespace observe_node_cross_section
{
    Real x_observe_node = x_observe_start;
    // By kernel weight.
    StdVec<Vecd> observation_location;
    StdVec<Vecd> observation_theoretical_locations;
    void get_observation_locations()
    {
        observation_location.push_back(Vec2d(x_observe_node, resolution_ref / 2.0 - observer_offset_distance));
        observation_location.push_back(Vec2d(x_observe_node, (DH - resolution_ref / 2.0) + observer_offset_distance));
        observation_theoretical_locations.push_back(Vec2d(x_observe_node, resolution_ref / 2.0));
        observation_theoretical_locations.push_back(Vec2d(x_observe_node, DH - resolution_ref / 2.0));
    }
    constexpr const char* namespace_prefix = "Node_";
    void output_observe_positions()
    {
        std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_positions.dat";
        std::ofstream outfile(filename);
        if (!outfile.is_open())
        {
            std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
            return;
        }
        for (const Vecd& position : observation_location)
        {
            outfile << position[0] << " " << position[1] << "\n";
        }
        outfile.close();
    }
    void output_observer_theoretical_y()
    {
        std::string filename = "../bin/output/" + std::string(namespace_prefix) + "observer_theoretical_y.dat";
        std::ofstream outfile(filename);
        if (!outfile.is_open())
        {
            std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
            return;
        }
        for (const Vecd& position : observation_theoretical_locations)
        {
            outfile << position[1] << "\n";
        }
        outfile.close();
    }
} // namespace observe_node_cross_section
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
        : u_ref_(U_inlet), t_ref_(time_gradually_increase_vel),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        //target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        //target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / half_channel_height / half_channel_height);
        target_velocity[0] = u_ave;
        if (type_velocity_inlet == 1)
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

        if (type_velocity_inlet == 2)
        {
            //** Impose fully-developed velocity from OF6-28 result */
            //** Calculate the distance to wall, Y. position[1] is the distance to the centerline */
            Real Y = half_channel_height - std::abs(position[1]);

            //** 2 segments *
            const Real y1 = 0.2;

            Real polynomial_value = 0.0;
            if (Y > 0.0 && Y <= y1)
            {
                static const std::vector<Real> coeff1 = {
                -1.658125e-04, 1.081005e+01, 3.522822e-01,
                -1.026430e+03, 6.581696e+04, -2.817724e+06,
                6.091954e+07, -7.664212e+08, 6.050634e+09,
                -3.059173e+10, 9.656963e+10, -1.738073e+11,
                1.364291e+11,
                };
                polynomial_value = polyEval(coeff1, Y);
                //std::cout << "polynomial_value=" << polynomial_value << " Y=" << Y << "========1=========\n";
            }
            else
            {
                static const std::vector<Real> coeff2 = {
                2.900603e-01, 6.951898e+00, -3.804471e+01,
                1.419110e+02, -3.399496e+02, 4.913445e+02,
                -3.197135e+02, -1.928404e+02, 5.883885e+02,
                -5.044169e+02, 1.846665e+02, -9.893014e+00,
                -7.546598e+00,
                };
                polynomial_value = polyEval(coeff2, Y);
                //std::cout << "polynomial_value=" << polynomial_value << " Y=" << Y << "========3=========\n";
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
        if (type_velocity_inlet == 3)
        {
            //** Impose fully-developed velocity from PY2-11 result */
            //** Calculate the distance to wall, Y. position[1] is the distance to the centerline */
            Real Y = half_channel_height - std::abs(position[1]);

            Real profile_value = inletProfileTable(Y);

            if (Y > half_channel_height || Y < 0.0)
            {
                std::cout << "position[1]=" << position[1] << std::endl;
                std::cout << "Y=" << Y << std::endl;
                std::cout << "profile_value=" << profile_value << std::endl;
                std::cout << "Stop" << std::endl;
                std::cout << "=================" << std::endl;
                std::cin.get();
            }

            //** Impose inlet velocity gradually */
            target_velocity[0] = current_time < t_ref_ ? 0.5 * profile_value * (1.0 - cos(Pi * current_time / t_ref_)) : profile_value;
            //target_velocity[0] = profile_value;
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
    Real polyEval(const std::vector<Real>& a, Real x)
    {
        const int n = static_cast<int>(a.size());
        if (n == 0)
        {
            std::cout << "size of coefficient = 0 " << '\n'
                << "=================\n";
            std::cin.get();
        }
        Real result = a[n - 1];
        for (int i = n - 2; i >= 0; --i)
        {
            result = result * x + a[i];
        }
        return result;
    }
    Real inletProfileChebyshev(Real Y)
    {
        const Real Y_min = 3.81517410e-03;
        const Real Y_max = 9.99451160e-01;

        Y = SMAX(Y_min, SMIN(Y, Y_max));

        //** X_min = std::log(Y_min), X_max = std::log(Y_max) *
        const Real xmin = -5.5687689796382367e+00;
        const Real xmax = -5.4899066780366194e-04;
        const Real inv_range = 1.0 / (xmax - xmin);

        const Real X = std::log(Y);
        const Real x = (2.0 * X - xmin - xmax) * inv_range;

        const Real c[] = {
            8.7874536325444819e-01,
            1.8682245615059126e-01,
            -1.8117146154221684e-03,
            -3.1827710521799830e-03,
            -2.5781840831446594e-03,
            -9.6544247521810548e-04,
            -1.2917847312668713e-03,
        };

        const int n = sizeof(c) / sizeof(c[0]);

        Real b1 = 0.0;
        Real b2 = 0.0;

        for (int k = n - 1; k >= 1; --k)
        {
            Real b0 = 2.0 * x * b1 - b2 + c[k];
            b2 = b1;
            b1 = b0;
        }

        return x * b1 - b2 + c[0];
    }
    Real inletProfileTable(Real Y)
    {
        const Real Y_min = 3.81517410e-03;
        const Real Y_max = 9.99451160e-01;

        constexpr int N = 2000;
        const Real inv_dY = Real(N - 1) / (Y_max - Y_min);

        static bool initialized = false;
        static std::array<Real, N> table;

        if (!initialized)
        {
            for (int i = 0; i < N; ++i)
            {
                Real Yi = Y_min + Real(i) / inv_dY;
                table[i] = inletProfileChebyshev(Yi);
            }

            initialized = true;
        }

        Y = SMAX(Y_min, SMIN(Y, Y_max));

        Real s = (Y - Y_min) * inv_dY;
        int i = static_cast<int>(s);

        if (i >= N - 1)
            i = N - 2;
        if (i < 0)
            i = 0;

        Real w = s - Real(i);

        return (1.0 - w) * table[i] + w * table[i + 1];
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

namespace SPH 
{
//=================================================================================================//
    class UpdateVolume : public LocalDynamics
    {
    public:
        explicit UpdateVolume(SPHBody& sph_body)
            : LocalDynamics(sph_body),
            rho_(particles_->getVariableDataByName<Real>("Density")),
            mass_(particles_->getVariableDataByName<Real>("Mass")),
            Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
            indicator_(particles_->getVariableDataByName<int>("Indicator")) {}
        virtual ~UpdateVolume() {};

        void update(size_t index_i, Real dt = 0.0) 
        {
            Vol_[index_i] = mass_[index_i] / rho_[index_i];
        };

    protected:
        Real* rho_, * mass_, * Vol_;
        int* indicator_;
    };
//=================================================================================================//
}
