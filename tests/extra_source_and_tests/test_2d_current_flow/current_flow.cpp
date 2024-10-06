/**
 * @file 	current_flow.cpp
 * @brief 	
 * @details 
 * @author 	
 */
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
Real DL = 10.0;              /**< Reference length. */
Real DH = 0.5;               /**< Reference and the height of main channel. */
Real resolution_ref = 0.015; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;
Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real gravity_g = 9.81;
Real wave_angular_freq = 2.0 * Pi;
Real wave_amplitude = 0.05;
Real wave_phase = 0.0;
Real Outlet_pressure = 0.0;
Real rho0_f = 1000.0;
Real emitter_velocity = 0.4;
Real emitter_time = 1; /**< Free stream time. */
Real U_f = 2.0 * sqrt(0.5 * gravity_g);
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f;
Real mu_f = 1.0e-3;

Vec2d emitter_buffer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;

Vec2d disposer_buffer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_buffer_translation = Vec2d(DL, 1.5 * DH) - disposer_buffer_halfsize;
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(-DL_sponge, 0.0)};
/** the bottom wall polygon. */
std::vector<Vecd> bottom_wall_shape{
    Vecd(-DL_sponge - 2.0 * BW, -BW), Vecd(-DL_sponge - 2.0 * BW, 0.0), Vecd(DL + 2.0 * BW, 0.0), Vecd(DL + 2.0 * BW, -BW), Vecd(-DL_sponge - 2.0 * BW, -BW)};

MultiPolygon createDampingBufferShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(DL - 3.0 - DL_sponge, 0.0));
    pnts.push_back(Vecd(DL - 3.0 - DL_sponge, DH));
    pnts.push_back(Vecd(DL - DL_sponge, DH));
    pnts.push_back(Vecd(DL - DL_sponge, 0.0));
    pnts.push_back(Vecd(DL - 3.0 - DL_sponge, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(bottom_wall_shape, ShapeBooleanOps::add);
    }
};

struct NonPrescribedPressure
{
    template <class BoundaryConditionType>
    NonPrescribedPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

struct OutletPressure
{
    template <class BoundaryConditionType>
    OutletPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real curent_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

class WaveCalculator
{
  private:
    Real gravity_;
    Real water_depth_;
    Real wave_angular_freq_;
    Real wave_amplitude_;
    Real wave_phase_;
    Real free_surface_height_;
    Real wave_number_;

    void computeWaveNumber()
    {
        int max_iterations = 20;
        Real tolerance = 1.0e-6;
        Real wave_number = 1.0;
        for (int i = 1; i <= max_iterations; ++i)
        {
            Real term1 = tanh(wave_number * water_depth_);
            Real term2 = (wave_angular_freq_ * wave_angular_freq_) / gravity_;
            Real term3 = wave_number * term1 - term2;
            Real term4 = term1 + wave_number * water_depth_ * (1.0 - term1 * term1);

            Real wave_number_old = wave_number;

            wave_number = wave_number_old - term3 / term4;

            Real error = abs(wave_number - wave_number_old) / abs(wave_number);

            if (error <= tolerance)
                break;
        }
        wave_number_ = wave_number;
    }

  public:
    WaveCalculator(Real gravity, Real water_depth, Real wave_angular_freq, Real wave_amplitude, Real wave_phase)
        : gravity_(gravity), water_depth_(water_depth), wave_angular_freq_(wave_angular_freq),
          wave_amplitude_(wave_amplitude), wave_phase_(wave_phase)
    {
        computeWaveNumber();
    }

    Real computeFreeSurfaceHeight(Real time, Real position_x = 0.0)
    {
        free_surface_height_ = wave_amplitude_ * cos(-wave_angular_freq_ * time + wave_phase_ + wave_number_ * position_x);
        return free_surface_height_;
    }

    Vecd computeEmitterVelocity(Real time, Real position_y, Real position_x = 0.0)
    {
        Vecd emitter_velocity = Vecd::Zero();
        emitter_velocity[0] = (gravity_ * wave_amplitude_ * wave_number_ / wave_angular_freq_) *
                              (cosh(wave_number_ * (position_y + water_depth_)) / cosh(wave_number_ * water_depth_)) *
                              cos(-wave_angular_freq_ * time + wave_phase_ + wave_number_ * position_x);
        return emitter_velocity;
    }
};

//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct EmitterVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;
    WaveCalculator wave_calculator;

    template <class BoundaryConditionType>
    EmitterVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(emitter_velocity), t_ref_(emitter_time),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()),
          wave_calculator(gravity_g, DH, wave_angular_freq, wave_amplitude, wave_phase) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();

        // Real current_free_surface_height = 0.0;

        // current_free_surface_height = DH - wave_calculator.computeFreeSurfaceHeight(current_time);
        // if (position[1] + halfsize_[1] > current_free_surface_height)
        // {
        //     target_velocity[0] = 0.0;
        //     target_velocity[1] = 0.0;
        // }
        // else
        // {
        //     target_velocity = wave_calculator.computeEmitterVelocity(current_time, position[1] + halfsize_[1]);
        // }
        Real u_ave = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;

        //Real normalized_y = (position[1] + halfsize_[1]) / (halfsize_[1] / 0.75);

        //target_velocity[0] = 0.3 * u_ave + (u_ave - 0.3 * u_ave) * normalized_y;
        target_velocity[0] = u_ave;

        target_velocity[1] = 0.0;

        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class EmitterBasedAcceleration : public Gravity
{
    Real target_time_;

  public:
    EmitterBasedAcceleration(Vecd target_velocity, Real target_time)
        : Gravity(target_velocity / target_time), target_time_(target_time) {}

    Vecd InducedAcceleration(const Vecd &position, Real physical_time) const
    {
        Real time_factor = physical_time / target_time_;
        Vecd acceleration = Vecd::Zero();

        if (time_factor < 1.0)
        {
            Real time_dependent_factor = 0.5 * Pi * sin(Pi * time_factor);

            Vecd base_acceleration = time_dependent_factor * Gravity::InducedAcceleration();

            // Real normalized_y = position[1] / DH;

            // Real interpolated_acceleration_x = 0.3 * base_acceleration[0] +
            //                                    (base_acceleration[0] - 0.3 * base_acceleration[0]) * normalized_y;

            //acceleration[0] = interpolated_acceleration_x;
            acceleration[0] = base_acceleration[0];
            acceleration[1] = 0.0;

            return acceleration;
        }
        else
        {
            acceleration[0] = 0.0;
            acceleration[1] = 0.0;

            return acceleration;
        }
    }
};

class CustomDampingBoundaryCondition : public fluid_dynamics::BaseFlowBoundaryCondition
{
  public:
    explicit CustomDampingBoundaryCondition(BodyRegionByCell &body_part, Real strength)
        : fluid_dynamics::BaseFlowBoundaryCondition(body_part), strength_(strength),
          damping_zone_bounds_(body_part.getBodyPartShape().getBounds()) {}

    virtual ~CustomDampingBoundaryCondition() {}

    void update(size_t index_particle_i, Real dt = 0.0)
    {
        // Real damping_factor = (pos_[index_particle_i][0] - damping_zone_bounds_.first_[0]) /
        //                       (damping_zone_bounds_.second_[0] - damping_zone_bounds_.first_[0]);

        vel_[index_particle_i][0] = 0.0;
        vel_[index_particle_i][1] = 0.0;
    }

  protected:
    Real strength_;
    BoundingBox damping_zone_bounds_;
};

//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -DH - BW), Vec2d(DL + 2.0 * BW, DH));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> emitter_disposer_buffer_particle(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(emitter_disposer_buffer_particle);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    EmitterBasedAcceleration emitter_acceleration(Vecd(emitter_velocity, 0.0), emitter_time);
    SimpleDynamics<GravityForce<EmitterBasedAcceleration>> apply_emitter_acceleration(water_block, emitter_acceleration);
    //InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> water_block_corrected_configuration(ConstructorArgs(water_block_inner, 0.5), water_wall_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    //InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_particle_indicator(water_block_inner, water_wall_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    // BodyRegionByCell damping_buffer(water_block, makeShared<MultiPolygonShape>(createDampingBufferShape()));
    // SimpleDynamics<CustomDampingBoundaryCondition> damping_wave(damping_buffer, 100.0);

    BodyAlignedBoxByCell emitter_buffer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<NonPrescribedPressure> emitter_bidirectional_buffer(emitter_buffer, emitter_disposer_buffer_particle);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<EmitterVelocity>> emitter_velocity_condition(emitter_buffer);

    BodyAlignedBoxByCell disposer_buffer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(disposer_buffer_translation)), disposer_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<NonPrescribedPressure> disposer_bidirectional_buffer(disposer_buffer, emitter_disposer_buffer_particle);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<EmitterVelocity>> disposer_velocity_condition(disposer_buffer);

    //fluid_dynamics::BidirectionalBuffer<OutletPressure> disposer_bidirectional_buffer(disposer_buffer, emitter_disposer_buffer_particle);
    //fluid_dynamics::BidirectionalBuffer<OutletPressure> disposer_bidirectional_buffer(disposer_buffer, emitter_disposer_buffer_particle);
    //SimpleDynamics<fluid_dynamics::PressureCondition<OutletPressure>> disposer_pressure_condition(disposer_buffer);

    //InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_density_by_summation(water_block_inner, water_wall_contact);
    //SimpleDynamics<fluid_dynamics::PressureCondition<NonPrescribedPressure>> emitter_pressure_condition(emitter_buffer);
    //SimpleDynamics<fluid_dynamics::PressureCondition<OutletPressure>> disposer_pressure_condition(disposer_buffer);

    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    write_body_states.addToWrite<Real>(water_block, "Pressure");               // output for debug
    write_body_states.addToWrite<int>(water_block, "Indicator");               // output for debug
    write_body_states.addToWrite<int>(water_block, "BufferParticleIndicator"); // output for debug
    ReducedQuantityRecording<TotalKineticEnergy> write_water_kinetic_energy(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_particle_indicator.exec();
    emitter_bidirectional_buffer.tag_buffer_particles.exec();
    disposer_bidirectional_buffer.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 100.0;
    Real output_interval = 0.1; /**< Time stamps for output of body states. */
    Real dt = 0.0;              /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            apply_emitter_acceleration.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();

            update_density_by_summation.exec();
            //water_block_corrected_configuration.exec();
            viscous_force.exec();
            //transport_velocity_correction.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            size_t inner_ite_dt = 0;
            while (relaxation_time < Dt)
            {
                dt = get_fluid_time_step_size.exec();
                pressure_relaxation.exec(dt);
                //kernel_summation.exec();
                //emitter_pressure_condition.exec(dt);
                //disposer_pressure_condition.exec(dt);
                density_relaxation.exec(dt);
                emitter_velocity_condition.exec(dt);
                disposer_velocity_condition.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "	Dt / dt = " << inner_ite_dt << "\n";
                write_water_kinetic_energy.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            //damping_wave.exec(Dt);
            emitter_bidirectional_buffer.injection.exec();
            //disposer_bidirectional_buffer.injection.exec();

            //emitter_bidirectional_buffer.deletion.exec();
            disposer_bidirectional_buffer.deletion.exec();
            //disposer_condition.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            boundary_particle_indicator.exec();
            emitter_bidirectional_buffer.tag_buffer_particles.exec();
            disposer_bidirectional_buffer.tag_buffer_particles.exec();
        }

        TickCount t2 = TickCount::now();
        write_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    return 0;
}
