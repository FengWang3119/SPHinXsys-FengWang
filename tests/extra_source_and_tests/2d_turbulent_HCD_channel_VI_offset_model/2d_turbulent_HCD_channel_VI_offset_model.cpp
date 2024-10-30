#include "2d_turbulent_HCD_channel_VI_offset_model.h"
using namespace SPH;

int main(int ac, char *av[])
{
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);

    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(true);

    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    IOEnvironment io_environment(sph_system);
    /**
     * @brief Material property, particles and body creation of fluid.
     */

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape();
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    //water_block.generateParticles<ParticleGeneratorLattice>();
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    //water_block.generateParticlesWithReserve<Lattice>(inlet_particle_buffer);
    //water_block.generateParticles<Lattice>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticlesWithReserve<Reload>(inlet_particle_buffer, water_block.getName())
        : water_block.generateParticlesWithReserve<Lattice>(inlet_particle_buffer);
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<Lattice>();

    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        using namespace relax_dynamics;
        /** body topology only for particle relaxation */
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(wall_boundary);
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles_water(water_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp({&wall_boundary});
        BodyStatesRecordingToVtp write_inserted_body_to_vtp_water({&water_block});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        ReloadParticleIO write_particle_reload_files_water(water_block);
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);
        relax_dynamics::RelaxationStepLevelSetCorrectionInner relaxation_step_inner_water(water_block_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        random_inserted_body_particles_water.exec(0.25);

        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_inner_water.SurfaceBounding().exec();

        write_inserted_body_to_vtp.writeToFile(0);
        write_inserted_body_to_vtp_water.writeToFile(0);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            relaxation_step_inner_water.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
                write_inserted_body_to_vtp_water.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of the wall_boundary finish !" << std::endl;
        std::cout << "The physics relaxation process of the water_block finish !" << std::endl;

        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        write_particle_reload_files_water.writeToFile(0);
        return 0;
    }

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** Turbulent standard wall function needs normal vectors of wall. */
    NearShapeSurface near_surface(water_block, makeShared<WallBoundary>("Wall"));

    /** Pressure relaxation algorithm with Riemann solver for viscous flows. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    /** Density relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::TurbuIntegration2ndHalfWithWallDissipativeRiemann> density_relaxation(water_block_inner, water_wall_contact);

    /** Turbulent.Note: When use wall function, K Epsilon calculation only consider inner */
    InteractionWithUpdate<fluid_dynamics::JudgeIsNearWall> update_near_wall_status(water_block_inner, water_wall_contact, near_surface);

    InteractionDynamics<fluid_dynamics::GetVelocityGradientInner> get_velocity_gradient(water_block_inner);
    //InteractionDynamics<fluid_dynamics::GetVelocityGradientComplex> get_velocity_gradient(water_block_inner, water_wall_contact);

    /** Turbulent.Note: Temporarily transfer parameters at this place. The 3rd parameter refers to extra dissipation for viscous */
    InteractionWithUpdate<fluid_dynamics::K_TurbulentModelInner> k_equation_relaxation(water_block_inner, initial_turbu_values, 1);
    InteractionWithUpdate<fluid_dynamics::E_TurbulentModelInner> epsilon_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::TKEnergyForceComplex> turbulent_kinetic_energy_force(water_block_inner, water_wall_contact);
    InteractionDynamics<fluid_dynamics::StandardWallFunctionCorrection> standard_wall_function_correction(water_block_inner, water_wall_contact, y_p_constant);

    /** Turbulent eddy viscosity calculation needs values of Wall Y start. */
    SimpleDynamics<fluid_dynamics::ConstrainNormalVelocityInRegionP> constrain_normal_velocity_in_P_region(water_block);

    SimpleDynamics<fluid_dynamics::GetTimeAverageCrossSectionData> get_time_average_cross_section_data(water_block_inner, num_observer_points, monitoring_bound, offset_distance);

    /** Choose one, ordinary or turbulent. Computing viscous force, */
    InteractionWithUpdate<fluid_dynamics::TurbulentViscousForceWithWall> turbulent_viscous_force(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);

    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);

    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);

    water_block.addBodyStateForRecording<Real>("Pressure"); // output for debug
    water_block.addBodyStateForRecording<int>("Indicator"); // output for debug
    water_block.addBodyStateForRecording<Real>("Density");  // output for debug

    /** Initialize particle acceleration. */
    TimeDependentAcceleration time_dependent_acceleration(Vec2d::Zero());
    SimpleDynamics<GravityForce> apply_gravity_force(water_block, time_dependent_acceleration);

    BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, inlet_particle_buffer, xAxis);
    BodyAlignedBoxByCell inlet_velocity_buffer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inlet_velocity_buffer_inflow_condition(inlet_velocity_buffer, 1.0);

    BodyAlignedBoxByCell disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, 0);

    /** Turbulent InflowTurbulentCondition.It needs characteristic Length to calculate turbulent length  */
    SimpleDynamics<fluid_dynamics::InflowTurbulentCondition> impose_turbulent_inflow_condition(inlet_velocity_buffer, characteristic_length, 0.8);

    /** Choose one, ordinary or turbulent. Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    //ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);

    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    /** Turbulent eddy viscosity calculation needs values of Wall Y start. */
    SimpleDynamics<fluid_dynamics::TurbulentEddyViscosity> update_eddy_viscosity(water_block);

    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);
    /**
     * @brief Setup geometry and initial conditions.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");

    /** Output the start states of bodies. */
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 200.0;              /**< End time. */
    Real Output_Time = end_time / 40.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                      /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    int num_output_file = 0;
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            apply_gravity_force.exec();

            //Real Dt = get_fluid_advection_time_step_size.exec();
            Real Dt = get_turbulent_fluid_advection_time_step_size.exec();

            inlet_outlet_surface_particle_indicator.exec();

            update_density_by_summation.exec();

            update_eddy_viscosity.exec();

            //viscous_force.exec();
            turbulent_viscous_force.exec();

            transport_velocity_correction.exec();
            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            int inner_itr = 0;
            while (relaxation_time < Dt)
            {
                //if (inner_itr == 8)
                //{
                //    std::cin.get();
                //}

                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                turbulent_kinetic_energy_force.exec();

                pressure_relaxation.exec(dt);

                constrain_normal_velocity_in_P_region.exec();

                inlet_velocity_buffer_inflow_condition.exec();

                impose_turbulent_inflow_condition.exec();

                density_relaxation.exec(dt);

                update_near_wall_status.exec();
                get_velocity_gradient.exec(dt);
                standard_wall_function_correction.exec();
                k_equation_relaxation.exec(dt);
                epsilon_equation_relaxation.exec(dt);
                k_equation_relaxation.update_prior_turbulent_value();

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                inner_itr++;
                //std::cout << "num_output_file=" << num_output_file << std::endl;
                //if (GlobalStaticVariables::physical_time_ >9.3)
                //{
                //    body_states_recording.writeToFile();
                //}
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** inflow injection*/
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();

            if (GlobalStaticVariables::physical_time_ > end_time * 0.5)
            {
                get_time_average_cross_section_data.exec();
                get_time_average_cross_section_data.output_time_history_data(end_time * 0.75);
            }
            //if (GlobalStaticVariables::physical_time_ > end_time * 10.0)
            //{
            //    body_states_recording.writeToFile();
            //    num_output_file++;
            //}
        }
        //TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        num_output_file++;
        if (num_output_file == 440)
            std::cin.get();
        //TickCount t3 = TickCount::now();
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    get_time_average_cross_section_data.get_time_average_data(end_time * 0.75);
    std::cout << "The time-average data is output " << std::endl;
    return 0;
}
