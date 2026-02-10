#include "19.h"
using namespace SPH;

int main(int ac, char *av[])
{
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);

    /** Restart. */
    bool is_write_restart_file = true;
    int restart_output_interval = 500;
    sph_system.setRestartStep(0); //% SPH

    /** Average. */
    bool is_write_average_contour_file = true;
    Real time_start_average_data = 80.0; //% Average, make sure time span is large engouth to achieve steady 
    Real time_output_contour_average_data = 90.0; //% Average
    int num_output_contour_average_file_limit = 40;
    Real magnify_ratio_avergae_contour = 10.0;

    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);

    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    /**
     * @brief Material property, particles and body creation of fluid.
     */

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    //water_block.defineBodyLevelSetShape();
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticlesWithReserve<BaseParticles, Reload>(inlet_particle_buffer, water_block.getName())
        : water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody observer_center_point(sph_system, "ObserverCenterPoint");
    observer_center_point.generateParticles<ObserverParticles>(observer_location_center_point);

    get_observation_locations();
    output_observer_theoretical_y();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_location);

    observe_nearwall::getObservingLineLengthAndEndPoints();
    observe_nearwall::getPositionsOfMultipleObserveLines();
    observe_nearwall::output_observe_positions();
    observe_nearwall::output_observe_theoretical_x();
    observe_nearwall::output_number_observe_points_on_lines();
    ObserverBody friction_velocity_observer(sph_system, "NearwallFrictionVelocityObserver");
    friction_velocity_observer.generateParticles<ObserverParticles>(observe_nearwall::observation_locations);

    ObserverBody observer_body(sph_system, makeShared<WaterBlock>("ObserverBody")); //% Average
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? observer_body.generateParticles<BaseParticles, Reload>(water_block.getName())
        : observer_body.generateParticles<BaseParticles, Lattice>();

    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    ContactRelation observer_centerpoint_contact(observer_center_point, {&water_block});
    ContactRelation friction_velocity_observer_contact(friction_velocity_observer, {&water_block});
    ContactRelation fluid_observer_contact2(observer_body, {&water_block}); //% Average
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
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(wall_boundary);
        BodyStatesRecordingToVtp write_inserted_body_to_vtp_water(water_block);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        ReloadParticleIO write_particle_reload_files_water(water_block);
        /** A  Physics relaxation step. */
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner_water(water_block_inner);
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
    /** For pressure outlet . */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

    /** Turbulent standard wall function needs normal vectors of wall. */
    //NearShapeSurface near_surface(water_block, makeShared<WallBoundary>("Wall"));

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(water_block_inner, water_wall_contact);

    //InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_fluid(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::udf::TurbulentLinearGradientCorrectionMatrixInner> corrected_configuration_fluid_only_inner(water_block_inner);

    /** Pressure relaxation algorithm with Riemann solver for viscous flows. */
    //Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionForOpenBoundaryFlowWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);

    /** Density relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerNoRiemann> density_relaxation(water_block_inner);
    //InteractionDynamics<fluid_dynamics::Integration2ndHalfOnlyWallAcousticRiemannAdjusted> density_relaxation_wall(water_wall_contact);
    //density_relaxation.post_processes_.push_back(&density_relaxation_wall);

    /** Turbulent.Note: When use wall function, K Epsilon calculation only consider inner */
    InteractionWithUpdate<fluid_dynamics::udf::JudgeIsNearWall> update_near_wall_status(water_block_inner, water_wall_contact, y_p_constant);

    //InteractionWithUpdate<fluid_dynamics::kOmega_GetVelocityGradientInner> get_velocity_gradient(water_block_inner, weight_vel_grad_sub_nearwall);
    InteractionWithUpdate<fluid_dynamics::udf::kOmega_GetVelocityGradientComplex> get_velocity_gradient(water_block_inner, water_wall_contact);


    SimpleDynamics<fluid_dynamics::udf::kOmega_kTransportEquationInner> k_equation_relaxation(water_block_inner, initial_turbu_values, is_AMRD, is_blended);
    InteractionDynamics<fluid_dynamics::udf::kOmega_TKE_Diffusion> compute_TKE_diffusion(water_block_inner);
    SimpleDynamics<fluid_dynamics::udf::kOmega_omegaTransportEquationInner> epsilon_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::udf::kOmega_TSDR_Diffusion_and_Gradient_Dot_Inner> compute_TSDR_diffusion_and_gradient_k_omega(water_block_inner);

    InteractionDynamics<fluid_dynamics::udf::TKEnergyForceComplex> turbulent_kinetic_energy_force(water_block_inner, water_wall_contact);
    InteractionDynamics<fluid_dynamics::udf::kOmega_WallFunctionCorrection> standard_wall_function_correction(water_block_inner, water_wall_contact);

    SimpleDynamics<fluid_dynamics::udf::ConstrainNormalVelocityInRegionP> constrain_normal_velocity_in_P_region(water_block);

    /** Choose one, ordinary or turbulent. Computing viscous force, */
    InteractionWithUpdate<fluid_dynamics::udf::TurbulentViscousForceWithWall> turbulent_viscous_force(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);

    /** Impose transport velocity, with or without limiter . */
    //InteractionWithUpdate<fluid_dynamics::TransportVelocityLimitedCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::udf::TVC_ModifiedLimited_RKGC_OBFCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);

    /** A temporarily test for the limiter . */
    SimpleDynamics<fluid_dynamics::udf::GetLimiterOfTransportVelocityCorrection> get_limiter_of_transport_velocity_correction(water_block, 1000);

    /** Evaluation of density by summation approach. */
    //InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);

    /** Initialize particle acceleration. */
    StartupAcceleration time_dependent_acceleration(Vec2d(U_f, 0.0), 2.0);
    SimpleDynamics<GravityForce<StartupAcceleration>> apply_gravity_force(water_block, time_dependent_acceleration);

    //----------------------------------------------------------------------
    // Left/Inlet buffer
    //----------------------------------------------------------------------
    AlignedBox left_emitter_shape(xAxis, Transform(Vec2d(left_buffer_translation)), left_buffer_halfsize);
    AlignedBoxByCell left_emitter(water_block, left_emitter_shape);
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, inlet_particle_buffer);

    //SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);

    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);

    /** Turbulent kEpsilon_InflowTurbulentCondition.It needs characteristic Length to calculate turbulent length  */
    SimpleDynamics<fluid_dynamics::udf::kOmega_InflowTurbulentCondition> impose_turbulent_inflow_condition(left_emitter, characteristic_length, relaxation_rate_turbulent_inlet, type_turbulent_inlet);

    //----------------------------------------------------------------------
    // Right/Outlet buffer
    //----------------------------------------------------------------------
    AlignedBox right_emitter_shape(xAxis, Transform(Rotation2d(Pi), Vec2d(right_buffer_translation)), right_buffer_halfsize);
    AlignedBoxByCell right_emitter(water_block, right_emitter_shape);
    fluid_dynamics::BidirectionalBuffer<RightOutflowPressure> right_bidirection_buffer(right_emitter, inlet_particle_buffer);

    //SimpleDynamics<fluid_dynamics::PressureCondition<RightOutflowPressure>> right_outflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<RightOutflowPressure>> right_outflow_pressure_condition(right_emitter);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density_pressure(water_block_inner, water_wall_contact);
    SimpleDynamics<UpdateVolume> update_volume(water_block);

    /** Choose one, ordinary or turbulent. Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::udf::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    //ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);

    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    /** Turbulent eddy viscosity calculation needs values of Wall Y start. */
    SimpleDynamics<fluid_dynamics::udf::kOmegaTurbulentEddyViscosity> update_eddy_viscosity(water_block);
    
    ObservingAQuantity<Real> observing_pressure(fluid_observer_contact2, "Pressure");          //% Average pressure
    SimpleDynamics<ParticleSnapshotAverage<Real>> average_pressure(observer_body, "Pressure"); //% Average pressure
    
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    
    /** Restart. */
    RestartIO restart_io(sph_system);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");            // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator");            // output for debug
    body_states_recording.addToWrite<Real>(water_block, "Density");             // output for debug
    body_states_recording.addToWrite<Vecd>(water_block, "KernelGradientIntegral"); // output for debug
    ObservedQuantityRecording<Vecd> write_recorded_water_velocity("Velocity", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_k("TurbulenceKineticEnergy", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_mut("TurbulentViscosity", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_omega("TurbulentSpecificDissipation", fluid_observer_contact);
    body_states_recording.addToWrite<int>(water_block, "BufferIndicator");
    //RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_centerpoint_quantity("TurbulentViscosity", observer_centerpoint_contact);
    ObservedQuantityRecording<Real> write_nearwall_friction_velocity("WallShearStress", friction_velocity_observer_contact);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");

    BodyStatesRecordingToVtp write_observation_states(observer_body);     //% Average
    write_observation_states.addToWrite<Real>(observer_body, "Pressure"); //% Average pressure

    /**
     * @brief Setup geometry and initial conditions.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");

    /** Restart. */
    if (sph_system.RestartStep() != 0)
    {
        physical_time = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedList();
        water_block_complex.updateConfiguration();
        observer_centerpoint_contact.updateConfiguration();
        fluid_observer_contact.updateConfiguration();
        friction_velocity_observer_contact.updateConfiguration();
        fluid_observer_contact2.updateConfiguration(); //** Average *
    }
    size_t number_of_iterations = sph_system.RestartStep();

    int screen_output_interval = 100;
    //int observation_sample_interval = screen_output_interval * 2;

    int num_output_contour_average_file = 0;  //** Average *

    Real end_time = 100.0;                      /**< End time. */
    Real cutoff_ratio = 0.9;                    //** cutoff_time should be a integral and the same as the PY script */
    Real cutoff_time = end_time * cutoff_ratio; //** cutoff_time should be a integral and the same as the PY script */
    
    Real num_output_files = 40.0 * (is_write_average_contour_file ? magnify_ratio_avergae_contour : 1.0);  //** Average *
    
    Real Output_Time = end_time / num_output_files; /**< Time stamps for output of body states. */
    Real index_check_file_fully_developed = num_output_files * cutoff_ratio;

    Real dt = 0.0;                      /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;

    //----------------------------------------------------------------------
    //	Preparation, if use restart, better to fullfill
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec();
    /** Tag inlet/outlet truncated particles */
    inlet_outlet_surface_particle_indicator.exec();
    /** Tag in/outlet buffer particles */
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    update_near_wall_status.exec();
    corrected_configuration_fluid.exec();
    corrected_configuration_fluid_only_inner.exec();
    get_velocity_gradient.exec();
    update_eddy_viscosity.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    //write_centerpoint_quantity.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    int num_output_file = 0;
    std::ofstream logfile("output/output.log");
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            apply_gravity_force.exec();

            //Real Dt = get_fluid_advection_time_step_size.exec();
            Real Dt = get_turbulent_fluid_advection_time_step_size.exec();

            //inlet_outlet_surface_particle_indicator.exec();

            //update_density_by_summation.exec();
            update_fluid_density_pressure.exec();
            //** This is to address the bug in density summation *
            update_volume.exec();

            corrected_configuration_fluid.exec();
            corrected_configuration_fluid_only_inner.exec();

            if (physical_time > turbulent_module_activate_time) //** A temporary treatment *
            {
                update_eddy_viscosity.exec();
            }

            //viscous_force.exec();
            turbulent_viscous_force.exec();

            if (physical_time > turbulent_module_activate_time) //** A temporary treatment *
            {
                get_velocity_gradient.exec();
                compute_TKE_diffusion.exec();
                compute_TSDR_diffusion_and_gradient_k_omega.exec();
                update_near_wall_status.exec();
                standard_wall_function_correction.exec();
            }

            transport_velocity_correction.exec();

            kernel_summation.exec();

            get_limiter_of_transport_velocity_correction.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            int inner_itr = 0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                if (physical_time > turbulent_module_activate_time) //** A temporary treatment *
                {
                    turbulent_kinetic_energy_force.exec();
                }

                pressure_relaxation.exec(dt);

                left_inflow_pressure_condition.exec(dt);
                right_outflow_pressure_condition.exec(dt);

                if (is_constrain_normal_velocity_in_P_region)
                    constrain_normal_velocity_in_P_region.exec();

                inflow_velocity_condition.exec();

                if (physical_time > turbulent_module_activate_time) //** A temporary treatment *
                {
                    impose_turbulent_inflow_condition.exec();
                }

                density_relaxation.exec(dt);

                if (physical_time > turbulent_module_activate_time) //** A temporary treatment *
                {
                    k_equation_relaxation.exec(dt);
                    epsilon_equation_relaxation.exec(dt);
                }

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_itr++;
                //std::cout << "num_output_file=" << num_output_file << std::endl;
                //if (GlobalStaticVariables::physical_time_ >9.3)
                //{
                //body_states_recording.writeToFile();
                //}
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                //if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                //{
                //    write_centerpoint_quantity.writeToFile(number_of_iterations);
                //}
                logfile << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                        << physical_time
                        << "	Dt = " << Dt << "	dt = " << dt << std::endl;
            }
            /** Restart. */
            if (is_write_restart_file)
            {
                if (number_of_iterations % restart_output_interval == 0)
                {
                    restart_io.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            // ** First do injection for all buffers *
            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();
            // ** Then do deletion for all buffers *
            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();

            /** Update cell linked list and configuration. */
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            friction_velocity_observer_contact.updateConfiguration();

            /** Tag truncated inlet/outlet particles*/
            inlet_outlet_surface_particle_indicator.exec();
            /** Tag in/outlet buffer particles that suffer pressure condition*/
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();

            if (physical_time > cutoff_time)
            {
                write_recorded_water_velocity.writeToFile(number_of_iterations);
                write_recorded_water_k.writeToFile(number_of_iterations);
                write_recorded_water_mut.writeToFile(number_of_iterations);
                write_recorded_water_omega.writeToFile(number_of_iterations);
                write_nearwall_friction_velocity.writeToFile(number_of_iterations);
            }
            //if (GlobalStaticVariables::physical_time_ > end_time * 0.5)
            //body_states_recording.writeToFile();
            if (is_write_average_contour_file) //** Average *
            {
                if (physical_time > time_start_average_data)
                {
                    fluid_observer_contact2.updateConfiguration(); //** Average *
                    //% Average pressure
                    observing_pressure.exec();
                    average_pressure.exec();
                }
            }
        }
        //TickCount t2 = TickCount::now();
        if (!is_write_average_contour_file)  //** Average *
        {
            body_states_recording.writeToFile();
        }
        observer_centerpoint_contact.updateConfiguration();
        num_output_file++;
        //if (num_output_file == 100)
        //    system("pause");
        //TickCount t3 = TickCount::now();
        if (is_write_average_contour_file) //** Average *
        {
            if (physical_time > time_output_contour_average_data)
            {
                if (num_output_contour_average_file < num_output_contour_average_file_limit)
                {
                    fluid_observer_contact2.updateConfiguration(); //% Average
                    //% Average pressure
                    observing_pressure.exec();
                    average_pressure.exec();
                    write_observation_states.writeToFile(); //% Average
                    num_output_contour_average_file++;
                    if (num_output_contour_average_file == num_output_contour_average_file_limit)
                    {
                        std::cout << "Finish outputing average contour files " << std::endl;
                        system("pause");
                    }
                }
            }
        }
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << "Cutoff_time: " << cutoff_time
              << " seconds." << std::endl;
    std::cout << "For checking fully-developed or not, index of the cutoff output file =  " << index_check_file_fully_developed << std::endl;
    logfile << "Total wall time for computation: " << tt.seconds()
            << " seconds." << std::endl;
    logfile.close();
    //if (sph_system.GenerateRegressionData())
    //{
    //    write_centerpoint_quantity.generateDataBase(1.0e-3);
    //}
    //else
    //{
    //    write_centerpoint_quantity.testResult();
    //}
    return 0;
}
