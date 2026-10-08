/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	zeroth-order_residue.h
 * @brief 	
 * @details     
 * @author Xiangyu Hu
 */

#ifndef ZEROTH_ORDER_RESIDUE_H
#define ZEROTH_ORDER_RESIDUE_H

#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
//=================================================================================================//
class CalculateAverageKGS : public LocalDynamicsReduce<ReduceSum<Real>>
{
public:
    explicit CalculateAverageKGS(SPHBody& sph_body, Real channel_length);
    virtual ~CalculateAverageKGS() {};

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
    size_t get_num_particle_in_domain(size_t input_num)
    {
        return num_particle_in_domain_ = input_num;
    }
    void output_average_kgs(Real physical_time, Real average_kgs)
    {
        const char* file_name = "average_kgs.dat";
        std::ifstream check_file(file_name, std::ios::binary | std::ios::ate);
        const bool write_header =
            !check_file.is_open() || check_file.tellg() == std::streampos(0);
        check_file.close();
        std::ofstream output_file(file_name, std::ios::app);
        if (!output_file)
        {
            throw std::runtime_error("Cannot open average_kgs.dat");
        }
        if (write_header)
        {
            output_file << "#VARIABLES = \"Time\", \"AverageKGS\"\n";
        }
        output_file << std::scientific << std::setprecision(15)
            << physical_time << "\t" << average_kgs << "\n";
    }

protected:
    Vecd* zero_gradient_residue_;
    Vecd* pos_;
    Real channel_length_;
    Real particle_spacing_;
    //
    Real mean_kgs_;
    size_t num_particle_in_domain_ = 0;
    size_t total_real_particles_;
};
//=================================================================================================//
class CalculateParticleInDomain : public LocalDynamicsReduce<ReduceSum<size_t>>
{
public:
    explicit CalculateParticleInDomain(SPHBody& sph_body, Real channel_length);
    virtual ~CalculateParticleInDomain() {};

    size_t reduce(size_t index_i, Real dt = 0.0);
    virtual size_t outputResult(size_t reduced_value) override;

protected:
    Vecd* pos_;
    Real channel_length_;
    Real particle_spacing_;
};
//=================================================================================================//
class CalculateAverageRiemannDissipation : public LocalDynamicsReduce<ReduceSum<Real>>
{
public:
    explicit CalculateAverageRiemannDissipation(SPHBody& sph_body, Real channel_length);
    virtual ~CalculateAverageRiemannDissipation() {};

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
    size_t get_num_particle_in_domain(size_t input_num)
    {
        return num_particle_in_domain_ = input_num;
    }
    void output_average_dissipation(Real physical_time, Real average_kgs)
    {
        const char* file_name = "average_riemann_dissipation.dat";
        std::ifstream check_file(file_name, std::ios::binary | std::ios::ate);
        const bool write_header =
            !check_file.is_open() || check_file.tellg() == std::streampos(0);
        check_file.close();
        std::ofstream output_file(file_name, std::ios::app);
        if (!output_file)
        {
            throw std::runtime_error("Cannot open average_riemann_dissipation.dat");
        }
        if (write_header)
        {
            output_file << "#VARIABLES = \"Time\", \"AverageDissipation\"\n";
        }
        output_file << std::scientific << std::setprecision(15)
            << physical_time << "\t" << average_kgs << "\n";
    }

protected:
    Real* dissipation_riemann_;
    Vecd* pos_;
    Real channel_length_;
    Real particle_spacing_;
    size_t num_particle_in_domain_ = 0;
};
} // namespace SPH
#endif // K_EPSILON_TURBULENT_MODEL_H