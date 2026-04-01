#ifndef UDF_P_REFINEMENT_H
#define UDF_P_REFINEMENT_H

#include "udf_k-omega_turbulent_model.h"
#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
namespace udf
{
//=================================================================================================//
    template <typename... InteractionTypes>
    class P_refinement_GetVelocityGradient;

    template <class DataDelegationType>
    class P_refinement_GetVelocityGradient<DataDelegationType>
        : public LocalDynamics, public DataDelegationType
    {
    public:
        template <class BaseRelationType>
        explicit P_refinement_GetVelocityGradient(BaseRelationType& base_relation);
        virtual ~P_refinement_GetVelocityGradient() {};

    protected:
        Matd* velocity_gradient_only_P_;
        Vecd* k_gradient_only_P_;
        Vecd* omega_gradient_only_P_;
        //
        Real* Vol_;
        Vecd* vel_;
        int* is_near_wall_P1_;
        int* is_near_wall_P2_;
        Real* turbu_k_;
        Real* turbu_omega_;
    };

    //** Inner part *
    template <>
    class P_refinement_GetVelocityGradient<Inner<>> : public P_refinement_GetVelocityGradient<DataDelegateInner>
    {
    public:
        explicit P_refinement_GetVelocityGradient(BaseInnerRelation& inner_relation);
        virtual ~P_refinement_GetVelocityGradient() {};
        void interaction(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);

    protected:
        Matd* turbu_B_;
        Matd* B_;
    };
    using P_refinement_GetVelocityGradientInner = P_refinement_GetVelocityGradient<Inner<>>;

    //** Wall part *
    template <>
    class P_refinement_GetVelocityGradient<Contact<Wall>> : public InteractionWithWall<P_refinement_GetVelocityGradient>
    {
    public:
        explicit P_refinement_GetVelocityGradient(BaseContactRelation& contact_relation);
        virtual ~P_refinement_GetVelocityGradient() {};
        void interaction(size_t index_i, Real dt = 0.0);

    protected:

    };

    //** Interface part *
    using P_refinement_GetVelocityGradientComplex = ComplexInteraction<P_refinement_GetVelocityGradient<Inner<>, Contact<Wall>>>;
//=================================================================================================//
    class P_refinement : public LocalDynamics, public kOmega_BaseTurbuClosureCoeff, public WallFunctionCoefficient
    {
    public:
        explicit P_refinement(SPHBody& sph_body);
        virtual ~P_refinement(){};

        void update(size_t index_i, Real dt = 0.0);
        Vec6d solve_1D_sublayer_Neumann(double kinematic_viscosity, double u_p_outer, double k_p_outer,
            double w_p_outer, double vel_grad_p_outer, double nut_p_outer, double h_sublayer, 
            double utau_outer, double Q_target, double k_grad_p_outer, double w_grad_p_outer, double& vel_nodeO, double& vel_nodeUM);
        Vec6d solve_1D_sublayer_Dirichlet(double kinematic_viscosity, double u_p_outer, double k_p_outer,
            double w_p_outer, double vel_grad_p_outer, double nut_p_outer, double h_sublayer,
            double utau_outer, double Q_target, double k_grad_p_outer, double w_grad_p_outer, double& vel_nodeO, double& vel_nodeUM);
        void tdma(int N, const double* a, const double* b, const double* c, const double* d, double* x);
        void tdma5(const double a[5], const double b[5], const double c[5], const double d[5], double x[5]);
        void tdma10(const double a[10], const double b[10], const double c[10], const double d[10], double x[10]);

        inline Real get_loacal_flow_rate(Real average_flow_rate_over_particle_P, Real SPH_vel_grad_P, Real nodeO_U, Real dp)
        {
            Real nodeOS_U = nodeO_U + SPH_vel_grad_P * 0.5 * dp;
            Real flow_rate_half = (nodeO_U + nodeOS_U) * dp / 4.0;
            Real flow_rate_whole = average_flow_rate_over_particle_P;
            Real flow_rate_local = flow_rate_whole - flow_rate_half;
            return flow_rate_local;
        }

        inline Real obtainTangentialComponent(const Vecd& vec, const Vecd& normal)
        {
            Real dot_un = vec.dot(normal);                  
            Real norm_sqr = vec.dot(vec) - dot_un * dot_un;   
            return std::sqrt(std::max(0.0, norm_sqr));          
        }

        void writeTecplotFromVec6d(
            const Vec6d& node_val,   // node_value_[index_i]
            double U_nodeO,
            double U_nodeUM,
            double distance_to_wall,
            int num_sub_node_,       // = 5
            int NF,
            int index
        )
        {
            int ny = num_sub_node_; // 5
            // ================== 提取数据 ==================
            double utau = node_val[0];
            std::vector<double> U;
            std::vector<double> y;
            U.reserve(ny + 2);
            y.reserve(ny + 2);
            // ================== 构造 y ==================
            double hy = distance_to_wall / double(num_sub_node_);
            double y_p = 0.5 * hy;
            // ===== 原有5个节点 =====
            for (int i = 0; i < ny; ++i) {
                y.push_back(y_p + i * hy);
                U.push_back(node_val[i + 1]);
            }
            // ===== nodeO =====
            double y_nodeO = y.back() + 0.5 * hy;
            y.push_back(y_nodeO);
            U.push_back(U_nodeO);
            // ===== nodeUM =====
            double y_nodeUM = y.back() + 0.5 * hy; // 再往上 0.5hy → 总共 +hy
            y.push_back(y_nodeUM);
            U.push_back(U_nodeUM);
            int n_total = y.size(); // = 7
            std::vector<double> K(n_total, 0.0);
            std::vector<double> OMEGA(n_total, 0.0);
            std::vector<double> NUT(n_total, 0.0);
            std::string header_line = "ZONE T=\"SPH(1D)-NODE NF="
                + std::to_string(NF)
                + " ("
                + std::to_string(index)
                + ")\"";
            std::string filename = "pipe_node_nf"
                + std::to_string(NF) + "_"
                + std::to_string(index) + ".dat";
            std::ofstream fout(filename);
            fout << "$VARIABLES = \"Y\", \"U\", \"K\", \"OMEGA\", \"NUT(k/omega)\"\n";
            fout << "$friction velocity = " << utau << "\n";
            fout << header_line << "\n";
            fout << std::scientific << std::setprecision(8);
            for (int i = 0; i < n_total; ++i) {
                fout << y[i] << " "
                    << U[i] << " "
                    << K[i] << " "
                    << OMEGA[i] << " "
                    << NUT[i] << "\n";
            }
            fout.close();
            std::cout << "Tecplot (with nodeO & nodeUM) created: " << filename << std::endl;
        }
        void writeTecplotFromVec6dDirichlet(
            const Vec6d& node_val,   // node_value_[index_i]
            double U_nodeO,
            double U_nodeUM,
            double distance_to_wall,
            int num_sub_node_,       // = 5
            int NF,
            int index
        )
        {
            int ny = num_sub_node_; // 5
            // ================== 提取数据 ==================
            double utau = node_val[0];
            std::vector<double> U;
            std::vector<double> y;
            U.reserve(ny + 2);
            y.reserve(ny + 2);
            // ================== 构造 y ==================
            double hy = distance_to_wall / (double(num_sub_node_) + 0.5);
            double y_p = 0.5 * hy;
            // ===== 原有5个节点 =====
            for (int i = 0; i < ny; ++i) {
                y.push_back(y_p + i * hy);
                U.push_back(node_val[i + 1]);
            }
            // ===== nodeO =====
            double y_nodeO = y.back() + 1.0 * hy;
            y.push_back(y_nodeO);
            U.push_back(U_nodeO);
            // ===== nodeUM =====
            double y_nodeUM = y.back() + 0.0 * hy; // nodeO nodeUM overlaps 
            y.push_back(y_nodeUM);
            U.push_back(U_nodeUM);
            int n_total = y.size(); // = 7
            std::vector<double> K(n_total, 0.0);
            std::vector<double> OMEGA(n_total, 0.0);
            std::vector<double> NUT(n_total, 0.0);
            std::string header_line = "ZONE T=\"SPH(1D)-NODE NF="
                + std::to_string(NF)
                + " ("
                + std::to_string(index)
                + ")\"";
            std::string filename = "pipe_node_nf"
                + std::to_string(NF) + "_"
                + std::to_string(index) + ".dat";
            std::ofstream fout(filename);
            fout << "$VARIABLES = \"Y\", \"U\", \"K\", \"OMEGA\", \"NUT(k/omega)\"\n";
            fout << "$friction velocity = " << utau << "\n";
            fout << header_line << "\n";
            fout << std::scientific << std::setprecision(8);
            for (int i = 0; i < n_total; ++i) {
                fout << y[i] << " "
                    << U[i] << " "
                    << K[i] << " "
                    << OMEGA[i] << " "
                    << NUT[i] << "\n";
            }
            fout.close();
            std::cout << "Tecplot (with nodeO & nodeUM) created: " << filename << std::endl;
        }

    protected:
        int num_sub_node_;
        Real* friction_velocity_from_sublayer_;
        Real* target_flow_rate_in_sublayer_;
        Real* vel_ps_magnitude_;
        Real* dudn_for_local_flow_rate_;
        Real* utau_node_;
        Vec6d* node_value_; // ** Temporary treatment only valid for 5-node configuration, first is utau, then velocity *
        Real* dUdn_P_sublayer_magnitude_;
        Matd* dUdn_P_sublayer_;
        Real* vel_nodeO_;
        Real* vel_nodeUM_;
        Real* dUdn_P_nodeU_;
        Real* global_flow_rate_over_P_;
        Real* half_flow_rate_over_P_;
        //
        int* is_near_wall_P1_;
        Real* y_p_;
        Vecd* vel_;
        Real* turbu_k_;
        Real* turbu_omega_;
        Real* rho_;
        Viscosity& viscosity_;
        Real mu_;
        Matd* velocity_gradient_only_P_;
        Real* turbu_mu_;
        Real* distance_to_dummy_interface_;
        Real* wall_shear_stress_;
        Vecd* e_nearest_normal_;
        Real fluid_particle_spacing_;
        Real* physical_time_;
        Matd* velocity_gradient_;
        Vecd* k_gradient_only_P_;
        Vecd* omega_gradient_only_P_;
    };
//=================================================================================================//
    class BodyStatesRecordingToVtpIncludeNode : public BodyStatesRecording
    {
    public:
        BodyStatesRecordingToVtpIncludeNode(SPHSystem& sph_system) : BodyStatesRecording(sph_system) {};
        virtual ~BodyStatesRecordingToVtpIncludeNode() {};

    protected:
        virtual void writeWithFileName(const std::string& sequence) override;

        //** Apart from the one line P refinment statement, others are directly copied from Library *
        template <typename OutStreamType>
        void writeParticlesToVtk(OutStreamType& output_stream, BaseParticles& particles)
        {
            size_t total_real_particles = particles.TotalRealParticles();
            ParticleVariables& variables_to_write = particles.VariablesToWrite();

            // write sorted particles ID
            output_stream
                << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
            output_stream << "    ";
            for (size_t i = 0; i != total_real_particles; ++i)
            {
                output_stream << i << " ";
            }
            output_stream << std::endl;
            output_stream << "    </DataArray>\n";

            // write particle IDs
            constexpr int type_index_UnsignedInt = DataTypeIndex<UnsignedInt>::value;
            for (DiscreteVariable<UnsignedInt>* variable : std::get<type_index_UnsignedInt>(variables_to_write))
            {
                UnsignedInt* data_field = variable->Data();
                output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Int32\" Format=\"ascii\">\n";
                output_stream << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    output_stream << std::fixed << std::setprecision(9) << data_field[i] << " ";
                }
                output_stream << std::endl;
                output_stream << "    </DataArray>\n";
            }

            // write integers
            constexpr int type_index_int = DataTypeIndex<int>::value;
            for (DiscreteVariable<int>* variable : std::get<type_index_int>(variables_to_write))
            {
                int* data_field = variable->Data();
                output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Int32\" Format=\"ascii\">\n";
                output_stream << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    output_stream << std::fixed << std::setprecision(9) << data_field[i] << " ";
                }
                output_stream << std::endl;
                output_stream << "    </DataArray>\n";
            }

            // write scalars
            constexpr int type_index_Real = DataTypeIndex<Real>::value;
            for (DiscreteVariable<Real>* variable : std::get<type_index_Real>(variables_to_write))
            {
                Real* data_field = variable->Data();
                output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\" Format=\"ascii\">\n";
                output_stream << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    output_stream << std::fixed << std::setprecision(9) << data_field[i] << " ";
                }
                output_stream << std::endl;
                output_stream << "    </DataArray>\n";
            }

            // write vectors
            constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
            for (DiscreteVariable<Vecd>* variable : std::get<type_index_Vecd>(variables_to_write))
            {
                Vecd* data_field = variable->Data();
                output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
                output_stream << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    Vec3d vector_value = upgradeToVec3d(data_field[i]);
                    output_stream << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
                }
                output_stream << std::endl;
                output_stream << "    </DataArray>\n";
            }

            //** P refinement *
            writeNodeValueToVtk(output_stream, particles);

            // write matrices
            constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
            for (DiscreteVariable<Matd>* variable : std::get<type_index_Matd>(variables_to_write))
            {
                Matd* data_field = variable->Data();
                output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type= \"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
                output_stream << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    Mat3d matrix_value = upgradeToMat3d(data_field[i]);
                    for (int k = 0; k != 3; ++k)
                    {
                        Vec3d col_vector = matrix_value.col(k);
                        output_stream << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
                    }
                }
                output_stream << std::endl;
                output_stream << "    </DataArray>\n";
            }
        };
        template <typename OutStreamType>
        void writeNodeValueToVtk(OutStreamType& output_stream, BaseParticles& particles)
        {

            size_t total_real_particles = particles.TotalRealParticles();
            ParticleVariables& variables_to_write = particles.VariablesToWrite();

            // write vectors
            constexpr int type_index_Vecd = DataTypeIndex<Vec6d>::value;
            for (DiscreteVariable<Vec6d>* variable : std::get<type_index_Vecd>(variables_to_write))
            {
                Vec6d* data_field = variable->Data();
                output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"6\" Format=\"ascii\">\n";
                output_stream << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    Vec6d vector_value = upgradeToVec6d(data_field[i]);
                    output_stream << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2]
                        << " " << vector_value[3] << " " << vector_value[4] << " " << vector_value[5] << " ";
                }
                output_stream << std::endl;
                output_stream << "    </DataArray>\n";
            }
        };
        inline Vec6d upgradeToVec6d(const Vec6d& input)
        {
            return input;
        };
     };
//=================================================================================================//
} // namespace udf
} // namespace fluid_dynamics
} // namespace SPH
#endif // UDF_COMMON_TURBULENCE_MODEL_H