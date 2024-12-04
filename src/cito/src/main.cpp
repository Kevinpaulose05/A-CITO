#include "cito/penalty_loop.h"
#include <fstream>

int main(int argc, char const *argv[])
{
    int external_force_start_step = 10;
    int external_force_end_step = 20;
    int joint_index = 2;
    double external_force_value = 10.0;

    double adm_mass = 1.0;
    double adm_damping = 5.0;
    double adm_stiffness = 50.0;

    // Define a sinusoidal force profile
    auto variable_force = [](int step, double max_force, double frequency) {
        return max_force * sin(2 * M_PI * frequency * step);
    };

    double max_force = 10.0; // Maximum force
    double frequency = 0.01; // Frequency in Hz



    const char *mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);

    YAML::Node params = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/params.yaml");
    std::string modelPathStr = paths::workspaceDir + "/src/cito/model/" + params["model"].as<std::string>();
    const char *modelPath = modelPathStr.c_str();
    mjModel *m = (strlen(modelPath) > 4 && !strcmp(modelPath + strlen(modelPath) - 4, ".mjb")) ? mj_loadModel(modelPath, NULL) : mj_loadXML(modelPath, NULL, NULL, 0);

    if (!m)
    {
        mju_error("Cannot load the model");
    }

    Params cp(m);
    Control cc(m, &cp);
    SCVX scvx(m, &cp, &cc);
    PenaltyLoop pl(m, &cp, &cc, &scvx);

    eigMd U0, U;
    U0.resize(cp.m, cp.N);
    U.resize(cp.m, cp.N);

    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/vscm.yaml");
    double kCon0 = vscm["kCon0"].as<double>();
    for (int i = 0; i < cp.N; i++)
    {
        U0.col(i).setZero();
        for (int j = 0; j < cp.nPair; j++)
        {
            U0.col(i)[m->nu + j] = kCon0;
        }
    }

    U = (params["penaltyLoop"].as<bool>()) ? pl.solve(U0) : scvx.solveSCVX(U0);

    trajectory traj = scvx.runSimulation(U, false, 2, 1);
    trajectory traj_with_forces = scvx.runSimulation(U, false, 2, 1, true, external_force_start_step, external_force_end_step, joint_index, external_force_value);
    trajectory traj_with_admittance = scvx.runSimulationWithAdmittance(U, adm_mass, adm_damping, adm_stiffness, external_force_start_step, external_force_end_step, joint_index, external_force_value);
    trajectory traj_variable_force = scvx.runVariableForceSimulation(
    U, 
    [&](int step) { return variable_force(step, max_force, frequency); },
    joint_index);


    // Open log file for writing
    std::ofstream logFile("trajectory_log.csv");
    logFile << "Time,Position_Nominal,Velocity_Nominal,Position_Without_Admittance,Velocity_Without_Admittance,"
            "Position_With_Admittance,Velocity_With_Admittance,Control_Without_Admittance,Control_With_Admittance,"
            "Forces_Applied,Position_VariableForce,Velocity_VariableForce,Control_VariableForce,Variable_Force\n";

    for (int i = 0; i <= cp.N; i++) {
        eigVd pos_nominal = traj.X.col(i).head(m->nv);
        eigVd vel_nominal = traj.X.col(i).tail(m->nv);

        eigVd pos_without_admittance = traj_with_forces.X.col(i).head(m->nv);
        eigVd vel_without_admittance = traj_with_forces.X.col(i).tail(m->nv);

        eigVd pos_with_admittance = traj_with_admittance.X.col(i).head(m->nv);
        eigVd vel_with_admittance = traj_with_admittance.X.col(i).tail(m->nv);

        eigVd pos_variable_force = traj_variable_force.X.col(i).head(m->nv);
        eigVd vel_variable_force = traj_variable_force.X.col(i).tail(m->nv);

        eigVd control_without_admittance = (i < cp.N) ? traj_with_forces.U.col(i) : Eigen::VectorXd::Zero(cp.m).eval();
        eigVd control_with_admittance = (i < cp.N) ? traj_with_admittance.U.col(i) : Eigen::VectorXd::Zero(cp.m).eval();
        eigVd control_variable_force = (i < cp.N) ? traj_variable_force.U.col(i) : Eigen::VectorXd::Zero(cp.m).eval();

        eigVd forces_applied = Eigen::VectorXd::Zero(m->nv);
        if (i >= external_force_start_step && i <= external_force_end_step) {
            forces_applied[joint_index] = external_force_value;
        }

        double variable_force_value = variable_force(i, max_force, frequency);

        // Write data to CSV
        logFile << i << ",";
        logFile << pos_nominal.transpose() << "," << vel_nominal.transpose() << ",";
        logFile << pos_without_admittance.transpose() << "," << vel_without_admittance.transpose() << ",";
        logFile << pos_with_admittance.transpose() << "," << vel_with_admittance.transpose() << ",";
        logFile << control_without_admittance.transpose() << "," << control_with_admittance.transpose() << ",";
        logFile << forces_applied.transpose() << ",";
        logFile << pos_variable_force.transpose() << "," << vel_variable_force.transpose() << ",";
        logFile << control_variable_force.transpose() << ",";
        logFile << variable_force_value << "\n";
    }



    logFile.close();

    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}
