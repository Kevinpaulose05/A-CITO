

#include "cito/penalty_loop.h"

int main(int argc, char const *argv[])
{
    int external_force_start_step = 10;
    int external_force_end_step = 20;
    int joint_index = 2;
    double external_force_value = 10.0;

    double adm_mass = 1.0;
    double adm_damping = 5.0;
    double adm_stiffness = 50.0;

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

    auto computeAndDisplayMetrics = [](const trajectory &nominal, const trajectory &test, const std::string &label, int N, int nv) {
        double total_pos_deviation = 0.0, total_vel_deviation = 0.0, max_pos_deviation = 0.0, max_vel_deviation = 0.0;
        for (int i = 0; i <= N; i++)
        {
            eigVd pos_deviation = test.X.col(i).head(nv) - nominal.X.col(i).head(nv);
            eigVd vel_deviation = test.X.col(i).tail(nv) - nominal.X.col(i).tail(nv);

            total_pos_deviation += pos_deviation.squaredNorm();
            total_vel_deviation += vel_deviation.squaredNorm();
            max_pos_deviation = std::max(max_pos_deviation, pos_deviation.norm());
            max_vel_deviation = std::max(max_vel_deviation, vel_deviation.norm());
        }

        std::cout << "\nMetrics (" << label << "):\n";
        std::cout << "- RMSE in position: " << std::sqrt(total_pos_deviation / (N + 1)) << "\n";
        std::cout << "- RMSE in velocity: " << std::sqrt(total_vel_deviation / (N + 1)) << "\n";
        std::cout << "- Max position deviation: " << max_pos_deviation << "\n";
        std::cout << "- Max velocity deviation: " << max_vel_deviation << "\n";
    };

    std::cout << "\n### Trajectory Metrics ###\n";
    computeAndDisplayMetrics(traj, traj_with_forces, "No Admittance", cp.N, m->nv);
    computeAndDisplayMetrics(traj, traj_with_admittance, "With Admittance", cp.N, m->nv);

    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}


