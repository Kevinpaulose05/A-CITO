#include "cito/scvx.h"

// ***** CONSTRUCTOR ***********************************************************
SCVX::SCVX(const mjModel *m_, Params *cp_, Control *cc_) : m(m_), cp(cp_), cc(cc_), nd(m_, cp_, cc_), sq(m_, cp_)
{
    // initialize Eigen variables
    finalPos.resize(6);
    // get SCVX parameters
    YAML::Node paramSCVX = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/scvx.yaml");
    beta_expand = paramSCVX["beta_expand"].as<double>();
    beta_shrink = paramSCVX["beta_shrink"].as<double>();
    maxIter = paramSCVX["maxIter"].as<int>();
    dLTol = paramSCVX["dLTol"].as<double>();
    rho0 = paramSCVX["rho0"].as<double>();
    rho1 = paramSCVX["rho1"].as<double>();
    rho2 = paramSCVX["rho2"].as<double>();
    rMin = paramSCVX["rMin"].as<double>();
    rMax = paramSCVX["rMax"].as<double>();
    r0 = paramSCVX["r0"].as<double>();
    J = new double[maxIter + 1];
    JTemp = new double[maxIter + 1];
    JTilde = new double[maxIter + 1];
    dJ = new double[maxIter + 1];
    dL = new double[maxIter + 1];
    rho = new double[maxIter + 1];
    r = new double[maxIter + 1];
    accept = new bool[maxIter];
    // set initial trust-region radius
    r[0] = r0;
    // set bounds
    cc->getBounds();
    // resize trajectories
    XSucc.resize(cp->n, cp->N + 1);
    dX.resize(cp->n, cp->N + 1);
    XTilde.resize(cp->n, cp->N + 1);
    USucc.resize(cp->m, cp->N);
    UTemp.resize(cp->m, cp->N);
    dU.resize(cp->m, cp->N);
    Fx.resize(cp->N);
    Fu.resize(cp->N);
    for (int i = 0; i < cp->N; i++)
    {
        Fx[i].resize(cp->n, cp->n);
        Fu[i].resize(cp->n, cp->m);
    }
}
// ***** DESTRUCTOR ************************************************************
SCVX::~SCVX()
{
    delete[] J;
    delete[] JTemp;
    delete[] JTilde;
    delete[] dJ;
    delete[] dL;
    delete[] rho;
    delete[] r;
    delete[] accept;
}

// ***** FUNCTIONS *************************************************************
// getCost: returns the nonlinear cost given control trajectory and final state
double SCVX::getCost(const eigMd &X, const eigMd &U)
{
    // final cost
    finalPos = X.col(cp->N).segment(cp->controlJointDOF0, 6);
    Jf = 0.5 * (cp->weight[0] * (cp->desiredPos.head(2) - finalPos.head(2)).squaredNorm() +
                cp->weight[1] * (cp->desiredPos.tail(4) - finalPos.tail(4)).squaredNorm());
    // integrated cost
    Ji = cp->weight[2] * X.leftCols(cp->N).bottomRows(m->nv).squaredNorm() +
         cp->weight[3] * U.bottomRows(cp->nPair).sum();
    // total cost
    Jt = Jf + Ji;
    return Jt;
}

// runSimulation: rolls-out and linearizes the dynamics given control trajectory
trajectory SCVX::runSimulation(
    const eigMd &U,
    bool linearize,
    int save,
    double compensateBias,
    bool applyExternalForces,
    int external_force_start_step,
    int external_force_end_step,
    int joint_index,
    double external_force_value)
{
    // Create mjData
    mjData *d = NULL;
    d = mj_makeData(m);

    // Initialize d
    mju_copy(d->qpos, m->key_qpos, m->nq);
    mj_forward(m, d);
    cc->setControl(d, U.col(0), compensateBias);

    // Rollout (and linearize) the dynamics
    for (int i = 0; i < cp->N; i++)
    {
        mj_forward(m, d);

        // Apply external forces if enabled
        if (applyExternalForces)
        {
            // Apply external forces at specific time steps
            if (i >= external_force_start_step && i <= external_force_end_step)
            {
                // Apply external force to specific joint
                d->qfrc_applied[joint_index] += external_force_value;
            }
            else
            {
                // Reset external forces when not applied
                d->qfrc_applied[joint_index] = 0.0;
            }
        }

        // Get the current state values
        XSucc.col(i).setZero();
        XSucc.col(i) = cp->getState(d);

        // Linearization
        if (linearize)
        {
            Fx[i].setZero();
            Fu[i].setZero();
            nd.linDyn(d, U.col(i), Fx[i].data(), Fu[i].data(), compensateBias);
        }

        // Take tc/dt steps
        cc->takeStep(d, U.col(i), save, compensateBias);
    }

    XSucc.col(cp->N).setZero();
    XSucc.col(cp->N) = cp->getState(d);

    // Delete data
    mj_deleteData(d);

    // Build trajectory
    traj.X = XSucc;
    traj.U = U;
    if (linearize)
    {
        traj.Fx = Fx;
        traj.Fu = Fu;
    }
    return traj;
}

// applyAdmittanceControl: Apply admittance control based on dynamics
void SCVX::applyAdmittanceControl(mjData *d, int joint_index, double external_force, double adm_mass, double adm_damping, double adm_stiffness)
{
    double pos_error = d->qpos[joint_index] - 0.0; // Desired position (assumed 0)
    double vel_error = d->qvel[joint_index] - 0.0; // Desired velocity (assumed 0)

    double corrective_torque = (-adm_stiffness * pos_error) +
                               (-adm_damping * vel_error) +
                               (external_force / adm_mass);

    d->qfrc_applied[joint_index] += corrective_torque;
}

// runSimulationWithAdmittance: Run simulation with admittance control
trajectory SCVX::runSimulationWithAdmittance(
    const eigMd &U,
    double adm_mass,
    double adm_damping,
    double adm_stiffness,
    int external_force_start_step,
    int external_force_end_step,
    int joint_index,
    double external_force_value)
{
    mjData *d = mj_makeData(m);

    mju_copy(d->qpos, m->key_qpos, m->nq);
    mju_copy(d->qvel, m->key_qvel, m->nv);
    mj_forward(m, d);

    for (int i = 0; i < cp->N; i++)
    {
        mj_forward(m, d);

        if (i >= external_force_start_step && i <= external_force_end_step)
        {
            applyAdmittanceControl(d, joint_index, external_force_value, adm_mass, adm_damping, adm_stiffness);
        }

        XSucc.col(i) = cp->getState(d);
        cc->takeStep(d, U.col(i), 1, true);
    }

    XSucc.col(cp->N) = cp->getState(d);
    mj_deleteData(d);

    traj.X = XSucc;
    traj.U = U;

    return traj;
}

trajectory SCVX::runVariableForceSimulation(
    const eigMd &U,
    std::function<double(int)> forceProfile,
    int joint_index)
{
    // Create mjData
    mjData *d = mj_makeData(m);

    // Initialize data
    mju_copy(d->qpos, m->key_qpos, m->nq);
    mju_copy(d->qvel, m->key_qvel, m->nv);
    mj_forward(m, d);

    // Rollout the dynamics with variable forces
    for (int i = 0; i < cp->N; i++) {
        mj_forward(m, d);

        // Apply variable force
        double variable_force = forceProfile(i);
        d->qfrc_applied[joint_index] = variable_force;

        // Record the state
        XSucc.col(i) = cp->getState(d);

        // Take one simulation step
        cc->takeStep(d, U.col(i), 1, true);
    }

    // Record the final state
    XSucc.col(cp->N) = cp->getState(d);

    // Clean up
    mj_deleteData(d);

    // Build trajectory
    traj.X = XSucc;
    traj.U = U;

    return traj;
}



// solveSCVX: executes the successive convexification algorithm
eigMd SCVX::solveSCVX(const eigMd &U0)
{
    // initialize USucc for the first succession
    USucc = U0;
    // start the SCVX algorithm
    int iter = 0;
    while (!stop)
    {
        std::cout << "Iteration " << iter + 1 << ":" << '\n';
        // simulation and convexification ======================================
        if (iter == 0 || accept[iter - 1])
        {
            std::cout << "INFO: convexification in progress\n";
            auto tDiffStart = std::chrono::system_clock::now();
            trajS = {};
            trajS = this->runSimulation(USucc, true, 0, 1);
            auto tDiffEnd = std::chrono::system_clock::now();
            std::cout << "INFO: convexification took " << std::chrono::duration<double>(tDiffEnd - tDiffStart).count() << " s \n";
        }
        // get the nonlinear cost if the first iteration
        if (iter == 0)
        {
            J[iter] = this->getCost(trajS.X, USucc);
        }
        // convex optimization =================================================
        double *dTraj = new double[cp->nTraj];
        std::cout << "INFO: QP solver in progress\n\n";
        auto tQPStart = std::chrono::system_clock::now();
        sq.solveCvx(dTraj, r[iter], trajS.X, USucc, trajS.Fx, trajS.Fu, cc->isJFree, cc->isAFree,
                    cc->qposLB, cc->qposUB, cc->tauLB, cc->tauUB);
        auto tQPEnd = std::chrono::system_clock::now();
        std::cout << "\nINFO: QP solver took " << std::chrono::duration<double>(tQPEnd - tQPStart).count() << " s \n\n";
        // apply the change
        for (int i = 0; i < cp->N + 1; i++)
        {
            // states
            dX.col(i).setZero();
            XTilde.col(i).setZero();
            for (int j = 0; j < cp->n; j++)
            {
                dX.col(i)[j] = dTraj[i * cp->n + j];
            }
            XTilde.col(i) = trajS.X.col(i) + dX.col(i);
            // controls
            if (i < cp->N)
            {
                dU.col(i).setZero();
                UTemp.col(i).setZero();
                for (int j = 0; j < cp->m; j++)
                {
                    dU.col(i)[j] = dTraj[(cp->N + 1) * cp->n + i * cp->m + j];
                }
                UTemp.col(i) = USucc.col(i) + dU.col(i);
            }
        }
        // evaluate the dynamics for the change and get the cost values ========
        trajTemp = {};
        trajTemp = this->runSimulation(UTemp, false, 0, 1);
        // get the linear and nonlinear costs
        JTilde[iter] = this->getCost(XTilde, UTemp);
        JTemp[iter] = this->getCost(trajTemp.X, UTemp);
        // similarity measure ==================================================
        dJ[iter] = J[iter] - JTemp[iter];
        dL[iter] = J[iter] - JTilde[iter];
        rho[iter] = dJ[iter] / dL[iter];
        if (fabs(dL[iter]) < dLTol)
        {
            dLTolMet = 1;
        }
        // accept or reject the solution =======================================
        // reject
        if (rho[iter] <= rho0 || (dL[iter] < 0 && dJ[iter] < 0))
        {
            accept[iter] = false;
            r[iter + 1] = r[iter] / beta_shrink;
            J[iter + 1] = J[iter];
        }
        else
        {
            accept[iter] = true;
        }
        // accept
        if (accept[iter])
        {
            J[iter + 1] = JTemp[iter];
            USucc = UTemp;
            if (rho[iter] < rho1)
            {
                r[iter + 1] = r[iter] / beta_shrink;
            }
            else if (rho[iter] >= rho1 && rho[iter] < rho2)
            {
                r[iter + 1] = r[iter];
            }
            else if (rho[iter] >= rho2)
            {
                r[iter + 1] = r[iter] * beta_expand;
            }
        }
        // bound the trust region radius r
        r[iter + 1] = std::max(r[iter + 1], rMin);
        r[iter + 1] = std::min(r[iter + 1], rMax);
        // stopping criteria check =============================================
        if (iter + 1 == maxIter)
        {
            stop = true;
            std::cout << "\n\n\tINFO: Maximum number of iterations reached.\n\n";
        }
        if (dLTolMet)
        {
            stop = true;
            std::cout << "\n\n\tINFO: |dL| = |" << dL[iter] << "| < dLTol = " << dLTol << "\n\n";
        }
        // screen output for the iteration =====================================
        std::cout << "Actual:\nFinal pos: " << trajTemp.X.col(cp->N).head(m->nv).transpose() << "\n";
        std::cout << "Final vel: " << trajTemp.X.col(cp->N).tail(m->nv).transpose() << "\n";
        std::cout << "Predicted:\nFinal pos: " << XTilde.col(cp->N).head(m->nv).transpose() << "\n";
        std::cout << "Final vel: " << XTilde.col(cp->N).tail(m->nv).transpose() << "\n";
        std::cout << "L = " << JTilde[iter] << ", J = " << JTemp[iter] << ", kmax = " << trajTemp.U.bottomRows(cp->nPair).maxCoeff() << ", kavg = " << trajTemp.U.bottomRows(cp->nPair).sum() / (cp->nPair * cp->N) << "\n\n\n";
        // next iteration ======================================================
        iter++;
        delete[] dTraj;
    }
    // summary screen output ===============================================
    std::cout << "\n\nSCVX Summary\nJ0=" << J[0] << "\n\n";
    for (int i = 0; i < iter; i++)
    {
        if (i % 10 == 0)
        {
            printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n",
                   "Iteration", "L", "J", "dL", "dJ", "rho", "r", "accept");
        }
        printf("%-12d%-12.6g%-12.6g%-12.3g%-12.3g%-12.3g%-12.3g%-12d\n",
               i + 1, JTilde[i], JTemp[i], dL[i], dJ[i], rho[i], r[i], accept[i]);
    }
    return USucc;
}

// refresh: refreshes SCVX variables for a new run
void SCVX::refresh()
{
    // delete old variables
    delete[] J;
    delete[] JTemp;
    delete[] JTilde;
    delete[] dJ;
    delete[] dL;
    delete[] rho;
    delete[] r;
    delete[] accept;
    // create new variables
    J = new double[maxIter + 1];
    JTemp = new double[maxIter + 1];
    JTilde = new double[maxIter + 1];
    dJ = new double[maxIter + 1];
    dL = new double[maxIter + 1];
    rho = new double[maxIter + 1];
    r = new double[maxIter + 1];
    accept = new bool[maxIter];
    // set initial trust-region radius
    r[0] = r0;
    // reset flags
    stop = false;
    dLTolMet = false;
}