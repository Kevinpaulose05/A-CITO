#ifndef ADMITTANCE_CONTROL_H
#define ADMITTANCE_CONTROL_H

#include "mujoco.h"

class AdmittanceControl
{
public:
    AdmittanceControl(const mjModel *m_, double mass, double damping, double stiffness)
        : m(m_), virtual_mass(mass), virtual_damping(damping), virtual_stiffness(stiffness) {}

    void applyControl(mjData *d, int joint_index, double external_force)
    {
        double pos_error = d->qpos[joint_index] - desired_position[joint_index];
        double vel_error = d->qvel[joint_index] - desired_velocity[joint_index];

        // Admittance control law
        double corrective_torque = (-virtual_stiffness * pos_error) +
                                   (-virtual_damping * vel_error) +
                                   (external_force / virtual_mass);

        d->qfrc_applied[joint_index] += corrective_torque;
    }

    void setDesiredState(const eigVd &pos, const eigVd &vel)
    {
        desired_position = pos;
        desired_velocity = vel;
    }

private:
    const mjModel *m;
    double virtual_mass;
    double virtual_damping;
    double virtual_stiffness;
    eigVd desired_position;
    eigVd desired_velocity;
};

#endif
