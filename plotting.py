import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV file
data = pd.read_csv("trajectory_log.csv")

# Extract columns and convert to numpy arrays
time = data['Time'].to_numpy()
pos_nominal = data['Position_Nominal'].to_numpy()
vel_nominal = data['Velocity_Nominal'].to_numpy()
pos_without_admittance = data['Position_Without_Admittance'].to_numpy()
vel_without_admittance = data['Velocity_Without_Admittance'].to_numpy()
pos_with_admittance = data['Position_With_Admittance'].to_numpy()
vel_with_admittance = data['Velocity_With_Admittance'].to_numpy()
control_without_admittance = data['Control_Without_Admittance'].to_numpy()
control_with_admittance = data['Control_With_Admittance'].to_numpy()
forces_applied = data['Forces_Applied'].to_numpy()
pos_variable_force = data['Position_VariableForce'].to_numpy()
vel_variable_force = data['Velocity_VariableForce'].to_numpy()
control_variable_force = data['Control_VariableForce'].to_numpy()
variable_force = data['Variable_Force'].to_numpy()

# Convert string data to numerical arrays
def process_series(series):
    return np.array([np.array(list(map(float, x.strip("[]").split()))) for x in series])

pos_nominal = process_series(pos_nominal)
vel_nominal = process_series(vel_nominal)
pos_without_admittance = process_series(pos_without_admittance)
vel_without_admittance = process_series(vel_without_admittance)
pos_with_admittance = process_series(pos_with_admittance)
vel_with_admittance = process_series(vel_with_admittance)
control_without_admittance = process_series(control_without_admittance)
control_with_admittance = process_series(control_with_admittance)
forces_applied = process_series(forces_applied)
pos_variable_force = process_series(pos_variable_force)
vel_variable_force = process_series(vel_variable_force)
control_variable_force = process_series(control_variable_force)

# Compute norms for deviations
def compute_norm(series):
    return [np.linalg.norm(v) for v in series]

norm_pos_nominal = compute_norm(pos_nominal)
norm_vel_nominal = compute_norm(vel_nominal)
norm_pos_without_admittance = compute_norm(pos_without_admittance)
norm_vel_without_admittance = compute_norm(vel_without_admittance)
norm_pos_with_admittance = compute_norm(pos_with_admittance)
norm_vel_with_admittance = compute_norm(vel_with_admittance)
norm_control_without_admittance = compute_norm(control_without_admittance)
norm_control_with_admittance = compute_norm(control_with_admittance)
norm_pos_variable_force = compute_norm(pos_variable_force)
norm_vel_variable_force = compute_norm(vel_variable_force)
norm_control_variable_force = compute_norm(control_variable_force)

# Compute trajectory deviations
deviation_without_admittance = np.linalg.norm(pos_without_admittance - pos_nominal, axis=1)
deviation_with_admittance = np.linalg.norm(pos_with_admittance - pos_nominal, axis=1)
deviation_variable_force = np.linalg.norm(pos_variable_force - pos_nominal, axis=1)

# Compute cumulative control effort
cumulative_control_without_admittance = np.cumsum(norm_control_without_admittance)
cumulative_control_with_admittance = np.cumsum(norm_control_with_admittance)
cumulative_control_variable_force = np.cumsum(norm_control_variable_force)

# Stability Metric
stability_without_admittance = deviation_without_admittance / norm_pos_nominal
stability_with_admittance = deviation_with_admittance / norm_pos_nominal
stability_variable_force = deviation_variable_force / norm_pos_nominal


# Plot position deviations
plt.figure(figsize=(12, 6))
plt.plot(time, norm_pos_nominal, label="Nominal")
plt.plot(time, norm_pos_without_admittance, label="With Uniform Force Pattern (Admittance OFF)")
plt.plot(time, norm_pos_with_admittance, label="With Uniform Force Pattern (Admittance ON)")
plt.plot(time, norm_pos_variable_force, label="With Variable Force Pattern (Admittance ON)", linestyle="--", alpha=0.8)
plt.title("Position Deviation Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Position Norm")
plt.legend()
plt.grid()
plt.show()

# Plot velocity deviations
plt.figure(figsize=(12, 6))
plt.plot(time, norm_vel_nominal, label="Nominal")
plt.plot(time, norm_vel_without_admittance, label="With Uniform Force Pattern (Admittance OFF)")
plt.plot(time, norm_vel_with_admittance, label="With Uniform Force Pattern (Admittance ON)")
plt.plot(time, norm_vel_variable_force, label="With Variable Force Pattern (Admittance ON)", linestyle="--", alpha=0.8)
plt.title("Velocity Deviation Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Velocity Norm")
plt.legend()
plt.grid()
plt.show()

# Plot trajectory deviations
plt.figure(figsize=(12, 6))
plt.plot(time, deviation_without_admittance, label="Deviation With Uniform Force Pattern (Admittance OFF)")
plt.plot(time, deviation_with_admittance, label="Deviation With Uniform Force Pattern (Admittance ON)")
plt.plot(time, deviation_variable_force, label="Deviation With Variable Force Pattern (Admittance ON)", linestyle="--")
plt.title("Trajectory Deviation from Nominal Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Trajectory Deviation Norm")
plt.legend()
plt.grid()
plt.show()

# Stability Metric graph
plt.figure(figsize=(12, 6))
plt.plot(time, stability_without_admittance, label="Stability With Uniform Force Pattern (Admittance OFF)", linestyle="--")
plt.plot(time, stability_with_admittance, label="Stability With Uniform Force Pattern (Admittance ON)", linestyle=":")
plt.plot(time, stability_variable_force, label="Stability With Variable Force Pattern (Admittance ON)", linestyle="-")
plt.title("System Stability Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Deviation Ratio")
plt.legend()
plt.grid()
plt.show()

# Plot cumulative control effort
plt.figure(figsize=(12, 6))
plt.plot(time, cumulative_control_without_admittance, label="Control With Uniform Force Pattern (Admittance OFF)")
plt.plot(time, cumulative_control_with_admittance, label="Control With Uniform Force Pattern (Admittance ON)")
plt.plot(time, cumulative_control_variable_force, label="Control With Variable Force Pattern (Admittance ON)", linestyle="--")
plt.title("Cumulative Control Effort Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Cumulative Control Effort")
plt.legend()
plt.grid()
plt.show()

# Plot force over time
plt.figure(figsize=(12, 6))
plt.plot(time, variable_force, label="Variable Force", color="purple")
plt.title("Variable Force Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Force Magnitude")
plt.legend()
plt.grid()
plt.show()