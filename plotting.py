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
norm_forces_applied = compute_norm(forces_applied)

# Compute trajectory deviations
deviation_without_admittance = np.linalg.norm(pos_without_admittance - pos_nominal, axis=1)
deviation_with_admittance = np.linalg.norm(pos_with_admittance - pos_nominal, axis=1)

# Compute stability metrics: convergence or divergence of deviations
stability_without_admittance = deviation_without_admittance / norm_pos_nominal
stability_with_admittance = deviation_with_admittance / norm_pos_nominal

# Compute energy consumption
energy_without_admittance = [u**2 for u in norm_control_without_admittance]
energy_with_admittance = [u**2 for u in norm_control_with_admittance]

# Plot position deviations
plt.figure(figsize=(12, 6))
plt.plot(time, norm_pos_nominal, label="Nominal")
plt.plot(time, norm_pos_without_admittance, label="Without Admittance")
plt.plot(time, norm_pos_with_admittance, label="With Admittance")
plt.title("Position Deviation Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Position Norm")
plt.legend()
plt.grid()
plt.show()

# Plot velocity deviations
plt.figure(figsize=(12, 6))
plt.plot(time, norm_vel_nominal, label="Nominal")
plt.plot(time, norm_vel_without_admittance, label="Without Admittance")
plt.plot(time, norm_vel_with_admittance, label="With Admittance")
plt.title("Velocity Deviation Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Velocity Norm")
plt.legend()
plt.grid()
plt.show()

# Plot control efforts
plt.figure(figsize=(12, 6))
plt.plot(time[:-1], norm_control_without_admittance[:-1], label="Control Without Admittance", linestyle="--", alpha=0.8, marker='o')
plt.plot(time[:-1], norm_control_with_admittance[:-1], label="Control With Admittance", linestyle=":", alpha=0.8, marker='x')
plt.title("Control Effort Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Control Effort Norm")
plt.legend()
plt.grid()
plt.show()

# Plot system stability
plt.figure(figsize=(12, 6))
plt.plot(time, stability_without_admittance, label="Stability Without Admittance", linestyle="--", alpha=0.8, marker='o')
plt.plot(time, stability_with_admittance, label="Stability With Admittance", linestyle=":", alpha=0.8, marker='x')
plt.title("System Stability Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Deviation Ratio (Deviation / Nominal)")
plt.legend()
plt.grid()
plt.show()

# Plot trajectory deviations
plt.figure(figsize=(12, 6))
plt.plot(time, deviation_without_admittance, label="Deviation Without Admittance")
plt.plot(time, deviation_with_admittance, label="Deviation With Admittance")
plt.title("Trajectory Deviation from Nominal Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Trajectory Deviation Norm")
plt.legend()
plt.grid()
plt.show()

# Plot energy consumption
plt.figure(figsize=(12, 6))
plt.plot(time[:-1], energy_without_admittance[:-1], label="Energy Without Admittance", linestyle="--", alpha=0.8)
plt.plot(time[:-1], energy_with_admittance[:-1], label="Energy With Admittance", linestyle=":", alpha=0.8)
plt.title("Energy Consumption Over Time")
plt.xlabel("Time (Steps)")
plt.ylabel("Energy (Control Effort Squared)")
plt.legend()
plt.grid()
plt.show()
