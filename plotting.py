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

# Convert string data to numerical arrays
def process_series(series):
    return np.array([np.array(list(map(float, x.strip("[]").split()))) for x in series])

pos_nominal = process_series(pos_nominal)
vel_nominal = process_series(vel_nominal)
pos_without_admittance = process_series(pos_without_admittance)
vel_without_admittance = process_series(vel_without_admittance)
pos_with_admittance = process_series(pos_with_admittance)
vel_with_admittance = process_series(vel_with_admittance)

# Compute norms for deviations
def compute_norm(series):
    return [np.linalg.norm(v) for v in series]

norm_pos_nominal = compute_norm(pos_nominal)
norm_vel_nominal = compute_norm(vel_nominal)
norm_pos_without_admittance = compute_norm(pos_without_admittance)
norm_vel_without_admittance = compute_norm(vel_without_admittance)
norm_pos_with_admittance = compute_norm(pos_with_admittance)
norm_vel_with_admittance = compute_norm(vel_with_admittance)

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
