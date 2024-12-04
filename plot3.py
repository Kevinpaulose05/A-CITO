import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Load the CSV file
data = pd.read_csv("trajectory_log.csv")

# Extract columns for positions and convert to numpy arrays
time = data['Time'].to_numpy()
pos_nominal = data['Position_Nominal'].to_numpy()
pos_without_admittance = data['Position_Without_Admittance'].to_numpy()
pos_with_admittance = data['Position_With_Admittance'].to_numpy()
pos_variable_force = data['Position_VariableForce'].to_numpy()

# Convert string data to numerical arrays (if positions are stored as lists in strings)
def process_series(series):
    return np.array([np.array(list(map(float, x.strip("[]").split()))) for x in series])

pos_nominal = process_series(pos_nominal)
pos_without_admittance = process_series(pos_without_admittance)
pos_with_admittance = process_series(pos_with_admittance)
pos_variable_force = process_series(pos_variable_force)

# Extract x, y, z components
x_nominal, y_nominal, z_nominal = pos_nominal[:, 0], pos_nominal[:, 1], pos_nominal[:, 2]
x_without_admittance, y_without_admittance, z_without_admittance = (
    pos_without_admittance[:, 0],
    pos_without_admittance[:, 1],
    pos_without_admittance[:, 2],
)
x_with_admittance, y_with_admittance, z_with_admittance = (
    pos_with_admittance[:, 0],
    pos_with_admittance[:, 1],
    pos_with_admittance[:, 2],
)
x_variable_force, y_variable_force, z_variable_force = (
    pos_variable_force[:, 0],
    pos_variable_force[:, 1],
    pos_variable_force[:, 2],
)

# Create a 3D plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot each trajectory
ax.plot(x_nominal, y_nominal, z_nominal, label="Nominal Trajectory", color="black", linestyle="--")
# ax.plot(
#     x_without_admittance,
#     y_without_admittance,
#     z_without_admittance,
#     label="With Uniform Force Pattern (Admittance OFF)",
#     color="blue",
# )
# ax.plot(
#     x_with_admittance,
#     y_with_admittance,
#     z_with_admittance,
#     label="With Uniform Force Pattern (Admittance ON)",
#     color="orange",
# )
# ax.plot(
#     x_variable_force,
#     y_variable_force,
#     z_variable_force,
#     label="With Variable Force Pattern (Admittance ON)",
#     color="green",
#     linestyle="-.",
# )

# Set labels and title
ax.set_title("3D Trajectory of the System", fontsize=16)
ax.set_xlabel("X Position", fontsize=12)
ax.set_ylabel("Y Position", fontsize=12)
ax.set_zlabel("Z Position", fontsize=12)
ax.legend()
plt.show()
