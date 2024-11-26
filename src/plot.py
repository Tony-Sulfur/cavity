import numpy as np
import matplotlib.pyplot as plt

# Load data, handling irregular spacing
data = np.genfromtxt('plot.txt', delimiter=None, invalid_raise=False)

# Filter out rows with invalid or missing data
data = data[~np.isnan(data).any(axis=1)]

# Split data into x and y columns
x = data[:, 0]
y = data[:, 1]

# Create a plot
plt.figure(figsize=(8, 6))
plt.plot(x, y, linestyle='-', color='b', label='y vs. x')

# Add labels, title, and legend
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Plot from Irregular .txt Data')
plt.legend()

# Show the plot
plt.grid(True)
plt.show()
