import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

maxY = 1e-13

# Function to plot rows from two dataframes
def plot_rows(df1, df2):
    plt.ion()  # Enable interactive mode
    fig, ax = plt.subplots()

    for i in range(len(df1)):
        # Clear the plot for the next iteration
        ax.clear()
        
        # Plot the rows at the current index
        ax.plot(df1.iloc[i], label='Electric field')
        ax.plot(df2.iloc[i], label='Magnetic field')
        
        # Add legend and labels
        ax.legend()
        ax.set_xlabel('z-axis')
        ax.set_ylabel('Field strength')
        #ax.set_xticks(np.linspace(0,183,20, endpoint=True))
        ax.set_ylim(bottom=-maxY, top=maxY)
        ax.set_title('E and H fields')
        
        # Draw and pause for a brief moment
        plt.draw()
        plt.pause(0.05)
    
    plt.ioff()  # Disable interactive mode
    plt.show()

# Read CSV files
file1 = 'E.txt'
file2 = 'H.txt'

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Ensure both files have the same length
if len(df1) != len(df2):
    raise ValueError("The two CSV files must have the same number of rows.")

# Plot rows
plot_rows(df1, df2)
