import pandas as pd
import matplotlib.pyplot as plt
import time

# Function to plot rows from two dataframes
def plot_rows(df1, df2):
    plt.ion()  # Enable interactive mode
    fig, ax = plt.subplots()

    for i in range(len(df1)):
        # Clear the plot for the next iteration
        ax.clear()
        
        # Plot the rows up to the current index
        ax.plot(df1.iloc[i], label='File 1')
        ax.plot(df2.iloc[i], label='File 2')
        
        # Add legend and labels
        ax.legend()
        ax.set_xlabel('Time')
        ax.set_ylabel('Values')
        ax.set_title('Row-wise Plotting of CSV Data')
        
        # Draw and pause for a brief moment
        plt.draw()
        plt.pause(0.1)
    
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