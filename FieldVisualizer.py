import matplotlib as plt
import numpy as np
import matplotlib.animation as plta
import pandas as pd

# Function to read data from a CSV file
def read_data(filename):
    df = pd.read_csv(filename, sep=",", header=None)
    return df

def main():
    # Read data from files
    data1 = read_data("H.txt")
    data2 = read_data("E.txt")
    print(data1)
    print(data2)
    """
    # Create figure and axes
    fig, ax = plt.subplots()
    line1, = ax.plot([], [], label='E Field')
    line2, = ax.plot([], [], label='H Field')
    lines = [line1, line2]

    # Setting the axes limits
    ax.set_xlim(0, max(data1.shape[1], data2.shape[1]) - 1)
    ax.set_ylim(min(np.min(data1), np.min(data2)), max(np.max(data1), np.max(data2)))

    # Adding legend
    ax.legend()

    # Creating the animation
    ani = plta.FuncAnimation(fig, update, frames=min(len(data1), len(data2)),
                                  fargs=(lines, data1, data2), blit=True)

    # Show plot
    plt.show()
    """
main()