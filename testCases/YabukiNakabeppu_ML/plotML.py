import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_microlayer.py <folder_name>")
        sys.exit(1)

    folder = sys.argv[1]
    x_file = os.path.join(folder, "xml.csv")
    h_file = os.path.join(folder, "hml.csv")

    try:
        x = np.loadtxt(x_file, skiprows=2)
        h = np.loadtxt(h_file, skiprows=2)
    except Exception as e:
        print(f"Error reading files: {e}")
        sys.exit(1)

    if x.shape != h.shape:
        raise ValueError("Shape mismatch between x and h")

    plt.figure(figsize=(8, 5))
    plt.plot(x, h, '-o', label='Microlayer thickness')
    plt.xlabel("x (m)")
    plt.ylabel("h (m)")
    plt.title(f"Microlayer thickness vs x ({folder})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
