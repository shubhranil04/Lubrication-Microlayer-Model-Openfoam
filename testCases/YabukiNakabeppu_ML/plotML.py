import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_h_vs_x.py <folder_name>")
        sys.exit(1)

    folder = sys.argv[1]
    x_path = os.path.join(folder, "xml.csv")
    h_path = os.path.join(folder, "hml.csv")

    try:
        x = np.loadtxt(x_path, delimiter=",", skiprows=2)
        h = np.loadtxt(h_path, delimiter=",", skiprows=2)
    except Exception as e:
        print(f"Error reading files: {e}")
        sys.exit(1)

    plt.figure(figsize=(8, 4))
    plt.plot(x, h, label="h(x)", color="blue")
    plt.xlabel("x")
    plt.ylabel("h(x)")
    plt.title(f"h vs x (folder: {folder})")
    plt.grid(True)
    #plt.ylim([-1e-06, 50e-06])
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
