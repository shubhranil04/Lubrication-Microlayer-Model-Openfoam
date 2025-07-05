import numpy as np
import matplotlib.pyplot as plt

def main():
    x_path = "x0.csv"
    h_path = "h0.csv"

    try:
        x = np.loadtxt(x_path, delimiter=",", skiprows=2)
        h = np.loadtxt(h_path, delimiter=",", skiprows=2)
    except Exception as e:
        print(f"Error reading files: {e}")
        return

    plt.figure(figsize=(8, 4))
    plt.plot(x, h, label="h(x)", color="blue")
    plt.xlabel("x")
    plt.ylabel("h(x)")
    plt.title("h vs x")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

