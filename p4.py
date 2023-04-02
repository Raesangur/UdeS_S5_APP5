import matplotlib.pyplot as plt
import numpy as np

def monte_carlo(D0, phi, r, theta):
    D0x = D0 * np.cos(phi)
    D0y = D0 * np.sin(phi)

    rx = r * np.cos(theta)
    ry = r * np.sin(theta)

    Dx = D0x + rx
    Dy = D0y + ry

    return np.sqrt(Dx**2 + Dy**2)

def _fr(r, sigma):
    if (r < 0):
        return 0
    else:
        return r / sigma * np.exp(-r**2 / (2 * sigma))

def _Fr(r, sigma):
    if (r < 0):
        return 0
    else:
        return -1 * np.exp(-r**2 / (2 * sigma))

def _r(Fr, sigma):
    if (Fr <= 0):
        return 0
    else:
        return 2 * sigma * np.log(Fr)

def main():
    N     = 10000
    U     = np.random.rand(N)
    theta = U * np.pi * 2

    # Displaying histogram
    plt.hist(U, align="left")[0]
    plt.show()

    # Calculating r
    sigmar = [0.25, 1, 4, 9, 16]

    fr = [[_fr(r, sigma) for r in np.linspace(0, 10, 1000)] for sigma in sigmar]
    [plt.plot(r) for r in fr]
    plt.show()

    Fr = [[_Fr(r, sigma) for r in np.linspace(0, 10, 1000)] for sigma in sigmar]
    [plt.plot(r) for r in Fr]
    plt.show()

    r = [[_r(Fr, sigma) for Fr in np.linspace(0, 1, 100)] for sigma in sigmar]
    [plt.plot(R) for R in r]
    plt.show()

    Ur = np.random.rayleigh(4, N)
    D  = []
    Dr = []
    phi = 15 * np.pi / 180
    for i in range(N):
        D. append(monte_carlo(50, phi, theta[i], _r(U[i], 4)))
        Dr.append(monte_carlo(50, phi, theta[i], Ur[i]))

    fig, axs = plt.subplots(2)
    axs[0].hist(D)
    axs[1].hist(Dr)
    axs[0].set_title("Distribution avec fonction r inverse")
    axs[1].set_title("Distribution avec fonction rayleigh")
    plt.show()

if __name__ == '__main__':
    main()