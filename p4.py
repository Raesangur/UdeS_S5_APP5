import matplotlib.pyplot as plt
import numpy as np

def distances_calc(D0, phi, r, theta):
    D0x = D0 * np.cos(phi)
    D0y = D0 * np.sin(phi)

    rx = r * np.cos(theta + phi)
    ry = r * np.sin(theta + phi)

    Dx = D0x + rx
    Dy = D0y + ry

    D = np.sqrt(Dx**2 + Dy**2)
    angle = np.arctan(Dy / Dx)
    #D     = np.sqrt((D0 ** 2) + (2 * D0 * r * np.cos(theta)) + (r ** 2))
    #angle = phi + np.arctan(r * np.sin(theta)/(D0+r * np.cos(theta)))

    return D, angle, Dx, Dy

def _fr(r, sigma):
    if (r < 0):
        return 0
    else:
        return r / sigma * np.exp(-r**2 / (2 * sigma))

def _Fr(r, sigma):
    if (r < 0):
        return 0
    else:
        return 1 - np.exp(-r**2 / (2 * sigma))

def _r(Fr, sigma2):
#    if (Fr <= 0):
#        return 0
#    else:
    return np.sqrt(-2 * sigma2 * np.log(1 - Fr))


def monte_carlo(D0, phi_deg, sigma, N = 10000, rayleigh = False, cartesian = False):
    phi  = phi_deg * np.pi / 180
    distances = []
    angles    = []
    Dx = []
    Dy = []

    r = []
    theta = np.random.rand(10000) * np.pi * 2

    if rayleigh == True:
        r = np.random.rayleigh(sigma, N)
    else:
        r = [_r(u, sigma) for u in np.random.rand(N)]

    for i in range(N):
        mc = distances_calc(D0,  phi, theta[i], r[i])
        distances.append(mc[0])
        angles.append(mc[1] * 180 / np.pi)
        Dx.append(mc[2])
        Dy.append(mc[3])

    if cartesian == True:
        return Dx, Dy
    else:
        return distances, angles


def main():
    N     = 10000
    U     = np.random.rand(N)
    U1    = np.random.rand(N)
    theta = U * np.pi * 2

    # Displaying histogram
    fig, axs = plt.subplots(2)
    axs[0].hist(theta, align="left")[0]
    axs[1].hist(U1, align="left")[0]
    axs[0].set_title("Histogramme de la distribution de θ")
    axs[1].set_title("Histogramme de la distribution de U[0, 1]")

    sigmar = [0.25, 1, 4, 9, 16]

    # f(r)
    fig, axs = plt.subplots(3)
    fr = [[_fr(r, sigma) for r in np.linspace(0, 10, 1000)] for sigma in sigmar]
    [axs[0].plot(r) for r in fr]
    axs[0].set_title("f(r)")

    # F(r)
    Fr = [[_Fr(r, sigma) for r in np.linspace(0, 10, 1000)] for sigma in sigmar]
    [axs[1].plot(r) for r in Fr]
    axs[1].set_title("F(r)")

    # r(Fr)
    r = [[_r(Fr, sigma) for Fr in np.linspace(0, 1, 100)] for sigma in sigmar]
    [axs[2].plot(R) for R in r]
    axs[2].set_title("r(Fr)")

    # Distance and angles values, generated from inverse r and standard Rayleigh
    D4,   A4   = monte_carlo(50,  15, 4, N = N)
    D16,  A16  = monte_carlo(100, 15, 16, N = N)
    Dr4,  Ar4  = monte_carlo(50,  15, 4,  N = N, rayleigh = True)
    Dr16, Ar16 = monte_carlo(100, 15, 16, N = N, rayleigh = True)

    fig, axs = plt.subplots(2, 2)
    axs[0][0].hist(D4)
    axs[0][1].hist(D16)
    axs[1][0].hist(Dr4)
    axs[1][1].hist(Dr16)
    axs[0][0].set_title("r inverse (σ = 4)")
    axs[0][0].set_xlabel("Distance D")
    axs[0][1].set_title("r inverse (σ = 16)")
    axs[0][1].set_xlabel("Distance D")
    axs[1][0].set_title("rayleigh (σ = 4)")
    axs[1][0].set_xlabel("Distance D")
    axs[1][1].set_title("rayleigh (σ = 16)")
    axs[1][1].set_xlabel("Distance D")

    # Scatter plot of theta vs r
    fig, ax = plt.subplots()
    ax.scatter(theta, [_r(u, 4) for u in np.random.rand(N)])
    ax.set_title("Nuage de point de r en fonction de θ")

    # Scatter plot of distance vs angle
    fig, axs = plt.subplots(2, 2)
    axs[0][0].scatter(*monte_carlo(50,  15, 4, N = N))
    axs[0][0].set_title("D0 = 50, phi = 15")
    axs[0][0].set_ylabel("Angle phi")
    axs[0][0].set_xlabel("Distance D")
    axs[0][1].scatter(*monte_carlo(100, 15, 16, N = N))
    axs[0][1].set_title("D0 = 100, phi = 15")
    axs[0][1].set_ylabel("Angle phi")
    axs[0][1].set_xlabel("Distance D")
    axs[1][0].scatter(*monte_carlo(50,  30, 4, N = N))
    axs[1][0].set_title("D0 = 50, phi = 30")
    axs[1][0].set_ylabel("Angle phi")
    axs[1][0].set_xlabel("Distance D")
    axs[1][1].scatter(*monte_carlo(100, 30, 16, N = N))
    axs[1][1].set_title("D0 = 100, phi = 30")
    axs[0][1].set_ylabel("Angle phi")
    axs[0][1].set_xlabel("Distance D")

    fig, axs1 = plt.subplots(4, 2)
    fig.tight_layout()
    fig, axs2 = plt.subplots(2, 2)
    fig, axs3 = plt.subplots(4, 2)
    fig.tight_layout()

    distances = [50, 100]
    phis = [15, 30]

    for i in range(0, len(distances)):
        if distances[i] == 50:
            sigma2 = 4
        else:
            sigma2 = 16

        for j in range(0, len(phis)):
            Dx, Dy = monte_carlo(distances[i], phis[j], sigma2, N = N, cartesian = True)
            axs1[i * 2 + j][0].plot(Dx)
            axs1[i * 2 + j][1].plot(Dy)
            axs2[i][j].scatter(Dx, Dy)

            n = axs3[i * 2 + j][0].hist(Dx, align="left")
            axs3_ = axs3[i * 2 + j][0].twinx()
            axs3_.plot(n[1][:-1], n[0] / N * 100, color="red")

            n = axs3[i * 2 + j][1].hist(Dy, align="left")
            axs3_ = axs3[i * 2 + j][1].twinx()
            axs3_.plot(n[1][:-1], n[0] / N * 100, color="red")

            cov  = np.cov(Dx, Dy)
            print(cov)

            axs1[i * 2 + j][0].set_title("Dx (D0 = " + str(distances[i]) +
                                         " & phi = " + str(phis[j]) + ")")
            axs1[i * 2 + j][1].set_title("Dx (D0 = " + str(distances[i]) +
                                         " & phi = " + str(phis[j]) + ")")
            axs1[i * 2 + j][0].set_ylabel("Distance D")
            axs1[i * 2 + j][1].set_ylabel("Distance D")

            axs2[i][j].set_title("Nuage de point de [Dx, Dy]\n" +
                                 "(D0 = " + str(distances[i]) +
                                 " & phi = " + str(phis[j]) + ")")
            axs2[i][j].set_ylabel("Angle phi")
            axs2[i][j].set_xlabel("Distance D")

            axs3[i * 2 + j][0].set_title("Histogramme de Dx" + 
                                         " (D0 = " + str(distances[i]) + " & phi = " + str(phis[j]) + ")" +
                                         "\n µ = " + str(np.average(Dx)) + " & σ = " + str(np.std(Dx)))
            axs3[i * 2 + j][1].set_title("Histogramme de Dy" +
                                         " (D0 = " + str(distances[i]) + " & phi = " + str(phis[j]) + ")" +
                                         "\n µ = " + str(np.average(Dy)) + " & σ = " + str(np.std(Dy)))
            axs3[i * 2 + j][0].set_ylabel("n")
            axs3_.set_ylabel("%")



if __name__ == '__main__':
    main()
    plt.show()
