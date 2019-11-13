#!/usr/bin/env python3

# Este codigo resuelve el problema de Schrödinger dependiente del tiempo
# con el metodo mostrado en el libro "Basic Concepts in Computational Physics"
# de Benjamin A. Stickler y Ewald Schachinger

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import ImageMagickWriter

masa = 1
hbarra = 1


def main():
    xi = 0.0
    xf = 500.0
    ti = 0.0
    tf = 300.0
    delta_x = 0.4 # 0.4
    delta_t = 0.15 # 0.1

    nx = int(np.rint((xf - xi) / delta_x))
    nt = int(np.rint((tf - ti) / delta_t))

    x = np.linspace(xi, xf, nx)
    v = np.zeros(nx)
    phi = np.zeros([nt, nx], dtype=np.complex_)

    # Parametros del problema, Potencial y Condicion Inicial
    parametros = [5.0, 100.0, 15.0]
    for k in range(0, nx-1 + 1):
        v[k] = potencial(x[k])
        phi[0] = condicion_inicial(x, parametros)

    # Condiciones de bordes impuestas
    phi[0][0] = 0.0
    phi[0][nx-1] = 0.0
    evolucionar(phi, nx, nt, delta_x, delta_t, v)

    # phi2 = np.absolute(phi[nt-1])
    # plt.plot(x, phi2)

    fig, ax = plt.subplots()
    ln, = plt.plot([], [])

    def init():
        ax.set_xlim(xi, xf)
        ax.set_ylim(0, 1)
        return ln,

    def update(frame):
        xdata = x
        ydata = np.absolute(phi[frame])
        # ydata = ydata / np.sqrt(np.sum(ydata))
        ln.set_data(xdata, ydata)
        return ln,

    ani = FuncAnimation(fig, update, frames=range(0, nt-1, int((nt-1)/100)),
                        init_func=init, blit=True)
    plt.show()


def evolucionar(phi, nx, nt, delta_x, delta_t, V):
    # Constantes para resumir calculos
    alpha = 2.0*masa*delta_x**2.0/(delta_t*hbarra)
    print(str(alpha/2.0))
    beta = masa*delta_x**2.0/hbarra**2.0

    # Inicializar variables
    a = np.zeros(nx, dtype=np.complex_)
    omega = np.zeros([nt, nx], dtype=np.complex_)
    b = np.zeros([nt, nx], dtype=np.complex_)

    # Empieza Algoritmo
    a[1] = 2.0*(1.0 + beta * V[1] - 1j*alpha)
    for k in range(2, nx-2 + 1):
        a[k] = 2.0*(1.0 + beta * V[k] - 1j*alpha) - 1.0/a[k-1]
        #
    for n in range(0, nt-1):
        # phi_norm = np.absolute(phi[n])
        # phi[n] = phi[n]/np.sum(phi_norm)
        normalizar(phi[n])
        for k in range(1, nx-2):
            omega[n][k] = -phi[n][k-1] + 2.0*(1j*alpha + 1.0 + beta*V[k])*phi[n][k] - phi[n][k+1]
            # print(str(omega[n][k]) + 'n:' + str(n) + ' k:' + str(k))

        b[n][1] = omega[n][1]
        for k in range(2, nx-2):
            b[n][k] = b[n][k-1]/a[k-1] + omega[n][k]

        # Condiciones de bordes impuestas
        for k in range(nx-2, 1, -1):
            phi[n+1][k] = 1.0/a[k]*(phi[n+1][k+1]-b[n][k])
        phi[n+1][0] = 0
        phi[n+1][nx-1] = 0


def potencial(x):
    V = 0
    if x >= 300 and x <= 350:
        V = 2.0
    return V


def condicion_inicial(x, parametros):
    q = parametros[0]
    x0 = parametros[1]
    sigma = parametros[2]
    phi0 = np.exp(1j*q*x)*np.exp(-(x-x0)**2/(2*sigma**2))
    return phi0


def normalizar(vector):
    phi_abs = np.absolute(vector)
    phi_sum = 0
    for k in phi_abs:
        phi_sum = k**2 + phi_sum
    vector = vector/phi_sum
    return vector

if __name__ == ("__main__"):
    main()
