#!/usr/bin/env python3
import numpy as np

from copy import deepcopy


def print_matrix(u: np.array) -> None:
    print("")
    for i in range(u.shape[0]):
        print("[", end="")
        for j in range(u.shape[1]):
            if j < u.shape[1] - 1:
                print("{:.9f}".format(u[i,j]), end="\t")
            else:
                print("{:.9f}".format(u[i,j]), end="]\n")
            

def calc_matrix(boundary_matrix: np.array, p: int, q: int) -> np.array:
    for i in range(1, p-1):
        for j in range(1, q-1):
            boundary_matrix[i,j] = (boundary_matrix[0,j] + boundary_matrix[i,0] +  
                    boundary_matrix[p-1,j] + boundary_matrix[i,q-1]
            ) / 4.

    return boundary_matrix


def calc_omega(p: int, q: int) -> float:
    a = (np.cos(np.pi/p) + np.cos(np.pi/q))**2
    b = -16.
    c = 16.

    omega_1 = (-b + np.sqrt(b**2 - 4.*a*c)) / (2.*a)
    omega_2 = (-b - np.sqrt(b**2 - 4.*a*c)) / (2.*a)

    return min(omega_1, omega_2)


def calc_norm(u: np.array, u_hat: np.array) -> float:
    norm = 0.

    for i in range(1, u.shape[0]-1):
        for j in range(1, u.shape[1]-1):
            new_norm = abs(u[i,j] - u_hat[i,j]) 
            if new_norm > norm:
                norm = new_norm

    return norm


def sor(boundary_matrix: np.array, tol: float, max_iter: int) -> np.array:
    print("This is the SOR Method for Laplace's Equation.")

    p = boundary_matrix.shape[0]
    q = boundary_matrix.shape[1]
    u = calc_matrix(boundary_matrix, p, q)
    print("\nU = ", end="")
    print_matrix(u)

    omega = calc_omega(p-1, q-1)
    print("Value of omega = {:.9f}".format(omega))

    for k in range(1, max_iter+1):
        u_hat = deepcopy(u)
        for i in range(1, p-1):
            for j in range(1, q-1):
                u_hat[i,j] = ((1.-omega)*u[i,j] + 
                        (omega/4.)*(u_hat[i,j-1] + u_hat[i,j+1] + 
                                    u_hat[i-1,j] + u_hat[i+1,j])
                )

        norm = calc_norm(u, u_hat)
        if norm < tol:
            print(f"\nThe procedure was successful after {k} iterations.")
            print("||U - U^|| = {:.9f} < TOL = {:.3f}".format(norm, tol))
            return u_hat

        u = deepcopy(u_hat)

    print("\nThe procedure was not successful, U = ", end="")
    print_matrix(u)


def main() -> None:
    boundary_matrix = np.array([
        [100., 90., 70., 50., 0.],
        [120., 0., 0., 0., 0.],
        [150., 0., 0., 0., 0.],
        [150., 0., 0., 0., 0.],
        [120., 0., 0., 0., 0.],
        [100., 70., 50., 30., 0.]
    ])
    tol = 1e-3
    max_iter = 20

    u = sor(boundary_matrix, tol, max_iter)

    if u is not None:
        print("\nResult: U = ", end="")
        print_matrix(u)


if __name__ == "__main__":
    main()
