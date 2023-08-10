#!/usr/bin/env python3
import numpy as np
from scipy import sparse

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


def calc_norm(u_hat: np.array, u: np.array) -> float:
    norm = 0.

    for i in range(1, u_hat.shape[0]-1):
        for j in range(1, u_hat.shape[1]-1):
            new_norm = abs(u_hat[i,j] - u[i,j]) 
            if new_norm > norm:
                norm = new_norm

    return norm


def get_horizontal_coeffs(a: sparse.csr_matrix, rho: float) -> sparse.csr_matrix:
    r = (1./rho) + 2.
    u = np.array([
            [r, -1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [-1., r, -1., 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, -1., r, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, r, -1., 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, -1., r, -1., 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -1., r, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, r, -1., 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, -1., r, -1., 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, -1., r, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, r, -1., 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, -1., r, -1.],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1., r]
    ])

    return sparse.csr_matrix(u)


def get_vertical_coeffs(a: sparse.csr_matrix, rho: float) -> sparse.csr_matrix:
    r = (1./rho) + 2.
    u = np.array([
            [r, -1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [-1., r, -1., 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, -1., r, -1., 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, -1., r, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, r, -1., 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -1., r, -1., 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, -1., r, -1., 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, -1., r, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, r, -1., 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, -1., r, -1., 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, -1., r, -1.],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1., r]
    ])

    return sparse.csr_matrix(u)


def get_vertical_values(a: np.array, p: int, q: int) -> np.array:
    v = np.zeros(12)
    k = -2
    for j in range(1, q-1):
        for i in range(1, p-1):
            v[i+j+k] = a[i,j]
        k += 3

    return v


def set_vertical_values(a: np.array, v: np.array, p: int, q: int) -> np.array:
    k = -2
    for j in range(1, q-1):
        for i in range(1, p-1):
            a[i,j] = v[i+j+k]
        k += 3

    return a


def horizontal_traverse(a: np.array, u: sparse.csr_matrix, 
        v: np.array, rho: float, p: int, q: int
) -> np.array:
    r = (1./rho) - 2.
    b = np.array([[
        a[1,0] + a[0,1] + r*v[0] + v[1],
        a[0,2] + r*v[4] + v[5],
        a[1,q-1] + a[0,3] + r*v[8] + v[9],
        a[2,0] + v[0] + r*v[2] + v[3],
        v[4] + r*v[5] + v[6],
        a[2,q-1] + v[8] + r*v[9] + v[10],
        a[3,0] + v[1] + r*v[2] + v[3],
        v[5] + r*v[6] + v[7],
        a[3,q-1] + v[9] + r*v[10] + v[11],
        a[4,0] + v[2] + r*v[3] + a[p-1,1],
        v[6] + r*v[7] + a[p-1,2],
        a[4,q-1] + v[10] + r*v[11] + a[p-1,3]
    ]])

    return solve(u, b.transpose())


def vertical_traverse(a: np.array, v: sparse.csr_matrix, 
        u: np.array, rho: float, p: int, q: int
) -> np.array:
    r = (1./rho) - 2.
    b = np.array([[
        a[0,1] + a[1,0] + r*u[0] + u[1],
        a[2,0] + r*u[3] + u[4],
        a[3,0] + r*u[6] + u[7],
        a[p-1,1] + a[4,0] + r*u[9] + u[10],
        a[0,2] + u[0] + r*u[1] + u[2],
        u[3] + r*u[4] + u[5],
        u[6] + r*u[7] + u[8],
        a[p-1,2] + u[9] + r*u[10] + u[11],
        a[0,3] + u[1] + r*u[2] + a[1,q-1],
        u[4] + r*u[5] + a[2,q-1],
        u[7] + r*u[8] + a[3,q-1],
        a[p-1,3] + u[10] + r*u[11] + a[4,q-1]
    ]])

    return solve(v, b.transpose())


def solve(u: sparse.csr_matrix, b: np.array) -> np.array:
    return sparse.linalg.spsolve(u, b)


def adi(boundary_matrix: np.array, 
        rho: float, tol: float, max_iter: int
) -> np.array:
    print("This is the A.D.I. Method for Laplace's Equation.")

    p = boundary_matrix.shape[0]
    q = boundary_matrix.shape[1]
    a = calc_matrix(boundary_matrix, p, q)
    print("\nA = ", end="")
    print_matrix(a)

    print("Value of rho = {:.1f}".format(rho))

    u_coeffs = get_horizontal_coeffs(a, rho)
    v_coeffs = get_vertical_coeffs(a, rho)

    for k in range(1, max_iter+1, 2):
        a_hat = deepcopy(a)

        v = get_vertical_values(a_hat, p, q)
        u = horizontal_traverse(a, u_coeffs, v, rho, p, q)
        v = vertical_traverse(a, v_coeffs, u, rho, p, q)
        a_hat = set_vertical_values(a_hat, v, p, q)

        norm = calc_norm(a_hat, a)
        if norm < tol:
            print(f"\nThe procedure was successful after {k} iterations.")
            print("||A^(k+2) - A^(k)|| = {:.9f} < TOL = {:.5f}".format(norm, tol))
            return a_hat

        a = deepcopy(a_hat)

    print("\nThe procedure was not successful, A = ", end="")
    print_matrix(a)


def main() -> None:
    boundary_matrix = np.array([
        [100., 90., 70., 50., 0.],
        [120., 0., 0., 0., 0.],
        [150., 0., 0., 0., 0.],
        [150., 0., 0., 0., 0.],
        [120., 0., 0., 0., 0.],
        [100., 70., 50., 30., 0.]
    ])
    rho = 1.
    tol = 1e-5
    max_iter = 30

    u = adi(boundary_matrix, rho, tol, max_iter)

    if u is not None:
        print("\nResult: U = ", end="")
        print_matrix(u)


if __name__ == "__main__":
    main()
