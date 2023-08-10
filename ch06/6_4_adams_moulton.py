#!/usr/bin/env python3
from math import pi, sin


def f(x: float, y: float) -> float:
    return y * sin(pi*x) 


def predictor(x_vals: list[float], y_vals: list[float], h: float) -> float:
    f_vals = 55.*f(x_vals[3], y_vals[3]) - 59.*f(x_vals[2], y_vals[2])
    f_vals += 37.*f(x_vals[1], y_vals[1]) - 9.*f(x_vals[0], y_vals[0])
    return y_vals[3] + (h/24.)*f_vals


def corrector(x_vals: list[float], y_vals: list[float], h: float) -> float:
    f_vals = 9.*f(x_vals[3], y_vals[3]) + 19.*f(x_vals[2], y_vals[2])
    f_vals += -5.*f(x_vals[1], y_vals[1]) + f(x_vals[0], y_vals[0])
    return y_vals[2] + (h/24.)*f_vals


def main():
    output_vals = [1., 1.0626782047, 1.2460111437, 1.5169078366] 
    h = .2
    x_vals = [h*float(i) for i in range(6)]
    
    for i in range(2):
        predicted = predictor(x_vals[i:i+5], output_vals[i:i+5], h)
        print(f"y_{i+4},p = {predicted}")
        output_vals.append(predicted)
        corrected = corrector(x_vals[i+1:i+6], output_vals[i+1:i+6], h)
        output_vals[i+4] = corrected
        print(f"y_{i+4},c = {corrected}")


if __name__ == "__main__":
    main()
