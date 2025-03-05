# Thermodynamic Properties of Real Fluids
# Liquid-Vapor Equilibrium (Saturation Curve)

# Importing libraries
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Provided table data
T = [320, 350, 370, 420, 470, 520, 570, 620, 650, 670, 690]
Pexp = [0.000283, 0.00173, 0.00473, 0.0348, 0.1501, 0.4552, 1.095, 2.2826, 3.4118, 4.4262, 5.7334]

# Constants
# Triple point
PPT = 0.000188
TPT = 314.06
# Critical point
PC = 6.13
TC = 694.25
# Acentric factor
w = 0.44346
# Gas constant
R = 8.3145

# Functions
# Function to calculate TR, m, and alpha
def calculate_factors(T):
    TR = T / TC    
    m = 0.37464 + (1.54226*w) - (0.26992*w*w)    
    alpha = (1 + m * (1 - math.sqrt(TR)))**2    
    return TR, m, alpha

# Function to calculate constants a and b
def calculate_constants(R, TC, alpha, PC):
    a = 0.45724 * (R**2 * TC**2 * alpha) / PC
    b = 0.0778 * R * TC / PC
    return a, b

# Function to calculate dimensionless constants A and B
def calculate_dimensionless_constants(a, b, R, T, Pcalc):
    A = (a * Pcalc) / ( R**2 * T**2)
    B = (b * Pcalc) / (R * T)
    return A, B

# Function to calculate Zl and Zv using the polynomial roots
def calculate_roots(A, B):
    x1 = 1
    x2 = -(1 - B)
    x3 = A - 2*B - 3 * B**2
    x4 = - A * B + B**2 + B**3
    
    coefficients = [x1, x2, x3, x4]
    roots = np.roots(coefficients)

    Zl = min(roots)
    Zv = max(roots)
    return Zl, Zv

# Function to calculate the fugacity coefficient of z
def fugacity_coefficient(z, A, B):
    term1 = (z - 1)
    term2 = math.log(z - B)
    term3_numerator = z + (1 - math.sqrt(2)) * B
    term3_denominator = z + (1 + math.sqrt(2)) * B
    term3 = (A / (2 * math.sqrt(2) * B)) * math.log(term3_numerator / term3_denominator)
    
    ln_phi = term1 - term2 + term3   
    phi = math.exp(ln_phi)  
    return phi

# Main function to find the calculated pressure
def equilibrium_pressure(T, Pexp):
    Pcalc = Pexp
    TR, m, alpha = calculate_factors(T)

    a, b = calculate_constants(R, TC, alpha, PC)
    A, B = calculate_dimensionless_constants(a, b, R, T, Pcalc)

    Zl, Zv = calculate_roots(A, B)
    phi_l = fugacity_coefficient(Zl, A, B)
    phi_v = fugacity_coefficient(Zv, A, B)

    if (((phi_v / phi_l) - 1) < 0.0001):
        Pcalc = Pexp
    else:
        Pcalc = Pcalc * (phi_l / phi_v)
        equilibrium_pressure(T, Pcalc)

    return Pcalc

# Calculating pressure for each instance in the table
Pcalc_vetor = []
for i in range(len(T)):
    value = equilibrium_pressure(T[i], Pexp[i])
    Pcalc_vetor.append(value)

# Calculated data
data = {'Temperature (K)': T, 'Experimental Pressure (MPa)': Pexp, 'Calculated Pressure (MPa)': Pcalc_vetor}
df = pd.DataFrame(data)
print(df.to_string(index=False))

# Graphs
# Experimental pressure graph
plt.figure(1)
plt.plot(T, Pexp, marker='o', label='Experimental Pressure', color='blue')
plt.title('Temperature vs Experimental Pressure')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (MPa)')
plt.legend()
plt.grid(True)

# Calculated pressure graph
plt.figure(2)
plt.plot(T, Pexp, marker='o', label='Calculated Pressure', color='red')
plt.title('Temperature vs Calculated Pressure')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (MPa)')
plt.legend()
plt.grid(True)

plt.show()
