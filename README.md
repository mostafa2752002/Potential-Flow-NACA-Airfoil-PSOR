# Potential Flow Analysis over NACA 23021 Airfoil using PSOR

**Author:** Mustafa Taha  
**Course:** AER 4110 - Computational Aerodynamics, Cairo University  
**Date:** November 2025

## 📝 Abstract
This project implements a Computational Fluid Dynamics (CFD) solution for steady, incompressible, inviscid flow over a **NACA 23021** airfoil at an Angle of Attack (AoA) of **7°**.

The solution solves the Laplace equation for the Stream Function ($\Psi$) using the **Finite Difference Method (FDM)** on a body-fitted H-grid topology. The linear system is solved iteratively using the **Point Successive Over-Relaxation (PSOR)** method.

## ⚙️ Key Features
* **Grid Generation:** Custom algebraic H-Grid generation algorithm adaptable for NACA 4 and 5-digit series.
* **Numerical Solver:** PSOR iterative solver for the discretized Laplace equation ($\nabla^2 \Psi = 0$).
* **Post-Processing:** Automated calculation of:
    * Velocity Vector Fields and Magnitude Contours.
    * Pressure Coefficient ($C_p$) distribution.
    * Aerodynamic Coefficients ($C_l, C_d, C_m$).
* **Validation:** Results are validated against XFLR5 data.

## 🧮 Mathematical Model

### Governing Equation
The flow is governed by the Laplace equation for the Stream Function ($\Psi$):
$$\Psi_{xx} + \Psi_{yy} = 0$$


### Coordinate Transformation
To handle the curved airfoil geometry, the physical domain $(x, y)$ is transformed into a computational rectangular domain $(\xi, \eta)$.The governing equation is discretized using Central Difference approximations.

### Numerical Method (PSOR)
The discretized equation is solved using the Point Successive Over-Relaxation method. The iterative update formula is:
$$\Psi_{i,j}^{n+1} = (1-\omega)\Psi_{i,j}^{n} + \frac{\omega}{2(1+\beta^2)} [\Psi_{i+1,j}^{n} + \Psi_{i-1,j}^{n+1} + \beta(\Psi_{i,j+1}^{n} + \Psi_{i,j-1}^{n+1})]$$
Where $\omega$ is the relaxation factor.

## 💻 Usage

1.  Clone the repository.
2.  Open `main.m` in MATLAB.
3.  Adjust parameters in the `%% Inputs` section if necessary:
    ```matlab
    NACA = [2 3 0 2 1];   % Airfoil selection
    Alpha = 7;            % Angle of Attack
    i_max = 200;          % Grid size X
    j_max = 200;          % Grid size Y
    ```
4.  Run the script. The code will generate the grid, solve for $\Psi$, and output plots to the active window.

## 📊 Results

### 1. Grid Generation
An H-grid topology was generated to discretize the domain, ensuring fine resolution near the airfoil surface and leading/trailing edges.

**

### 2. Streamlines & Velocity
The streamlines illustrate the flow physically turning around the airfoil at $\alpha=7^{\circ}$.
**

### 3. Aerodynamic Coefficients
The Lift ($C_l$), Drag ($C_d$), and Pitching Moment ($C_m$) coefficients were calculated by integrating the Pressure Coefficient over the surface.

| Parameter | Current Code (PSOR) | XFLR5 (Validation) |
| :--- | :--- | :--- |
| **Lift ($C_l$)** | 1.2698 | 0.9734 |
| **Drag ($C_d$)** | 0.1559* | 0.0060 |
| **Moment ($C_{m_{0.25c}}$)** | -0.0469 | -0.0162 |



> **Note on Drag:** The discrepancy in $C_d$ is expected. The potential flow assumption neglects viscosity ($Re \to \infty$). In theoretical potential flow, $C_d$ should be zero (d'Alembert's paradox). The calculated value represents numerical error/induced drag from the grid, whereas XFLR5 includes viscous corrections.

## 📉 Convergence
The solution was iterated until the maximum Root Mean Square (RMS) error fell below the tolerance of $10^{-4}$.


## 📚 References
1.  Anderson, J. D., *Computational Fluid Dynamics: The Basics with Applications*.
2.  Cairo University, Faculty of Engineering, *AER 4110 Project Problem Statement*, Nov 2025.
3.  Mustafa Taha, *Solution of Flow Over NACA 23021 Airfoil Using PSOR*, Project Report, 2025.

---
*This project is part of the Aerospace Engineering curriculum at Cairo University.*
