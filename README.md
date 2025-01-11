# Linear Programming Midterm Project 
---
Project 1: Implement Simplex Methods

Author: Group 4

Date: 2024-10-07

Course: Linear Programming

---
## 1. Main Features

- Reads and parses LP data from a CSV file.
- Determines the appropriate simplex method (Primal, Dual, Bland's Rule, Two-Phase).
- Implements Bland's Rule to handle degeneracy and Two-Phase Method for infeasibility.
- Tracks iteration count manually (without built-in methods).
- Outputs all non-zero solution variables (including slack variables).


## 2. Simplex Method Selection Logic

 1. Check for unboundedness:
    - If the LP is unbounded, terminate and return an unbounded solution.

 2. Check primal/dual feasibility to determine the simplex method:
   - **Case 1:** If all \( b_i \geq 0 \) (basic variables) and all \( c_N \leq 0 \) (non-basic variables):  
     → The initial dictionary is optimal. Terminate with the optimal solution.

   - **Case 2:** If all \( b_i \geq 0 \) (basic variables) and some \( c_N > 0 \) (non-basic variables):  
     → The initial dictionary is primal feasible. Use the **Primal Simplex Method**.

   - **Case 3:** If some \( b_i < 0 \) (basic variables) and all \( c_N \leq 0 \) (non-basic variables):  
     → The initial dictionary is dual feasible. Use the **Dual Simplex Method**.

   - **Case 4:** If both the \( b \) vector (basic variables) and \( c \) vector (non-basic variables)  
     have negative components (neither primal nor dual feasible):  
     → Use the **Two-Phase Simplex Method** to handle infeasibility.


 3. Handle degeneracy:
    - If degeneracy is detected, apply **Bland's Rule** to prevent cycling.

## 3. Contributions
| Teammates | Department | Contribution |
|-----------|------------|--------------|
| **周佳萱 (ME)** | **Institute of Statistics** |  **Simplex Method Selection Logic,  Two phase method**  |
|楊雅茗| Department of Industrial Engineering and Management| Simplex Method Selection Logic, Two phase method |
|吳和晏| Department of Industrial Engineering and Management | Simplex Method Selection Logic, Bland's Rule | 
|王睿言| Department of Industrial Engineering and Management | Simplex Method Selection Logic, Bland's Rule |



---
Project 2: Implement Simplex Methods, sensitivity analysis, path-following method

Author: Group 4

Date: 2024-11-11

Course: Linear Programming

---

## 1. Main Features

We have added the following features from Project 1:

- Calculate the shadow price.
- Perform sensitivity analysis.
- Solve the LP again using the path-following method.
- Compare the running time of the Simplex Method and the Path-Following Method.


## 2. Selection Logic

### **Part 1: Simplex Method Selection Logic**

Logic for choosing between the various simplex methods based on the problem's feasibility and degeneracy.(Same as project 1)

### **Part 2: Printing Slack Variables, Dual Variables, and Dual Surplus Variables**

- **Slack Variables**:
    - Calculated during the Primal Simplex Method, Bland's Rule, and Two-Phase Methods.
    - Values for slack variables are displayed after computing the final solution.
    - Defined as the difference between the LHS and RHS of each constraint.

- **Dual Surplus**:
    - Calculated during the Dual Simplex Method.
    - Represents how much the right-hand side of a constraint can increase without violating dual feasibility.

- **Dual Variables/Shadow Prices**:
    - Derived from the objective coefficients of the basic variables.
    - Represent the marginal value of increasing the RHS of a constraint by one unit.


### **Part 3: Sensitivity Analysis**

- Examines the robustness of the solution by calculating allowable ranges for:
    1. Right-hand side (RHS) constants.
    2. Objective function coefficients.
    3. Coefficients in the left-hand side (LHS) matrix.

- Determines the extent to which these values can change without altering the optimal basis of the solution.
- Outputs the allowable ranges for each constraint and variable, providing insights into the solution's sensitivity to data changes in the model.


### Part 4: Path-Following Method

1. **Initialization**:
    - Primal variables (`x`, `ω`) and dual variables (`y`, `z`) are initialized with positive feasible values.
    - Variables must satisfy linear constraints defined by matrix `A` and vectors `b`, `c`:
        - `x`, `ω`: Primal variables.
        - `y`, `z`: Dual variables corresponding to constraints and the objective function.

2. **Residual Calculation**:
    - **Primal Residual (`ρ`)**:
      ```
      ρ = b - A * x - ω
      ```
      Ensures primal variables satisfy constraints.
    - **Dual Residual (`σ`)**:
      ```
      σ = c - A^T * y + z
      ```
      Ensures dual variables satisfy constraints.

3. **Duality Measure (`μ`)**:
    - Measure of convergence, calculated as:
      ```
      μ = (δ * γ) / (n + m)
      ```
      Where:
        - `δ`: Convergence rate.
        - `γ = x^T * z + ω^T * y`: Duality gap, reflecting the difference between primal and dual objective values.

4. **Solve Linear System**:
    - A system of equations is solved to find updates (`Δx`, `Δy`, `Δz`, `Δω`):
        - Primal residual reduction:
          ```
          A * Δx + Δω = ρ
          ```
        - Dual residual reduction:
          ```
          A^T * Δy - Δz = σ
          ```
        - Centrality maintenance:
          ```
          Z * Δx + X * Δz = μ - x * z
          W * Δy + Y * Δω = μ - y * ω
          ```

5. **Variable Updates**:
    - Adjust variables iteratively using step length (`θ`):
      ```
      x = x + θ * Δx
      y = y + θ * Δy
      z = z + θ * Δz
      ω = ω + θ * Δω
      ```
    - Recalculate `μ` and repeat until `μ < tolerance` or the maximum iteration limit is reached.

## 3. Comparison

- Compare the running times of the Simplex Method and Path-Following Method.

## 4. Contributions
| Teammates | Department | Contribution |
|-----------|------------|--------------|
| **周佳萱 (ME)** | **Institute of Statistics** | **Path-following method** |
|楊雅茗| Department of Industrial Engineering and Management| Sensitive analysis |
|吳和晏| Department of Industrial Engineering and Management | Path-following method | 
|王睿言| Department of Industrial Engineering and Management | Printing out slack variables, dual variables and dual surplus variables |







