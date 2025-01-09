# Linear Programming Midterm Project 
---
Project: Implement Simplex Methods

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
    - If all b_i >= 0 (basic variables) and all c_N <= 0 (non-basic variables):
      -> The initial dictionary is optimal. Terminate with the optimal solution.

    - If all b_i >= 0 (basic variables) and some c_N > 0 (non-basic variables):
      -> The initial dictionary is primal feasible. Use the Primal Simplex Method.

    - If some b_i < 0 (basic variables) and all c_N <= 0 (non-basic variables):
      -> The initial dictionary is dual feasible. Use the Dual Simplex Method.

    - If both the b vector (basic variables) and c vector (non-basic variables) 
      have negative components (neither primal nor dual feasible):
      -> Use the Two-Phase Simplex Method to handle infeasibility.

 3. Handle degeneracy:
    - If degeneracy is detected, apply Bland's Rule to prevent cycling.

## 3. Contributions
| Teammates | Department | Contribution |
|-----------|------------|--------------|
| **周佳萱 (ME)** | **Institute of Statistics** |  **Simplex method selection Logic, two phase method**  |
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

We add some functions from project1
 - Calculates the shadow price
 - Do the sensitivity analysis
 -Solves the LP again by using path following method
 - Compares the running time: simplex method & path following method 


## 2. Selection Logic
**Part 1**
 - Simplex Method Selection Logic

**Part 2**
 Printing out slack variables, dual variables and dual surplus variables:
    - Slack variables: 
      Calculated during the Primal Simplex method, Bland's Rule and Two-Phase Methods. 
      The values for slack variables are displayed after the final solution is computed. 
      It is calculated based on the difference between the LHS and RHS for each constraint.
    - Dual surplus: 
      Calculated during the Dual Simplex method. 
      Dual surplus measures how much the right-hand side of a constraint can increase without violating dual feasibility. 
    - Dual variables/Shadow prices: 
      The dual variables (or shadow prices) are calculated from the objective coefficients of the basic variables. 
      These represent the marginal value of increasing each constraint's RHS by one unit. 

**Part 3**
 Sensitive analysis:
    - Examines the robustness of the solution by calculating allowable ranges for:
        1. the right-hand side (RHS) constants
        2. objective function coefficients
        3. coefficients in the left-hand side (LHS) matrix 
    - Determines the extent to which these values can change without altering the optimal basis of the solution. 
    - For each constraint and variable, it outputs the allowable range, enabling an understanding of the solution's sensitivity to data changes in the model.

**Part 4**
 Path-following method: 
    - Initializing primal (x, ω) and dual (y, z) variables with positive feasible values:
      These variables must satisfy the linear constraints defined by matrix A and vectors b, c 
      where:
        - x and ω are primal variables for the linear program.
        - y and z are dual variables corresponding to the constraints and objective function.

    - Calculating residuals to measure the feasibility of the current solution:
        - Primal residual, ρ = b - A  x - ω, checks if the primal variables satisfy the constraints.
        - Dual residual, σ = c - A^T  y + z, verifies the constraints on the dual side.
      These residuals indicate how far the current solution is from being feasible in both the primal and dual spaces.
  
    - Computing the duality measure (μ) to assess convergence:
      The duality measure is calculated as μ = δ  γ / (n + m), 
      where:
        - δ (delta) controls the convergence rate.
        - γ = x^T  z + ω^T  y is the duality gap, reflecting the difference between primal 
          and dual objective values.
      This measure μ determines the progress toward the optimal solution; the method stops once μ is below a predefined tolerance level.
  
    - Constructing and solving a linear system to find adjustments (deltas) for primal and dual variables:
      The system of equations includes:
        - A * Δx + Δω = ρ, ensuring primal residual reduction.
        - A^T * Δy - Δz = σ, ensuring dual residual reduction.
        - Z * Δx + X * Δz = μ - x * z, enforcing centrality by maintaining a positive product of x and z.
        - W * Δy + Y * Δω = μ - y * ω, maintaining positive products for y and ω.
     Solving this system provides updates (Δx, Δy, Δz, Δω) to improve feasibility and reduce the duality gap.
 
    - Updating variables and recalculating the duality measure iteratively:
      Using the step length θ, calculated to keep updates within the feasible region, variables are adjusted as:
        - x = x + θ * Δx, y = y + θ * Δy, z = z + θ * Δz, and ω = ω + θ * Δω.
        - After each update, the duality measure μ is recalculated to check convergence.
      The process repeats until μ < tolerance or the maximum iteration limit is reached.


## 3. Contributions
| Teammates | Department | Contribution |
|-----------|------------|--------------|
| **周佳萱 (ME)** | **Institute of Statistics** | **Path-following method** |
|楊雅茗| Department of Industrial Engineering and Management| Sensitive analysis |
|吳和晏| Department of Industrial Engineering and Management | Path-following method | 
|王睿言| Department of Industrial Engineering and Management | Printing out slack variables, dual variables and dual surplus variables |







