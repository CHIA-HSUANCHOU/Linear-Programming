# Linear Programming Midterm Project 

## **Project 1: Implement Simplex Methods**
---
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
   - **Case 1:** If all `b[i] >= 0` (basic variables) and all `c_N[j] <= 0` (non-basic variables): 
     → The initial dictionary is optimal. Terminate with the optimal solution.

   - **Case 2:** If all `b[i] >= 0` (basic variables) and some `c_N[j] > 0` (non-basic variables): 
     → The initial dictionary is primal feasible. Use the **Primal Simplex Method**.

   - **Case 3:** If some `b[i] < 0` (basic variables) and all `c_N[j] <= 0` (non-basic variables): 
     → The initial dictionary is dual feasible. Use the **Dual Simplex Method**.

   - **Case 4:**  If both the `b` vector (basic variables) and `c` vector (non-basic variables)  
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





##**Project 2: Implement Simplex Methods, sensitivity analysis, path-following method**
---
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



##**Project 3: Farmer’s Problem with Bender Decomposition**

---

Author: Group 4

Date: 2024-12-16

Course: Linear Programming

---

## 1. Overviews

This project solves a farmer's problem using **Bender Decomposition** to optimize planting decisions while accounting for uncertainties in yield.

## 2. Problem Description

- **Crops**: Wheat, Corn, Sugar Beets  
- **Available Land**: 500 acres  

| Crop         | Yield (T/acre) | Planting Cost ($/acre) | Selling Price ($/T)                 | Purchase Price ($/T) | Minimum Demand (T) |
|--------------|----------------|------------------------|-------------------------------------|-----------------------|--------------------|
| **Wheat**    | 2.5            | 150                    | 170                                 | 238                   | 200                |
| **Corn**     | 3.0            | 230                    | 150                                 | 210                   | 240                |
| **Sugar Beets** | 20.0         | 260                    | 36 (up to 6000 T), 10 (above 6000 T) | Not Applicable    | Not Applicable    |

## 3. Objective Function

(Taking the number of scenarios = 3 and choosing C1 as first-stage constraints)

### First-Stage Objective Function
```
150x_1 + 230x_2 + 260x_3
```
### Second-Stage Objective Function
```
1/3*(238y_11 + 210y_21 - 170w_11 - 150w_21 - 36w_31 - 10w_41)

+1/3*(238y_12 + 210y_22 - 170w_12 - 150w_22 - 36w_32 - 10w_42)

+1/3*(238y_13 + 210y_23 - 170w_13 - 150w_23 - 36w_33 - 10w_43)
```
## 4. Constraints

### First-Stage Constraints
```
C1:x_1 + x_2 + x_3 <= 500
```
### Second-Stage Constraints
```
C2_S1: -2.0*x_1 + -1.0*y_1_1 + 1.0*w_1_1 ≤ -200.0

C3_S1: -2.4*x_2 + -1.0*y_2_1 + 1.0*w_2_1 ≤ -240.0

C4_S1: -16.0*x_3 + 1.0*w_3_1 + 1.0*w_4_1 ≤ 0.0

C5_S1: 1.0*w_3_1 ≤ 6000.0

C2_S2: -2.5*x_1 + -1.0*y_1_2 + 1.0*w_1_2 ≤ -200.0

C3_S2: -3.0*x_2 + -1.0*y_2_2 + 1.0*w_2_2 ≤ -240.0

C4_S2: -20.0*x_3 + 1.0*w_3_2 + 1.0*w_4_2 ≤ 0.0

C5_S2: 1.0*w_3_2 ≤ 6000.0

C2_S3: -3.0*x_1 + -1.0*y_1_3 + 1.0*w_1_3 ≤ -200.0

C3_S3: -3.6*x_2 + -1.0*y_2_3 + 1.0*w_2_3 ≤ -240.0

C4_S3: -24.0*x_3 + 1.0*w_3_3 + 1.0*w_4_3 ≤ 0.0

C5_S3: 1.0*w_3_3 ≤ 6000.0
```
## 5. Flow:

1. **User Input**:
   - Select the number of scenarios and specify constraints for the master problem.

2. **Solve Master Problem**:
   - Solve the master problem using the primal simplex method. If infeasible, switch to the dual simplex method.
   - Minimize:
     ```
     150x_1 + 230x_2 + 260x_3
     ```
   - Constraint:
     ```
     C1: x_1 + x_2 + x_3 <= 500
     ```
   - If the master problem is a complete recourse:
     - No need to add feasibility cuts. Only optimal cuts are required.
 
3. **Solve subproblems**
- Check feasibility /Add Feasibility Cut

  -If the problem is complete recourse => No need to do feasibility test
  
- Solve new master problem and resolve subproblem
  
  -If feasibility cut added, solve the master problem and subproblem again
  
- Check  optimality/Add  optimality cut
  
  -After all scenarios of the subproblem  have been run through and optimality  has not been reached, add optimality cut to the master problem 

4. Iterate until convergence
   - Stop when the **Upper Bound (UB)** and **Lower Bound (LB)** converge:
     ```
     UB = min(UB, objval)
     LB = max(LB, thetaVal)
     ```
   - Or the maximum number of iterations (500) is reached.

## 6. Contributions
| Teammates | Department | Contribution |
|-----------|------------|--------------|
| **周佳萱 (ME)** | **Institute of Statistics** |  **Code-Structure, L-shape method, Optimility Cut, Feasibility Cut**  |
|楊雅茗| Department of Industrial Engineering and Management| Code Stucture,Toolbox Implementation, Optimility Cut|
|吳和晏| Department of Industrial Engineering and Management | Code Structure,Toolbox Implementation, Optimility Cut| 
|王睿言| Department of Industrial Engineering and Management | L-shape method, Optimility Cut, Feasibility Cut|






