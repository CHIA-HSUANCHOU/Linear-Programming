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


---

Project 3: Farmer’s Problem with Bendor Decomposition

Author: Group 4

Date: 2024-12-16

Course: Linear Programming

---

## 1. Overviews

This project solves a farmer's crop allocation problem using **Bender Decomposition** to optimize planting decisions while accounting for uncertainties in yield.

## 2. Problem Description

- **Crops**: Wheat, Corn, Sugar Beets  
- **Available Land**: 500 acres  

| Crop         | Yield (T/acre) | Planting Cost ($/acre) | Selling Price ($/T)                 | Purchase Price ($/T) | Minimum Demand (T) |
|--------------|----------------|------------------------|-------------------------------------|-----------------------|--------------------|
| **Wheat**    | 2.5            | 150                    | 170                                 | 238                   | 200                |
| **Corn**     | 3.0            | 230                    | 150                                 | 210                   | 240                |
| **Sugar Beets** | 20.0         | 260                    | 36 (up to 6000 T), 10 (above 6000 T) | Not Applicable    | Not Applicable    |

## Objective Function
(Take number of scenario = 3 and choose C1 as first-stage constraints)

The objective is to **minimize total costs**, which include:
1.**First-Stage Planting Cost**
150x_1 + 230x_2 + 260x_3

2. **Second-Stage Shortfall and Surplus Costs**:
238y_1 + 210y_2 - 170w_1 - 150w_2 - 36w_3 - 10w_4

## Constraints

- **Total Land Constraint**:x_1 + x_2 + x_3 <= 500

- **Wheat Yield Shortfall**:2.5x_1 + y_1 - w_1 >= 200 (rearranged as: -2.5x_1 - y_1 + w_1 <= -200)

- **Corn Yield Shortfall**:3x_2 + y_2 - w_2 >= 240 (rearranged as: -3x_2 - y_2 + w_2 <= -240)

- **Sugar Beet Selling Capacity**:w_3 + w_4 <= 20x_3 (rearranged as: w_3 + w_4 - 20x_3 <= 0)

- **Sugar Beet Tier Limit**: w_3 <= 6000

## 3. Flow:
1.User input:
Choose scenario and Constraints in Master


→ Build master problem → Solve master problem → Solve subproblems → Add cuts → Iterate until convergence.





## 5. Contributions
| Teammates | Department | Contribution |
|-----------|------------|--------------|
| **周佳萱 (ME)** | **Institute of Statistics** |  **Code-Structure, L-shape method, Optimility Cut, Feasibility Cut**  |
|楊雅茗| Department of Industrial Engineering and Management| Code Stucture,Toolbox Implementation, Optimility Cut, ppt |
|吳和晏| Department of Industrial Engineering and Management | Code Structure,Toolbox Implementation, Optimility Cut, ppt | 
|王睿言| Department of Industrial Engineering and Management | L-shape method, Optimility Cut, Feasibility Cut, ppt|






