# Linear Programming Midterm Project 
---
Project: Implement Simplex Methods

Author: Group 4

Date: 2024-10-07

Course: Linear Programming

---
## 1. Objective

- Reads and parses LP data from a CSV file.
- Determines the appropriate simplex method (Primal, Dual, Bland's Rule, Two-Phase).
- Implements Bland's Rule to handle degeneracy and Two-Phase Method for infeasibility.
- Tracks iteration count manually (without built-in methods).
- Outputs all non-zero solution variables (including slack variables).


## 2. Simplex Method Selection Logic:

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

