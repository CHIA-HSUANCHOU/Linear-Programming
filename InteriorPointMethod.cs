/*
Midterm2 project
組員：王睿言、吳和晏、楊雅茗、周佳萱

 Main Features:
 - Reads and parses LP data from a CSV file.
 - Determines the appropriate simplex method (Primal, Dual, Bland's Rule, Two-Phase).
 - Implements Bland's Rule to handle degeneracy and Two-Phase Method for infeasibility.
 - Tracks iteration count manually (without built-in methods).
 - Outputs all non-zero solution variables (including slack variables).
 - Calculates the shadow price
 - Do the sensitivity analysis
 - Solves the LP again by using path following method
 - Compares the running time: simplex method & path following method 

Part 1
 Simplex Method Selection Logic:

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

Part 2
 Printing out slack variables, dual varaibles and dual surplus variables:
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

Part 3
 Sensitive analysis:
    - Examines the robustness of the solution by calculating allowable ranges for:
        1. the right-hand side (RHS) constants
        2. objective function coefficients
        3. coefficients in the left-hand side (LHS) matrix 
    - Determines the extent to which these values can change without altering the optimal basis of the solution. 
    - For each constraint and variable, it outputs the allowable range, enabling an understanding of the solution's sensitivity to data changes in the model.

Part 4
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

    - Test data:
        1. data.csv
            -> initial value：x = 200, z = 800, omega = 650, y = 400 
            -> d = 0.15, gamma = 0.85, theta = 1, r = 0.965

        2. 20-Klee-Minty Instance.csv
            -> Since the range of the lhs and rhs variables is very large, from 1 to around 10^28, the primal solution for x also spans a large range.
            -> It challenging to initialize values and ensure convergence within a finite number of iterations in the path-following method.
*/

using System;
using System.IO;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using static System.Runtime.InteropServices.JavaScript.JSType;

class ResourceAllocationProblem
{
    static double[] objCoeffs;

    // Solve Using Primal Simplex Method
    static void SolveUsingPrimalSimplex(double[,] lhs, double[] rhs,
                                    double[] objectiveCoefficients, int totalVars, int numConstraints)
    {
        Stopwatch stopwatch = new Stopwatch();
        stopwatch.Start();

        bool shouldContinue = true;
        int iterationCount = 0;
        int originalVars = totalVars - numConstraints;

        // Store the original objective coefficients
        double[] objCoeffs = new double[totalVars];
        Array.Copy(objectiveCoefficients, objCoeffs, totalVars);

        // Initialize basic variables (slack variables)
        int[] basicVariables = new int[numConstraints];
        for (int i = 0; i < numConstraints; i++)
        {
            basicVariables[i] = originalVars + i;
        }

        while (shouldContinue)
        {
            iterationCount++;
            Console.WriteLine($"\n--- Iteration {iterationCount} ---");

            // Calculate the dual variables (pi)
            double[] pi = new double[numConstraints];
            for (int i = 0; i < numConstraints; i++)
            {
                pi[i] = objectiveCoefficients[basicVariables[i]];
            }

            // Compute reduced costs for all variables
            double[] reducedCosts = new double[totalVars];
            for (int j = 0; j < totalVars; j++)
            {
                reducedCosts[j] = objectiveCoefficients[j];
                for (int i = 0; i < numConstraints; i++)
                {
                    reducedCosts[j] -= pi[i] * lhs[i, j];
                }
            }

            // Identify the entering variable (most positive reduced cost)
            int enteringVarIndex = -1;
            double maxReducedCost = 0;
            for (int j = 0; j < totalVars; j++)
            {
                if (reducedCosts[j] > maxReducedCost)
                {
                    maxReducedCost = reducedCosts[j];
                    enteringVarIndex = j;
                }
            }

            // If no entering variable is found, the solution is optimal
            if (enteringVarIndex == -1)
            {
                Console.WriteLine("Optimal solution found.");
                shouldContinue = false;
                break;
            }

            Console.WriteLine($"Chosen entering variable: {(enteringVarIndex < originalVars ? "x" : "s")}{enteringVarIndex + 1}");

            // Perform ratio test to select leaving variable
            int leavingVarIndex = -1;
            double minRatio = double.MaxValue;
            for (int i = 0; i < numConstraints; i++)
            {
                if (lhs[i, enteringVarIndex] > 0)
                {
                    double ratio = rhs[i] / lhs[i, enteringVarIndex];
                    if (ratio < minRatio)
                    {
                        minRatio = ratio;
                        leavingVarIndex = i;
                    }
                }
            }

            // If no leaving variable is found, the problem is unbounded
            if (leavingVarIndex == -1)
            {
                Console.WriteLine("Problem is unbounded.");
                shouldContinue = false;
                break;
            }

            Console.WriteLine($"Chosen leaving variable: {(basicVariables[leavingVarIndex] < originalVars ? "x" : "s")}{basicVariables[leavingVarIndex] + 1}");

            // Perform pivot operation
            double pivotElement = lhs[leavingVarIndex, enteringVarIndex];

            // Update the leaving row in-place
            for (int j = 0; j < totalVars; j++)
            {
                lhs[leavingVarIndex, j] /= pivotElement;
            }
            rhs[leavingVarIndex] /= pivotElement;

            // Update all other rows in-place
            for (int i = 0; i < numConstraints; i++)
            {
                if (i != leavingVarIndex)
                {
                    double factor = lhs[i, enteringVarIndex];
                    for (int j = 0; j < totalVars; j++)
                    {
                        lhs[i, j] -= factor * lhs[leavingVarIndex, j];
                    }
                    rhs[i] -= factor * rhs[leavingVarIndex];
                }
            }

            // Update the objective coefficients in-place
            double objectiveFactor = objectiveCoefficients[enteringVarIndex];
            for (int j = 0; j < totalVars; j++)
            {
                objectiveCoefficients[j] -= objectiveFactor * lhs[leavingVarIndex, j];
            }

            // Update the basic variable index in-place
            basicVariables[leavingVarIndex] = enteringVarIndex;

            // Print the current tableau
            //PrintCurrentDictionaryForPrimal(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);
        }

        // Output the final solution with shadow prices and slack variables
        OutputFinalSolutionWithShadowPricesAndSlacks(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);
        PerformSensitivityAnalysis(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);

        stopwatch.Stop();
        Console.WriteLine($"Total runtime for Primal method: {stopwatch.Elapsed.TotalSeconds} seconds");
    }


    // Shared output function to display the final solution, shadow prices, and slack variables
    static void OutputFinalSolutionWithShadowPricesAndSlacks(double[] objCoeffs, double[,] lhs, double[] rhs, int[] basicVariables, int totalVars, int numConstraints)
    {
        int originalVars = totalVars - numConstraints;
        Console.WriteLine("\nFinal Solution:");
        double optimalValue = 0;
        HashSet<int> basicSet = new HashSet<int>(basicVariables);

        // Display values for basic variables
        for (int i = 0; i < numConstraints; i++)
        {
            int varIndex = basicVariables[i];
            double value = rhs[i];
            if (varIndex < originalVars)
            {
                Console.WriteLine($"x{varIndex + 1} = {value:F2}");
            }
            else
            {
                Console.WriteLine($"s{varIndex - originalVars + 1} = {value:F2}");
            }
            optimalValue += objCoeffs[varIndex] * value; // Calculate optimal value
        }

        // Display values for non-basic slack variables (not in basic set)
        for (int i = originalVars; i < totalVars; i++)
        {
            if (!basicSet.Contains(i))
            {
                Console.WriteLine($"s{i - originalVars + 1} = 0.00");
            }
        }

        Console.WriteLine($"\nOptimal Value: {optimalValue:F2}");

        // Calculate and display shadow prices
        Console.WriteLine("\nShadow Prices (Dual Variables):");
        for (int i = 0; i < numConstraints; i++)
        {
            double shadowPrice = objCoeffs[basicVariables[i]];
            Console.WriteLine($"Dual variable for constraint {i + 1}: {shadowPrice:F2}");
        }

        // Calculate and display slack variables correctly
        Console.WriteLine("\nSlack Variables:");
        for (int i = 0; i < numConstraints; i++)
        {
            double slack = rhs[i];
            if (basicVariables[i] < originalVars)
            {
                // The slack is zero for basic variables corresponding to original variables
                slack = 0;
            }
            Console.WriteLine($"Slack for constraint {i + 1}: {slack:F2}");
        }
    }

    // Shared function to display the current tableau
    static void PrintCurrentDictionaryForPrimal(double[] objectiveCoefficients, double[,] lhs, double[] rhs, int[] basicVariables, int totalVars, int numConstraints)
    {
        int originalVars = totalVars - numConstraints;

        Console.WriteLine("\nCurrent Dictionary:");
        Console.WriteLine("Objective Coefficients:");
        for (int j = 0; j < totalVars; j++)
        {
            if (j < originalVars)
                Console.Write($"x{j + 1}: {objectiveCoefficients[j]:F2} ");
            else
                Console.Write($"s{j - originalVars + 1}: {objectiveCoefficients[j]:F2} ");
        }
        Console.WriteLine();

        Console.WriteLine("\nConstraints:");
        for (int i = 0; i < numConstraints; i++)
        {
            Console.Write($"Basic Variable (Row {i + 1}): ");
            if (basicVariables[i] < originalVars)
                Console.Write($"x{basicVariables[i] + 1} = ");
            else
                Console.Write($"s{basicVariables[i] - originalVars + 1} = ");

            Console.Write($"{rhs[i]:F2} ");

            for (int j = 0; j < totalVars; j++)
            {
                if (lhs[i, j] != 0)
                {
                    if (j < originalVars)
                        Console.Write($"+ ({lhs[i, j]:F2}) x{j + 1} ");
                    else
                        Console.Write($"+ ({lhs[i, j]:F2}) s{j - originalVars + 1} ");
                }
            }
            Console.WriteLine();
        }
        Console.WriteLine();
    }


    static void SolveUsingDualSimplex(double[,] lhs, double[] rhs, double[] objectiveCoefficients, int totalVars, int numConstraints)
    {
        Console.WriteLine("\n--- Starting Dual Simplex Method ---");
        bool shouldContinue = true;
        int iterationCount = 0;
        int originalVars = totalVars - numConstraints;

        // Copy the original objective coefficients
        double[] objCoeffs = new double[totalVars];
        Array.Copy(objectiveCoefficients, objCoeffs, totalVars);

        // Initialize basic variables (slack variables should be basic initially)
        int[] basicVariables = new int[numConstraints];
        for (int i = 0; i < numConstraints; i++)
        {
            basicVariables[i] = originalVars + i; // Indices of slack variables
        }

        while (shouldContinue)
        {
            iterationCount++;
            Console.WriteLine($"\n--- Iteration {iterationCount} ---");

            // Step 1: Identify the most negative RHS entry (leaving variable)
            int leavingVarIndex = -1;
            double mostNegativeRHS = 0;
            for (int i = 0; i < numConstraints; i++)
            {
                if (rhs[i] < mostNegativeRHS)
                {
                    mostNegativeRHS = rhs[i];
                    leavingVarIndex = i;
                }
            }

            // If no negative RHS, optimal solution is found
            if (leavingVarIndex == -1)
            {
                Console.WriteLine("Optimal solution found.");
                shouldContinue = false;
                break;
            }

            // Step 2: Select entering variable by dual feasibility ratio test
            int enteringVarIndex = -1;
            double minRatio = double.MaxValue;
            for (int j = 0; j < totalVars; j++)
            {
                if (lhs[leavingVarIndex, j] < 0)  // Only consider negative entries in leaving row
                {
                    double reducedCost = objectiveCoefficients[j];
                    for (int i = 0; i < numConstraints; i++)
                    {
                        reducedCost -= objectiveCoefficients[basicVariables[i]] * lhs[i, j];
                    }
                    double ratio = reducedCost / lhs[leavingVarIndex, j];
                    if (ratio < minRatio)
                    {
                        minRatio = ratio;
                        enteringVarIndex = j;
                    }
                }
            }

            // If no entering variable is found, the problem is infeasible
            if (enteringVarIndex == -1)
            {
                Console.WriteLine("Problem is infeasible.");
                shouldContinue = false;
                break;
            }

            // Print chosen entering and leaving variables with corrected indexing
            Console.WriteLine($"Chosen leaving variable: {(basicVariables[leavingVarIndex] < originalVars ? "x" : "s")}{(basicVariables[leavingVarIndex] < originalVars ? basicVariables[leavingVarIndex] + 1 : basicVariables[leavingVarIndex] - originalVars + 1)}");
            Console.WriteLine($"Chosen entering variable: {(enteringVarIndex < originalVars ? "x" : "s")}{(enteringVarIndex < originalVars ? enteringVarIndex + 1 : enteringVarIndex - originalVars + 1)}");

            // Step 3: Perform pivot operation on the selected row and column
            double pivotElement = lhs[leavingVarIndex, enteringVarIndex];
            if (Math.Abs(pivotElement) < 1e-9)
            {
                Console.WriteLine("Pivot element is too small, leading to numerical instability.");
                shouldContinue = false;
                break;
            }

            // Normalize the pivot row
            for (int j = 0; j < totalVars; j++)
            {
                lhs[leavingVarIndex, j] /= pivotElement;
            }
            rhs[leavingVarIndex] /= pivotElement;

            // Update other rows to zero out the entering column
            for (int i = 0; i < numConstraints; i++)
            {
                if (i != leavingVarIndex)
                {
                    double factor = lhs[i, enteringVarIndex];
                    for (int j = 0; j < totalVars; j++)
                    {
                        lhs[i, j] -= factor * lhs[leavingVarIndex, j];
                    }
                    rhs[i] -= factor * rhs[leavingVarIndex];
                }
            }

            // Update the objective function coefficients
            double objectiveFactor = objectiveCoefficients[enteringVarIndex];
            for (int j = 0; j < totalVars; j++)
            {
                objectiveCoefficients[j] -= objectiveFactor * lhs[leavingVarIndex, j];
            }

            // Update the basic variables
            basicVariables[leavingVarIndex] = enteringVarIndex;

            // Print the current tableau
            PrintCurrentTableau1(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);
        }

        // Output the final solution
        OutputFinalSolutionWithShadowPricesAndDualSurplus(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);
        PerformSensitivityAnalysis(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);
    }

    // Utility method to print the current tableau for dual simplex
    static void PrintCurrentTableau1(double[] objectiveCoefficients, double[,] lhs, double[] rhs, int[] basicVariables, int totalVars, int numConstraints)
    {
        Console.WriteLine("\nCurrent Tableau:");
        Console.WriteLine("Objective Coefficients:");
        for (int j = 0; j < totalVars; j++)
        {
            Console.Write($"{(j < totalVars - numConstraints ? "x" : "s")}{(j < totalVars - numConstraints ? j + 1 : j - (totalVars - numConstraints) + 1)}: {objectiveCoefficients[j]:F2} ");
        }
        Console.WriteLine();

        Console.WriteLine("Constraints:");
        for (int i = 0; i < numConstraints; i++)
        {
            Console.Write($"Basic Variable (Row {i + 1}): {(basicVariables[i] < totalVars - numConstraints ? "x" : "s")}{(basicVariables[i] < totalVars - numConstraints ? basicVariables[i] + 1 : basicVariables[i] - (totalVars - numConstraints) + 1)} = {rhs[i]:F2} ");
            for (int j = 0; j < totalVars; j++)
            {
                Console.Write($"+ ({lhs[i, j]:F2}) {(j < totalVars - numConstraints ? "x" : "s")}{(j < totalVars - numConstraints ? j + 1 : j - (totalVars - numConstraints) + 1)} ");
            }
            Console.WriteLine();
        }
    }

    // Utility method to output the final solution, including shadow prices and dual surplus variables
    static void OutputFinalSolutionWithShadowPricesAndDualSurplus(double[] objCoeffs, double[,] lhs, double[] rhs, int[] basicVariables, int totalVars, int numConstraints)
    {
        Console.WriteLine("\nFinal Solution:");
        double optimalValue = 0;
        int originalVars = totalVars - numConstraints;
        HashSet<int> basicSet = new HashSet<int>(basicVariables);

        // Display values for basic variables and calculate the optimal value
        for (int i = 0; i < numConstraints; i++)
        {
            int varIndex = basicVariables[i];
            double value = rhs[i];
            Console.WriteLine($"{(varIndex < originalVars ? "x" : "s")}{(varIndex < originalVars ? varIndex + 1 : varIndex - originalVars + 1)} = {value:F2}");
            optimalValue += objCoeffs[varIndex] * value;
        }

        Console.WriteLine($"\nOptimal Value: {optimalValue:F2}");

        // Calculate and display shadow prices (dual variables) only for binding constraints
        Console.WriteLine("\nShadow Prices (Dual Variables):");
        for (int i = 0; i < numConstraints; i++)
        {
            // Calculate dual surplus to check if the constraint is binding
            double dualSurplus = rhs[i];
            for (int j = 0; j < totalVars; j++)
            {
                if (!basicSet.Contains(j)) // Only non-basic variables
                {
                    dualSurplus -= lhs[i, j] * objCoeffs[j];
                }
            }

            if (dualSurplus == 0) // Binding constraint
            {
                int basicVarIndex = basicVariables[i];
                double shadowPrice = objCoeffs[basicVarIndex];
                Console.WriteLine($"Shadow price for constraint {i + 1}: {shadowPrice:F2}");
            }
            else // Non-binding constraint
            {
                Console.WriteLine($"Shadow price for constraint {i + 1}: 0.00");
            }
        }

        // Display dual surplus variables for all constraints
        Console.WriteLine("\nDual Surplus Variables:");
        for (int i = 0; i < numConstraints; i++)
        {
            double dualSurplus = rhs[i];
            for (int j = 0; j < totalVars; j++)
            {
                if (!basicSet.Contains(j)) // Only non-basic variables
                {
                    dualSurplus -= lhs[i, j] * objCoeffs[j];
                }
            }
            Console.WriteLine($"Dual surplus for constraint {i + 1}: {dualSurplus:F2}");
        }
    }

    // Solve Using Bland's Rule
    static void SolveUsingBlandsRule(double[,] lhs, double[] rhs,
                                     double[] originalObjCoeffs, int totalVars, int numConstraints)
    {
        bool shouldContinue = true;
        int iterationCount = 0;

        // Copy the original objective coefficients
        double[] objCoeffs = new double[totalVars];
        double[] objectiveCoefficients = new double[totalVars]; // This will be modified during iterations

        Array.Copy(originalObjCoeffs, objCoeffs, totalVars); // Keep a copy of the original coefficients
        Array.Copy(originalObjCoeffs, objectiveCoefficients, totalVars); // This array will be updated

        // Initialize basic variables (slack variables)
        int[] basicVariables = new int[numConstraints];
        for (int i = 0; i < numConstraints; i++)
        {
            basicVariables[i] = totalVars - numConstraints + i; // Indices of slack variables
        }

        while (shouldContinue)
        {
            iterationCount++;
            Console.WriteLine($"Starting iteration {iterationCount} of Bland's Rule.");

            // Perform a single iteration of Bland's Rule
            shouldContinue = ApplyBlandsRule(lhs, rhs, objectiveCoefficients, totalVars, numConstraints, basicVariables, iterationCount);

            if (!shouldContinue)
            {
                Console.WriteLine("Optimal solution reached or no valid pivot available.");
                break;
            }

            // Print the current tableau
            PrintCurrentTableau1(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);
        }

        // Output the final solution
        OutputFinalSolutionWithShadowPricesAndSlacks(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);
        PerformSensitivityAnalysis(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);
    }

    static bool ApplyBlandsRule(double[,] lhs, double[] rhs, double[] objectiveCoefficients, int totalVars, int numConstraints, int[] basicVariables, int iterationCount)
    {
        // Print the current tableau (Dictionary)
        PrintCurrentTableau1(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);

        // Calculate the dual variables (pi)
        double[] pi = new double[numConstraints];
        for (int i = 0; i < numConstraints; i++)
        {
            pi[i] = objectiveCoefficients[basicVariables[i]]; // Use the basic variable's objective coefficient
        }

        // Use Bland's Rule to find the entering variable by objective coefficient
        int enteringVarIndex = -1;
        for (int j = 0; j < totalVars; j++)
        {
            double reducedCost = objectiveCoefficients[j];
            for (int i = 0; i < numConstraints; i++)
            {
                reducedCost -= pi[i] * lhs[i, j];
            }
            if (reducedCost > 0) // Check for positive reduced cost
            {
                enteringVarIndex = j; // Choose the smallest index with reducedCost > 0
                break; // Stop as we want the first (smallest index) variable with reducedCost > 0
            }
        }

        // If no entering variable is found, the current solution is optimal
        if (enteringVarIndex == -1)
        {
            Console.WriteLine("No entering variable found. Optimal solution may have been reached.");
            return false; // No more pivots required
        }

        // Print the chosen entering variable
        Console.WriteLine($"Chosen entering variable: x{enteringVarIndex + 1}");

        // Perform the ratio test to find the leaving variable
        int leavingVarIndex = -1;
        double minRatio = double.MaxValue;
        for (int i = 0; i < numConstraints; i++)
        {
            if (lhs[i, enteringVarIndex] > 0) // Positive coefficient for the entering variable
            {
                double ratio = rhs[i] / lhs[i, enteringVarIndex];
                if (ratio < minRatio)
                {
                    minRatio = ratio;
                    leavingVarIndex = i;
                }
                else if (ratio == minRatio) // Apply Bland's Rule to choose the leaving variable
                {
                    if (leavingVarIndex == -1 || basicVariables[i] < basicVariables[leavingVarIndex])
                    {
                        leavingVarIndex = i; // Choose the smallest index for the leaving variable
                    }
                }
            }
        }

        // If no leaving variable is found, the problem is unbounded
        if (leavingVarIndex == -1)
        {
            Console.WriteLine("No leaving variable found. The problem may be unbounded.");
            return false; // No pivot possible, likely an unbounded problem
        }

        // Print the chosen leaving variable
        Console.WriteLine($"Chosen leaving variable: x{basicVariables[leavingVarIndex] + 1}");

        // Perform the pivot operation
        double pivotElement = lhs[leavingVarIndex, enteringVarIndex];

        // Update the leaving row
        for (int j = 0; j < totalVars; j++)
        {
            lhs[leavingVarIndex, j] /= pivotElement;
        }
        rhs[leavingVarIndex] /= pivotElement;

        // Update all other rows
        for (int i = 0; i < numConstraints; i++)
        {
            if (i != leavingVarIndex)
            {
                double factor = lhs[i, enteringVarIndex];
                for (int j = 0; j < totalVars; j++)
                {
                    lhs[i, j] -= factor * lhs[leavingVarIndex, j];
                }
                rhs[i] -= factor * rhs[leavingVarIndex];
            }
        }

        // Update the objective coefficients in-place using the original objective coefficient
        double objectiveFactor = objectiveCoefficients[enteringVarIndex];
        for (int j = 0; j < totalVars; j++)
        {
            objectiveCoefficients[j] -= objectiveFactor * lhs[leavingVarIndex, j];
        }

        // Update the basic variable index
        basicVariables[leavingVarIndex] = enteringVarIndex;

        // Print updated objective coefficients and RHS after pivot
        Console.WriteLine("\nUpdated Objective Function:");
        for (int j = 0; j < totalVars; j++)
        {
            Console.Write($"{objectiveCoefficients[j]:F4} ");
        }
        Console.WriteLine();

        // Print BFS after pivot
        Console.WriteLine("\nUpdated Basic Feasible Solution (BFS) at this iteration:");
        for (int i = 0; i < numConstraints; i++)
        {
            Console.WriteLine($"x{basicVariables[i] + 1} = {rhs[i]:F4}");
        }

        // Output pivot operation
        Console.WriteLine($"Pivot completed: Entering variable x{enteringVarIndex + 1}, Leaving variable is x{basicVariables[leavingVarIndex] + 1}.\n");

        // Print the updated tableau
        PrintCurrentTableau1(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);

        return true; // Successful pivot operation, continue the process
    }


    // Two-Phase Method
    static void TwoPhaseMethod(double[] objCoeffs, double[,] lhs, double[] rhs, int numVars, int numConstraints)
    {
        Console.WriteLine("\n--- Starting Two-Phase Method ---");

        // Record start time
        var watch = Stopwatch.StartNew();

        // Phase 1: Add artificial variables
        int totalVars = numVars + numConstraints; // Original variables + slack variables
        int totalVarsPhase1 = totalVars + numConstraints; // + artificial variables

        double[] cPhase1 = new double[totalVarsPhase1];
        double[,] lhsPhase1 = new double[numConstraints, totalVarsPhase1];
        double[] rhsPhase1 = new double[numConstraints];

        // Objective coefficients for artificial variables in Phase 1 (Minimize sum of artificial variables)
        for (int j = totalVars; j < totalVarsPhase1; j++)
        {
            cPhase1[j] = -1; // Maximization problem
        }

        // Copy original lhs and rhs
        for (int i = 0; i < numConstraints; i++)
        {
            rhsPhase1[i] = rhs[i];
            for (int j = 0; j < totalVars; j++)
            {
                lhsPhase1[i, j] = lhs[i, j];
            }
            // Add artificial variables coefficients
            lhsPhase1[i, totalVars + i] = 1;
        }

        // Initialize basic variables (artificial variables)
        int[] basicVariables = new int[numConstraints];
        for (int i = 0; i < numConstraints; i++)
        {
            basicVariables[i] = totalVars + i; // Indices of artificial variables
        }

        // Start Phase 1
        bool feasible = true;
        while (true)
        {
            // Compute reduced costs
            double[] pi = new double[numConstraints];
            for (int i = 0; i < numConstraints; i++)
            {
                pi[i] = cPhase1[basicVariables[i]];
            }

            double[] reducedCosts = new double[totalVarsPhase1];
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                reducedCosts[j] = cPhase1[j];
                for (int i = 0; i < numConstraints; i++)
                {
                    reducedCosts[j] -= pi[i] * lhsPhase1[i, j];
                }
            }

            // Find entering variable (most positive reduced cost)
            int enteringVarIndex = -1;
            double maxReducedCost = 0;
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                if (reducedCosts[j] > maxReducedCost)
                {
                    maxReducedCost = reducedCosts[j];
                    enteringVarIndex = j;
                }
            }

            // If no entering variable, optimal for Phase 1
            if (enteringVarIndex == -1)
            {
                // Check if the artificial variables are zero
                double objValuePhase1 = 0;
                for (int i = 0; i < numConstraints; i++)
                {
                    objValuePhase1 += cPhase1[basicVariables[i]] * rhsPhase1[i];
                }
                if (Math.Abs(objValuePhase1) > 1e-6)
                {
                    feasible = false;
                }
                break;
            }

            // Perform ratio test to find leaving variable
            int leavingVarIndex = -1;
            double minRatio = double.MaxValue;
            for (int i = 0; i < numConstraints; i++)
            {
                if (lhsPhase1[i, enteringVarIndex] > 1e-6)
                {
                    double ratio = rhsPhase1[i] / lhsPhase1[i, enteringVarIndex];
                    if (ratio < minRatio)
                    {
                        minRatio = ratio;
                        leavingVarIndex = i;
                    }
                }
            }

            // If no leaving variable, problem is unbounded
            if (leavingVarIndex == -1)
            {
                Console.WriteLine("Problem is unbounded in Phase 1.");
                feasible = false;
                break;
            }

            // Pivot operation
            double pivotElement = lhsPhase1[leavingVarIndex, enteringVarIndex];

            // Update the leaving row
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                lhsPhase1[leavingVarIndex, j] /= pivotElement;
            }
            rhsPhase1[leavingVarIndex] /= pivotElement;

            // Update all other rows
            for (int i = 0; i < numConstraints; i++)
            {
                if (i != leavingVarIndex)
                {
                    double factor = lhsPhase1[i, enteringVarIndex];
                    for (int j = 0; j < totalVarsPhase1; j++)
                    {
                        lhsPhase1[i, j] -= factor * lhsPhase1[leavingVarIndex, j];
                    }
                    rhsPhase1[i] -= factor * rhsPhase1[leavingVarIndex];
                }
            }

            // Update the objective coefficients
            double objectiveFactor = cPhase1[enteringVarIndex];
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                cPhase1[j] -= objectiveFactor * lhsPhase1[leavingVarIndex, j];
            }

            // Update basic variables
            basicVariables[leavingVarIndex] = enteringVarIndex;
        }

        if (!feasible)
        {
            Console.WriteLine("The problem is infeasible.");
            return;
        }

        // Remove artificial variables from the basis
        for (int i = 0; i < numConstraints; i++)
        {
            if (basicVariables[i] >= totalVars)
            {
                // Artificial variable is in the basis
                // Try to remove it by finding a non-artificial variable with non-zero coefficient
                int enteringVarIndex = -1;
                for (int j = 0; j < totalVars; j++)
                {
                    if (Math.Abs(lhsPhase1[i, j]) > 1e-6)
                    {
                        enteringVarIndex = j;
                        break;
                    }
                }
                if (enteringVarIndex == -1)
                {
                    Console.WriteLine("Cannot remove artificial variable from basis. Problem is infeasible.");
                    feasible = false;
                    break;
                }
                else
                {
                    // Perform pivot operation to replace artificial variable
                    double pivotElement = lhsPhase1[i, enteringVarIndex];

                    // Normalize pivot row
                    for (int j = 0; j < totalVarsPhase1; j++)
                    {
                        lhsPhase1[i, j] /= pivotElement;
                    }
                    rhsPhase1[i] /= pivotElement;

                    // Update other rows
                    for (int k = 0; k < numConstraints; k++)
                    {
                        if (k != i)
                        {
                            double factor = lhsPhase1[k, enteringVarIndex];
                            for (int j = 0; j < totalVarsPhase1; j++)
                            {
                                lhsPhase1[k, j] -= factor * lhsPhase1[i, j];
                            }
                            rhsPhase1[k] -= factor * rhsPhase1[i];
                        }
                    }

                    // Update the objective function coefficients
                    double objFactor = cPhase1[enteringVarIndex];
                    for (int j = 0; j < totalVarsPhase1; j++)
                    {
                        cPhase1[j] -= objFactor * lhsPhase1[i, j];
                    }

                    // Update basic variable
                    basicVariables[i] = enteringVarIndex;
                }
            }
        }

        if (!feasible)
        {
            Console.WriteLine("The problem is infeasible.");
            return;
        }

        // Prepare for Phase 2 by removing artificial variables
        double[,] lhsPhase2 = new double[numConstraints, totalVars];
        double[] cPhase2 = new double[totalVars];
        double[] rhsPhase2 = new double[numConstraints];

        // Copy the relevant parts from Phase 1 arrays
        for (int i = 0; i < numConstraints; i++)
        {
            rhsPhase2[i] = rhsPhase1[i];
            for (int j = 0; j < totalVars; j++)
            {
                lhsPhase2[i, j] = lhsPhase1[i, j];
            }
        }

        // Set the original objective coefficients
        for (int j = 0; j < totalVars; j++)
        {
            cPhase2[j] = objCoeffs[j];
        }

        // Ensure basic variables are within the correct range
        for (int i = 0; i < numConstraints; i++)
        {
            if (basicVariables[i] >= totalVars)
            {
                Console.WriteLine("Artificial variable still in basis after removal. Problem is infeasible.");
                feasible = false;
                break;
            }
        }

        if (!feasible)
        {
            Console.WriteLine("The problem is infeasible.");
            return;
        }

        // Proceed with Phase 2 using the updated arrays
        lhsPhase1 = lhsPhase2;
        rhsPhase1 = rhsPhase2;
        cPhase1 = cPhase2;
        totalVarsPhase1 = totalVars; // Update totalVarsPhase1 to reflect the removal

        // Start Phase 2
        while (true)
        {
            // Compute reduced costs
            double[] pi = new double[numConstraints];
            for (int i = 0; i < numConstraints; i++)
            {
                pi[i] = cPhase1[basicVariables[i]];
            }

            double[] reducedCosts = new double[totalVarsPhase1];
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                reducedCosts[j] = cPhase1[j];
                for (int i = 0; i < numConstraints; i++)
                {
                    reducedCosts[j] -= pi[i] * lhsPhase1[i, j];
                }
            }

            // Find entering variable (most positive reduced cost)
            int enteringVarIndex = -1;
            double maxReducedCost = 0;
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                if (reducedCosts[j] > maxReducedCost)
                {
                    maxReducedCost = reducedCosts[j];
                    enteringVarIndex = j;
                }
            }

            // If no entering variable, optimal for Phase 2
            if (enteringVarIndex == -1)
            {
                break;
            }

            // Perform ratio test to find leaving variable
            int leavingVarIndex = -1;
            double minRatio = double.MaxValue;
            for (int i = 0; i < numConstraints; i++)
            {
                if (lhsPhase1[i, enteringVarIndex] > 1e-6)
                {
                    double ratio = rhsPhase1[i] / lhsPhase1[i, enteringVarIndex];
                    if (ratio < minRatio)
                    {
                        minRatio = ratio;
                        leavingVarIndex = i;
                    }
                }
            }

            // If no leaving variable, problem is unbounded
            if (leavingVarIndex == -1)
            {
                Console.WriteLine("Problem is unbounded in Phase 2.");
                return;
            }

            // Pivot operation
            double pivotElement = lhsPhase1[leavingVarIndex, enteringVarIndex];

            // Update the leaving row
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                lhsPhase1[leavingVarIndex, j] /= pivotElement;
            }
            rhsPhase1[leavingVarIndex] /= pivotElement;

            // Update all other rows
            for (int i = 0; i < numConstraints; i++)
            {
                if (i != leavingVarIndex)
                {
                    double factor = lhsPhase1[i, enteringVarIndex];
                    for (int j = 0; j < totalVarsPhase1; j++)
                    {
                        lhsPhase1[i, j] -= factor * lhsPhase1[leavingVarIndex, j];
                    }
                    rhsPhase1[i] -= factor * rhsPhase1[leavingVarIndex];
                }
            }

            // Update the objective coefficients
            double objectiveFactor = cPhase1[enteringVarIndex];
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                cPhase1[j] -= objectiveFactor * lhsPhase1[leavingVarIndex, j];
            }

            // Update basic variables
            basicVariables[leavingVarIndex] = enteringVarIndex;
        }

        // Output the final solution
        Console.WriteLine("Optimal Solution:");
        double[] solution = new double[totalVarsPhase1];

        for (int i = 0; i < numConstraints; i++)
        {
            solution[basicVariables[i]] = rhsPhase1[i];
        }

        for (int j = 0; j < totalVars; j++)
        {
            if (Math.Abs(solution[j]) > 1e-6)
            {
                Console.WriteLine($"x{j + 1} = {solution[j]:F4}");
            }
        }

        // Compute the objective value
        double objectiveValue = 0;
        for (int j = 0; j < totalVars; j++)
        {
            objectiveValue += objCoeffs[j] * solution[j];
        }
        Console.WriteLine($"Objective Value: {objectiveValue:F4}");

        // Record end time
        watch.Stop();
        TimeSpan ts = watch.Elapsed;
        Console.WriteLine($"Total runtime for Two-Phase Method: {ts.TotalSeconds} seconds");
    }

    // Solve Using Interior Point Method
    // Implementing the Path-Following Method
    static void SolveUsingPathFollowing(double[,] A, double[] b, double[] c, int numVars, int numConstraints)
    {
        Console.WriteLine("\n--- Starting Path-Following Method ---");

        int maxIterations = 100;
        double tolerance = 1e-6;
        double d = 0.15; // This is delta , (but can't be named delta )Parameter controlling convergence rate
        double gamma = 0.85;  // Duality measure

        double theta = 1.0;
        double r = 0.965; // Parameter less than 1 to ensure steps stay within the feasible region

        // Initialize variables
        int n = numVars;           // Number of variables
        int m = numConstraints;    // Number of constraints

        double[] x = new double[n];
        double[] omega = new double[m];
        double[] y = new double[m];
        double[] z = new double[n];

        // Start with feasible positive values
        for (int i = 0; i < n; i++)
        {
            x[i] = 200.0;
            z[i] = 800.0;
        }
        for (int i = 0; i < m; i++)
        {
            omega[i] = 650.0;
            y[i] = 400.0;
        }

        // Start timing
        Stopwatch stopwatch = new Stopwatch();
        stopwatch.Start();

        for (int iteration = 0; iteration < maxIterations; iteration++)
        {
            // Compute residuals
            // ρ = b - A x - ω
            double[] rho = new double[m];
            for (int i = 0; i < m; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                {
                    sum += A[i, j] * x[j];
                }
                rho[i] = b[i] - sum - omega[i];
            }

            // σ = c - A^T y + z
            double[] sigma = new double[n];
            for (int j = 0; j < n; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < m; i++)
                {
                    sum += A[i, j] * y[i];
                }
                sigma[j] = c[j] - sum + z[j];
            }

            // γ = z^T x + y^T ω
            gamma = 0.0;
            for (int j = 0; j < n; j++)
            {
                gamma += z[j] * x[j];
            }
            for (int i = 0; i < m; i++)
            {
                gamma += y[i] * omega[i];
            }

            // Duality measure μ = δ * γ / (n + m)
            double mu = d * gamma / (n + m);

            // Check for convergence
            if (mu < tolerance)
            {
                Console.WriteLine("Convergence achieved.");
                break;
            }

            // Compute RHS for equations (3) and (4)
            double[] rhs3 = new double[n];
            double[] rhs4 = new double[m];
            for (int j = 0; j < n; j++)
            {
                rhs3[j] = mu - x[j] * z[j];
            }
            for (int i = 0; i < m; i++)
            {
                rhs4[i] = mu - y[i] * omega[i];
            }

            // Build the linear system to solve for Δx, Δω, Δy, Δz
            int size = n + m + n + m; // Δx, Δω, Δy, Δz
            double[,] M = new double[size, size];
            double[] RHS = new double[size];

            // Equation 1: A Δx + Δω = ρ
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    M[i, j] = A[i, j]; // A
                }
                M[i, n + i] = 1.0; // Δω
                RHS[i] = rho[i];
            }

            // Equation 2: A^T Δy - Δz = σ
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    M[m + j, n + m + i] = A[i, j]; // A^T
                }
                M[m + j, n + m + m + j] = -1.0; // -Δz
                RHS[m + j] = sigma[j];
            }

            // Equation 3: Z Δx + X Δz = rhs3
            for (int j = 0; j < n; j++)
            {
                M[m + n + j, j] = z[j]; // Z
                M[m + n + j, n + m + m + j] = x[j]; // X
                RHS[m + n + j] = rhs3[j];
            }

            // Equation 4: W Δy + Y Δω = rhs4
            for (int i = 0; i < m; i++)
            {
                M[m + n + n + i, n + m + i] = omega[i]; // W
                M[m + n + n + i, n + i] = y[i]; // Y
                RHS[m + n + n + i] = rhs4[i];
            }

            // Solve the linear system M * Δ = RHS
            double[] delta = SolveLinearSystem(M, RHS);

            if (delta == null)
            {
                Console.WriteLine("Failed to solve linear system. Terminating path-following method.");
                return;
            }

            // Extract Δx, Δω, Δy, Δz
            double[] deltaX = new double[n];
            double[] deltaOmega = new double[m];
            double[] deltaY = new double[m];
            double[] deltaZ = new double[n];

            int index = 0;
            for (int j = 0; j < n; j++)
            {
                deltaX[j] = delta[index++];
            }
            for (int i = 0; i < m; i++)
            {
                deltaOmega[i] = delta[index++];
            }
            for (int i = 0; i < m; i++)
            {
                deltaY[i] = delta[index++];
            }
            for (int j = 0; j < n; j++)
            {
                deltaZ[j] = delta[index++];
            }

            // Compute step length θ


            for (int i = 0; i < n; i++)
            {
                if (deltaX[i] < 0)
                {
                    double ratio = -x[i] / deltaX[i];
                    if (ratio < theta)
                    {
                        theta = ratio;
                    }
                }
                if (deltaZ[i] < 0)
                {
                    double ratio = -z[i] / deltaZ[i];
                    if (ratio < theta)
                    {
                        theta = ratio;
                    }
                }
            }
            for (int i = 0; i < m; i++)
            {
                if (deltaY[i] < 0)
                {
                    double ratio = -y[i] / deltaY[i];
                    if (ratio < theta)
                    {
                        theta = ratio;
                    }
                }
                if (deltaOmega[i] < 0)
                {
                    double ratio = -omega[i] / deltaOmega[i];
                    if (ratio < theta)
                    {
                        theta = ratio;
                    }
                }
            }

            // Multiply theta by r
            theta = Math.Min(theta * r, 1.0);

            // Update variables
            for (int i = 0; i < n; i++)
            {
                x[i] += theta * deltaX[i];
                z[i] += theta * deltaZ[i];
                if (x[i] <= 0) x[i] = 1e-8; // Ensure positivity
                if (z[i] <= 0) z[i] = 1e-8;
            }
            for (int i = 0; i < m; i++)
            {
                y[i] += theta * deltaY[i];
                omega[i] += theta * deltaOmega[i];
                if (y[i] <= 0) y[i] = 1e-8;
                if (omega[i] <= 0) omega[i] = 1e-8;
            }

            // Optionally, print iteration info
            Console.WriteLine($"Iteration {iteration + 1}, Duality Measure μ = {mu:E6}, Step Length θ = {theta:F4}");
        }

        // Stop timing
        stopwatch.Stop();

        // Compute objective value c^T x
        double objectiveValue = 0.0;
        for (int j = 0; j < n; j++)
        {
            objectiveValue += c[j] * x[j];
        }

        // Output results
        Console.WriteLine("\nOptimal Solution (Path-Following Method):");
        for (int j = 0; j < n; j++)
        {
            Console.WriteLine($"x{j + 1} = {x[j]:F6}");
        }
        for (int i = 0; i < m; i++)
        {
            Console.WriteLine($"ω{i + 1} = {omega[i]:F6}");
        }
        for (int i = 0; i < m; i++)
        {
            Console.WriteLine($"y{i + 1} = {y[i]:F6}");
        }
        for (int j = 0; j < n; j++)
        {
            Console.WriteLine($"z{j + 1} = {z[j]:F6}");
        }
        Console.WriteLine($"Objective Value: {objectiveValue:F6}");
        Console.WriteLine($"Runtime for Path-Following Method: {stopwatch.Elapsed.TotalSeconds:F6} seconds");
    }

    // SolveLinearSystem using Gaussian Elimination with Partial Pivoting
    static double[] SolveLinearSystem(double[,] A, double[] b)
    {
        int n = b.Length;
        double[] x = new double[n];
        double[,] Aug = new double[n, n + 1];

        // Build augmented matrix
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Aug[i, j] = A[i, j];
            }
            Aug[i, n] = b[i];
        }

        // Gaussian elimination with partial pivoting
        for (int k = 0; k < n; k++)
        {
            // Find the pivot row
            double max = Math.Abs(Aug[k, k]);
            int pivotRow = k;
            for (int i = k + 1; i < n; i++)
            {
                if (Math.Abs(Aug[i, k]) > max)
                {
                    max = Math.Abs(Aug[i, k]);
                    pivotRow = i;
                }
            }

            // Swap rows if necessary
            if (pivotRow != k)
            {
                for (int j = 0; j < n + 1; j++)
                {
                    double temp = Aug[k, j];
                    Aug[k, j] = Aug[pivotRow, j];
                    Aug[pivotRow, j] = temp;
                }
            }

            // Check for zero pivot
            if (Math.Abs(Aug[k, k]) < 1e-12)
            {
                Console.WriteLine("Warning: Matrix is singular or nearly singular.");
                return null;
            }

            // Eliminate entries below pivot
            for (int i = k + 1; i < n; i++)
            {
                double factor = Aug[i, k] / Aug[k, k];
                for (int j = k; j < n + 1; j++)
                {
                    Aug[i, j] -= factor * Aug[k, j];
                }
            }
        }

        // Back substitution
        for (int i = n - 1; i >= 0; i--)
        {
            x[i] = Aug[i, n];
            for (int j = i + 1; j < n; j++)
            {
                x[i] -= Aug[i, j] * x[j];
            }
            x[i] /= Aug[i, i];
        }

        return x;
    }

    // Method to print the current dictionary (tableau)
    static void PrintCurrentDictionary(double[] objectiveCoefficients, double[,] lhs, double[] rhs, int[] basicVariables, int totalVars, int numConstraints)
    {
        Console.WriteLine("Current Dictionary:");

        // Print objective function
        Console.Write("η = ");
        bool first = true;
        for (int j = 0; j < totalVars; j++)
        {
            if (Math.Abs(objectiveCoefficients[j]) > 1e-6)
            {
                if (!first && objectiveCoefficients[j] > 0)
                    Console.Write("+ ");
                Console.Write($"{objectiveCoefficients[j]:F2}x{j + 1} ");
                first = false;
            }
        }
        Console.WriteLine();

        // Print constraints
        for (int i = 0; i < numConstraints; i++)
        {
            Console.Write($"x{basicVariables[i] + 1} = {rhs[i]:F2} ");
            for (int j = 0; j < totalVars; j++)
            {
                if (Math.Abs(lhs[i, j]) > 1e-6)
                {
                    if (lhs[i, j] > 0)
                        Console.Write("+ ");
                    Console.Write($"{-lhs[i, j]:F2}x{j + 1} ");
                }
            }
            Console.WriteLine();
        }
        Console.WriteLine();
    }

    // Method to output the final solution
    static void OutputFinalSolution(double[] originalObjCoeffs, double[,] lhs, double[] rhs, int[] basicVariables, int totalVars, int numConstraints)
    {
        Console.WriteLine("Optimal Solution:");
        double[] solution = new double[totalVars];

        for (int i = 0; i < numConstraints; i++)
        {
            solution[basicVariables[i]] = rhs[i];
        }

        for (int j = 0; j < totalVars; j++)
        {
            if (Math.Abs(solution[j]) > 1e-6)
            {
                Console.WriteLine($"x{j + 1} = {solution[j]:F4}");
            }
        }

        // Compute the objective value
        double objectiveValue = 0;
        for (int j = 0; j < totalVars; j++)
        {
            objectiveValue += originalObjCoeffs[j] * solution[j];
        }
        Console.WriteLine($"Objective Value: {objectiveValue:F4}");
    }

    // Determine which simplex method to use based on the conditions
    static string DetermineSimplexMethod(double[] rhs, double[] objCoeffs)
    {
        bool all_bi_nonnegative = true;
        bool all_cN_nonpositive = true;
        bool some_bi_zero = false;
        bool some_bi_nonnegative = false;
        bool some_bi_negative = false;

        foreach (double bi in rhs)
        {
            if (bi < 0) { all_bi_nonnegative = false; some_bi_negative = true; }
            if (bi >= 0) some_bi_nonnegative = true;
            if (bi == 0) some_bi_zero = true;
        }

        foreach (double cN in objCoeffs)
        {
            if (cN > 0) all_cN_nonpositive = false;
        }

        if (all_bi_nonnegative && all_cN_nonpositive)
        {
            return "Optimal";  // The problem is optimal
        }
        else if (some_bi_zero)
        {
            return "BlandsRule";  // Use Bland's Rule
        }
        else if (all_bi_nonnegative && !all_cN_nonpositive)
        {
            return "Primal";  // Use Primal Simplex
        }
        else if (some_bi_negative && all_cN_nonpositive)
        {
            return "Dual";  // Use Dual Simplex
        }
        else
        {
            return "TwoPhase";  // Use Two-phase method
        }
    }

    // Check for unboundedness
    static bool CheckUnboundedness(double[] objCoeffs, double[,] lhs, double[] rhs, int totalVars, int numConstraints)
    {
        // Check if any variable can increase the objective function indefinitely
        for (int j = 0; j < totalVars; j++)
        {
            if (objCoeffs[j] > 1e-6)
            {
                bool isBounded = false;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (lhs[i, j] > 1e-6)
                    {
                        isBounded = true;
                        break;
                    }
                }
                if (!isBounded)
                {
                    return true; // Unbounded in this direction
                }
            }
        }
        return false; // No unbounded direction found
    }

    static void PerformSensitivityAnalysis(double[] objectiveCoefficients, double[,] lhs, double[] rhs,
                                       int[] basicVariables, int totalVars, int numConstraints)
    {
        Console.WriteLine("\n--- Sensitivity Analysis ---");

        // Step 1: Compute the allowable range for RHS (right-hand side) constants
        Console.WriteLine("\nRange for RHS constants:");
        for (int i = 0; i < numConstraints; i++)
        {
            double minChange = double.MinValue;
            double maxChange = double.MaxValue;

            // Only consider positive entries in the basic columns for the i-th basic variable
            for (int j = 0; j < numConstraints; j++)
            {
                if (lhs[j, basicVariables[i]] != 0)
                {
                    double ratio = rhs[j] / lhs[j, basicVariables[i]];
                    if (lhs[j, basicVariables[i]] > 0)
                    {
                        maxChange = Math.Min(maxChange, ratio);
                    }
                    else if (lhs[j, basicVariables[i]] < 0)
                    {
                        minChange = Math.Max(minChange, ratio);
                    }
                }
            }

            Console.WriteLine($"Constraint {i + 1} (RHS): Allowable range = [{minChange}, {maxChange}]");
        }

        // Step 2: Compute the allowable range for objective function coefficients
        Console.WriteLine("\nRange for Objective Coefficients:");
        for (int j = 0; j < totalVars; j++)
        {
            if (Array.IndexOf(basicVariables, j) != -1) continue; // Skip basic variables

            double minChange = double.MinValue;
            double maxChange = double.MaxValue;

            // Iterate through all constraints to find the allowable range for each non-basic variable
            for (int i = 0; i < numConstraints; i++)
            {
                double reducedCost = objectiveCoefficients[j];
                for (int k = 0; k < numConstraints; k++)
                {
                    reducedCost -= objectiveCoefficients[basicVariables[k]] * lhs[k, j];
                }

                if (reducedCost < 0)
                {
                    maxChange = Math.Min(maxChange, -reducedCost);
                }
                else if (reducedCost > 0)
                {
                    minChange = Math.Max(minChange, -reducedCost);
                }
            }

            Console.WriteLine($"Objective coefficient for x{j + 1}: Allowable range = [{minChange}, {maxChange}]");
        }

        // Step 3: Compute the allowable range for coefficients in the LHS matrix
        Console.WriteLine("\nRange for LHS coefficients:");
        for (int i = 0; i < numConstraints; i++)
        {
            for (int j = 0; j < totalVars; j++)
            {
                if (Array.IndexOf(basicVariables, j) != -1) continue; // Skip coefficients associated with basic variables

                double minChange = double.MinValue;
                double maxChange = double.MaxValue;

                // Iterate through all constraints to find allowable range for each non-basic coefficient
                for (int k = 0; k < numConstraints; k++)
                {
                    if (lhs[k, j] == 0) continue;

                    double ratio = rhs[k] / lhs[k, j];
                    if (lhs[k, j] > 0)
                    {
                        maxChange = Math.Min(maxChange, ratio);
                    }
                    else
                    {
                        minChange = Math.Max(minChange, ratio);
                    }
                }

                Console.WriteLine($"Coefficient for LHS[{i}, {j}]: Allowable range = [{minChange}, {maxChange}]");
            }
        }

        Console.WriteLine("\n--- Sensitivity Analysis Completed ---");
    }


    static void Main(string[] args)
    {
        // -----------------Data Processing-----------------

        // Read CSV file
        var data = File.ReadAllLines("data.csv"); //699.25, 699.246(680), r=0.965(0.938)
        //var data = File.ReadAllLines("deg z=21 x1=3 x2=3.csv"); //21, 20.998(14.29), r=0.965(0.938)
        //var data = File.ReadAllLines("20-Klee-Minty Instance.csv");
        //var data = File.ReadAllLines("simplex_obj13_x1_2_x2_0_x3_1.csv"); //13, 12.93(4.31), r=0.965(0.938)
        //var data = File.ReadAllLines("dual_obj_-7_x1_7_x2_0_x3_0.csv"); //-7, -6.76, r=0.938

        // Determine number of variables and constraints based on data
        int numVars = data[0].Split(',').Length - 3;
        int numConstraints = data.Length - 1;
        int totalVars = numVars + numConstraints; // Including slack variables

        objCoeffs = new double[totalVars];
        double[] rhs = new double[numConstraints];
        double[,] lhs = new double[numConstraints, totalVars]; // Adjusted for slack variables

        // Initialize variables to track constraints
        int constraintIndex = 0;

        // Parse CSV file
        for (int i = 0; i < data.Length; i++)
        {
            var line = data[i].Split(',');

            // Check if this is the objective function row
            if (line[1].Trim() == "1")
            {
                // This row is the objective function
                for (int j = 2; j < line.Length - 1; j++)
                {
                    objCoeffs[j - 2] = double.Parse(line[j]);
                }
                // Coefficients for slack variables are zero
                for (int j = numVars; j < totalVars; j++)
                {
                    objCoeffs[j] = 0;
                }
            }
            else
            {
                // This row is a constraint
                for (int j = 2; j < line.Length - 1; j++)
                {
                    lhs[constraintIndex, j - 2] = double.Parse(line[j]);
                }
                // Coefficients for slack variables
                lhs[constraintIndex, numVars + constraintIndex] = 1; // Slack variable coefficient
                // RHS value
                rhs[constraintIndex] = double.Parse(line[^1]);
                constraintIndex++;
            }
        }

        // ----------------- Model -----------------

        // Unboundedness check
        if (CheckUnboundedness(objCoeffs, lhs, rhs, totalVars, numConstraints))
        {
            Console.WriteLine("The problem is unbounded.");
            return;
        }
        else
        {
            Console.WriteLine("The problem is bounded.");
        }

        // Determine which simplex method to use
        string method = DetermineSimplexMethod(rhs, objCoeffs);
        Console.WriteLine($"Using the {method} method.");

        // Measure performance for the chosen method
        Stopwatch stopwatchChosenMethod = new Stopwatch();
        stopwatchChosenMethod.Start();

        // Before calling the method, make a copy of the original objective coefficients
        double[] originalObjCoeffs = new double[totalVars];
        for (int j = 0; j < totalVars; j++)
        {
            originalObjCoeffs[j] = objCoeffs[j];
        }

        // Solve based on the chosen method
        switch (method)
        {
            case "Optimal":
                Console.WriteLine("The initial solution is optimal.");
                break;

            case "BlandsRule":
                SolveUsingBlandsRule(lhs, rhs, originalObjCoeffs, totalVars, numConstraints);
                break;

            case "Primal":
                SolveUsingPrimalSimplex(lhs, rhs, objCoeffs, totalVars, numConstraints);
                break;

            case "Dual":
                SolveUsingDualSimplex(lhs, rhs, objCoeffs, totalVars, numConstraints);
                break;

            case "TwoPhase":
                TwoPhaseMethod(objCoeffs, lhs, rhs, numVars, numConstraints);
                break;
        }

        // Total run time for the chosen method
        stopwatchChosenMethod.Stop();
        Console.WriteLine($"Total runtime for {method} method: {stopwatchChosenMethod.Elapsed.TotalSeconds} seconds");

        // Measure performance for the Path-Following method
        Stopwatch stopwatchPathFollowing = new Stopwatch();
        stopwatchPathFollowing.Start();

        if (method == "Optimal")
        {
            Console.WriteLine("optimal");

        }
        else
        {
            SolveUsingPathFollowing(lhs, rhs, originalObjCoeffs, numVars, numConstraints);
        }

        // Total run time for the Path-Following method
        stopwatchPathFollowing.Stop();
        Console.WriteLine($"Total runtime for Path-Following method: {stopwatchPathFollowing.Elapsed.TotalSeconds} seconds");

        // Compare runtimes
        Console.WriteLine("\n--- Runtime Comparison ---");
        Console.WriteLine($"{method} Method Runtime: {stopwatchChosenMethod.Elapsed.TotalSeconds} seconds");
        Console.WriteLine($"Path-Following Method Runtime: {stopwatchPathFollowing.Elapsed.TotalSeconds} seconds");
    }
}
