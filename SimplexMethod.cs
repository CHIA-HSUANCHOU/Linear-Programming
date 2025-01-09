/*
 Main Features:
 - Reads and parses LP data from a CSV file.
 - Determines the appropriate simplex method (Primal, Dual, Bland's Rule, Two-Phase).
 - Implements Bland's Rule to handle degeneracy and Two-Phase Method for infeasibility.
 - Tracks iteration count manually (without built-in methods).
 - Outputs all non-zero solution variables (including slack variables).

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
*/

using System;
using System.IO;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

class ResourceAllocationProblem
{
    static double[] objCoeffs;

    // Solve Using Primal Simplex Method
    static void SolveUsingPrimalSimplex(double[,] lhs, double[] rhs,
                                        double[] objectiveCoefficients, int totalVars, int numConstraints)
    {
        bool shouldContinue = true;
        int iterationCount = 0;

        // Store the original objective coefficients
        objCoeffs = new double[totalVars];
        for (int j = 0; j < totalVars; j++)
        {
            objCoeffs[j] = objectiveCoefficients[j];
        }

        // Initialize basic variables (slack variables)
        int[] basicVariables = new int[numConstraints];
        for (int i = 0; i < numConstraints; i++)
        {
            basicVariables[i] = totalVars - numConstraints + i; // Indices of slack variables
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

            // Identify non-basic variables
            HashSet<int> basicSet = new HashSet<int>(basicVariables);
            List<int> nonBasicVariables = new List<int>();
            for (int j = 0; j < totalVars; j++)
            {
                if (!basicSet.Contains(j))
                    nonBasicVariables.Add(j);
            }

            // Select entering variable (most positive reduced cost)
            int enteringVarIndex = -1;
            double maxReducedCost = 0;
            foreach (int j in nonBasicVariables)
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

            // Print the chosen entering variable
            Console.WriteLine($"Chosen entering variable: x{enteringVarIndex + 1}");

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

            // If no leaving variable is found, problem is unbounded
            if (leavingVarIndex == -1)
            {
                Console.WriteLine("Problem is unbounded.");
                shouldContinue = false;
                break;
            }

            // Print the chosen leaving variable
            Console.WriteLine($"Chosen leaving variable: x{basicVariables[leavingVarIndex] + 1}");

            // Perform pivot operation
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

            // Update the objective coefficients
            double objectiveFactor = objectiveCoefficients[enteringVarIndex];
            for (int j = 0; j < totalVars; j++)
            {
                objectiveCoefficients[j] -= objectiveFactor * lhs[leavingVarIndex, j];
            }

            // Update the basic variable index
            basicVariables[leavingVarIndex] = enteringVarIndex;

            // Print the current tableau
            PrintCurrentDictionary(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);
        }

        // Output the final solution
        OutputFinalSolution(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);
    }

    // Solve Using Dual Simplex Method
    static void SolveUsingDualSimplex(double[,] lhs, double[] rhs,
                                      double[] objectiveCoefficients, int totalVars, int numConstraints)
    {
        bool shouldContinue = true;
        int iterationCount = 0;

        // Store the original objective coefficients
        objCoeffs = new double[totalVars];
        for (int j = 0; j < totalVars; j++)
        {
            objCoeffs[j] = objectiveCoefficients[j];
        }

        // Initialize basic variables (slack variables)
        int[] basicVariables = new int[numConstraints];
        for (int i = 0; i < numConstraints; i++)
        {
            basicVariables[i] = totalVars - numConstraints + i; // Indices of slack variables
        }

        while (shouldContinue)
        {
            iterationCount++;
            Console.WriteLine($"\n--- Iteration {iterationCount} ---");

            // Identify leaving variable (most negative rhs)
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

            // If no leaving variable is found, the solution is optimal
            if (leavingVarIndex == -1)
            {
                Console.WriteLine("Optimal solution found.");
                shouldContinue = false;
                break;
            }

            // Print the chosen leaving variable
            Console.WriteLine($"Chosen leaving variable: x{basicVariables[leavingVarIndex] + 1}");

            // Compute ratios for entering variable selection
            int enteringVarIndex = -1;
            double minRatio = double.MaxValue;
            for (int j = 0; j < totalVars; j++)
            {
                if (lhs[leavingVarIndex, j] < 0)
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

            // If no entering variable is found, problem is infeasible
            if (enteringVarIndex == -1)
            {
                Console.WriteLine("Problem is infeasible.");
                shouldContinue = false;
                break;
            }

            // Print the chosen entering variable
            Console.WriteLine($"Chosen entering variable: x{enteringVarIndex + 1}");

            // Perform pivot operation
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

            // Update the objective coefficients
            double objectiveFactor = objectiveCoefficients[enteringVarIndex];
            for (int j = 0; j < totalVars; j++)
            {
                objectiveCoefficients[j] -= objectiveFactor * lhs[leavingVarIndex, j];
            }

            // Update the basic variable index
            basicVariables[leavingVarIndex] = enteringVarIndex;

            // Print the current tableau
            PrintCurrentDictionary(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);
        }

        // Output the final solution
        OutputFinalSolution(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);
    }

    // Solve Using Bland's Rule
    static void SolveUsingBlandsRule(double[,] lhs, double[] rhs,
                                     double[] originalObjCoeffs, int totalVars, int numConstraints)
    {
        bool shouldContinue = true;
        int iterationCount = 0;

        // Copy the original objective coefficients
        objCoeffs = new double[totalVars];
        double[] objectiveCoefficients = new double[totalVars]; // This will be modified during iterations

        for (int j = 0; j < totalVars; j++)
        {
            objCoeffs[j] = originalObjCoeffs[j]; // Keep a copy of the original coefficients
            objectiveCoefficients[j] = originalObjCoeffs[j]; // This array will be updated
        }

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

            shouldContinue = ApplyBlandsRule(lhs, rhs, objectiveCoefficients, totalVars, numConstraints, basicVariables, iterationCount);

            if (!shouldContinue)
            {
                Console.WriteLine("Optimal solution reached or no valid pivot available.");
                break;
            }

            // Print the current tableau
            PrintCurrentDictionary(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);
        }

        // Output the final solution
        OutputFinalSolution(objCoeffs, lhs, rhs, basicVariables, totalVars, numConstraints);
    }

    static bool ApplyBlandsRule(double[,] lhs, double[] rhs, double[] objectiveCoefficients, int totalVars, int numConstraints, int[] basicVariables, int iterationCount)
    {
        // Print the current tableau (Dictionary)
        PrintCurrentDictionary(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);

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

        // Update the objective coefficients
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
        PrintCurrentDictionary(objectiveCoefficients, lhs, rhs, basicVariables, totalVars, numConstraints);

        return true; // Successful pivot operation, continue the process
    }

    // Two-Phase Method
    static void TwoPhaseMethod(double[] objCoeffs, double[,] lhs, double[] rhs, int numVars, int numConstraints)
    {
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
            cPhase1[j] = -1; // Minimization problem
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

            // Find entering variable
            int enteringVarIndex = -1;
            double maxReducedCost = 0;
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                if (reducedCosts[j] < maxReducedCost)
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
                if (lhsPhase1[i, enteringVarIndex] > 0)
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

        // Phase 2: Remove artificial variables and proceed with original objective
        double[] cPhase2 = new double[totalVarsPhase1];
        for (int j = 0; j < totalVars; j++)
        {
            cPhase2[j] = objCoeffs[j];
        }
        // Objective coefficients for artificial variables are zero in Phase 2

        // Start Phase 2
        while (true)
        {
            // Compute reduced costs
            double[] pi = new double[numConstraints];
            for (int i = 0; i < numConstraints; i++)
            {
                pi[i] = cPhase2[basicVariables[i]];
            }

            double[] reducedCosts = new double[totalVarsPhase1];
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                reducedCosts[j] = cPhase2[j];
                for (int i = 0; i < numConstraints; i++)
                {
                    reducedCosts[j] -= pi[i] * lhsPhase1[i, j];
                }
            }

            // Find entering variable
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
                if (lhsPhase1[i, enteringVarIndex] > 0)
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
            double objectiveFactor = cPhase2[enteringVarIndex];
            for (int j = 0; j < totalVarsPhase1; j++)
            {
                cPhase2[j] -= objectiveFactor * lhsPhase1[leavingVarIndex, j];
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
            if (solution[j] != 0)
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
        Console.WriteLine($"Total runtime: {ts.TotalSeconds} seconds");
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

        foreach (double bi in rhs)
        {
            if (bi < 0) all_bi_nonnegative = false;
            if (bi > 0 || bi == 0) some_bi_nonnegative = true;
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
        else if (some_bi_nonnegative && all_cN_nonpositive)
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
        // Check if all RHS values are non-negative
        for (int i = 0; i < numConstraints; i++)
        {
            if (rhs[i] < 0) return false; // If any RHS is negative, problem may not be unbounded
        }

        // Check if we can increase the objective function indefinitely
        for (int j = 0; j < totalVars; j++)
        {
            // Look for a variable with a positive objective coefficient
            if (objCoeffs[j] > 0)
            {
                // Check if there is no constraint limiting this variable
                bool isBounded = false;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (lhs[i, j] > 0 && rhs[i] > 0)
                    {
                        isBounded = true; // There is a constraint that bounds this variable
                        break;
                    }
                }
                if (!isBounded) return true; // Unbounded in this direction
            }
        }
        return false; // No unbounded direction found
    }

    static void Main(string[] args)
    {
        // -----------------Data Processing-----------------

        // Read CSV file
        var data = File.ReadAllLines("20-Klee-Minty Instance.csv");

        // Determine number of variables and constraints based on data
        int numVars = data[0].Split(',').Length - 2;
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

        // Measure performance
        Stopwatch stopwatch = new Stopwatch();
        stopwatch.Start();

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

        // Total run time
        stopwatch.Stop();
        Console.WriteLine($"Total runtime: {stopwatch.Elapsed.TotalSeconds} seconds");
    }
}
