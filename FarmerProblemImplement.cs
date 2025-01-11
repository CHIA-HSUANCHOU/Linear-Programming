using Microsoft.VisualBasic;
using System;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Linq;
using System.Reflection;
using System.Security.Cryptography;
public static class ArrayExtensions
{
    public static double[] GetRow(this double[,] array, int row)
    {
        int cols = array.GetLength(1);
        double[] result = new double[cols];
        for (int col = 0; col < cols; col++)
        {
            result[col] = array[row, col];
        }
        return result;
    }
}


namespace FarmerLShapedNoILOG
{
    class Program
    {

        /*
         * Updated Farmer's Problem with L-Shaped Decomposition
         * 
         * Problem Context:
         * - Total Land: 500 acres
         * - Crops: Wheat, Corn, Sugarbeets
         * 
         * Crop Yields, Costs, and Prices:
         *   1. **Wheat**:
         *      - Yield: 2.5 T/acre
         *      - Planting Cost: $150/acre
         *      - Selling Price: $170/T
         *      - Purchase Price: $238/T
         *      - Minimum Demand: 200 T
         *   2. **Corn**:
         *      - Yield: 3 T/acre
         *      - Planting Cost: $230/acre
         *      - Selling Price: $150/T
         *      - Purchase Price: $210/T
         *      - Minimum Demand: 240 T
         *   3. **Sugarbeets**:
         *      - Yield: 20 T/acre
         *      - Planting Cost: $260/acre
         *      - Selling Price:
         *          - Up to 6000 T: $36/T
         *          - Above 6000 T: $10/T
         *      - No Minimum Demand or Purchase Requirement.
         * 
         * Variables:
         * - **First-Stage Variables** (Master Problem):
         *   - x_1: Land allocation for wheat
         *   - x_2: Land allocation for corn
         *   - x_3: Land allocation for sugarbeets
         *   - θ  : Second-stage cost approximation variable (L-shaped)
         *
         * - **Second-Stage Variables** (Subproblem):
         *   - y_1, y_2: Shortfall of wheat and corn (purchase)
         *   - w_1, w_2: Surplus of wheat and corn (sell)
         *   - w_3, w_4: Sugarbeet sales (two-tier pricing)
         *
         * Objective Function:
         * - Minimize the total cost including:
         *   - First-stage planting cost: `150x1 + 230x2 + 260x3`
         *   - Second-stage cost (shortfalls/surplus): `238y1 + 210y2 - 170w1 - 150w2 - 36w3 - 10w4`
         * 
         * Constraints:
         *  1. **C1**: Total land constraint: `x_1 + x_2 + x_3 <= 500`
         *  2. **C2**: Wheat yield shortfall: `2.5x_1 + y_1 - w_1 >= 200` → `-2.5x1 - y1 + w1 <= -200`
         *  3. **C3**: Corn yield shortfall: `3x_2 + y_2 - w_2 >= 240` → `-3x2 - y2 + w2 <= -240`
         *  4. **C4**: Sugarbeet selling capacity: `w_3 + w_4 <= 20x3` → `w3 + w4 - 20x3 <= 0`
         *  5. **C5**: Sugarbeet tier limit: `w3 <= 6000`
         * 
         * L-Shaped Decomposition:
         * - **Master Problem**:
         *   - Focuses on first-stage variables and uses θ to represent second-stage costs.
         *   - User selects constraints for inclusion.
         * - **Subproblem**:
         *   - Accounts for second-stage costs and feasibility under different scenarios.
         *   - Includes yield uncertainties ranging from **80% to 120%**.
         * 
         * Key Features:
         * - Implements primal simplex, dual simplex, and Bland’s Rule for solving.
         * - Dynamically builds and solves master and subproblems.
         * - Supports the addition of feasibility and optimality cuts to refine solutions.
         * 
         * Flow:
         * - User input → Build master problem → Solve master problem → Solve subproblems → Add cuts → Iterate until convergence.
         * 
         * Utility Functions:
         * - **BuildMasterProblemData**: Constructs the master problem with user-selected constraints.
         * - **SolveSubproblemScenario**: Solves subproblems for given yield scenarios.
         * - **AddOptimalityCut**: Adds cuts to refine second-stage approximation.
         * - **SolveWithDualAndPrimal**: Determines the appropriate simplex method and solves.
         * - **ComputeScenarioMultipliers**: Generates yield multipliers for scenario-based analysis.
         * 
         * Additional Updates:
         * - Handles degenerate cases with perturbation.
         * - Tracks pivot history to prevent cycling in simplex.
         * - Modular functions ensure separation between problem phases.
         */

        // Original farmer problem variables in order: [x1, x2, x3, y1, y2, w1, w2, w3, w4, theta]
        // We add 'theta' for L-shaped method

        // Master constraints: user picks which constraints go in master.
        // The rest + scenario yield logic goes to subproblem.

        // Store the original problem constraints in a data structure:
        static double[,] originalLHS = new double[5, 10];  // 5 constraints, 10 vars
        static double[] originalRHS = new double[5];
        static string[] consLabels = { "C1", "C2", "C3", "C4", "C5" };
        static Dictionary<string, int> consIndexMap = new Dictionary<string, int> { { "C1", 0 }, { "C2", 1 }, { "C3", 2 }, { "C4", 3 }, { "C5", 4 } };
        static List<string> allVars = new List<string> { "x_1", "x_2", "x_3", "y_1", "y_2", "w_1", "w_2", "w_3", "w_4", "theta" };

        // Objective coefficients (original):
        static Dictionary<string, double> varObj = new Dictionary<string, double>{
            {"x_1",-150},{"x_2",-230},{"x_3",-260},
            {"y_1",238},{"y_2",210},
            {"w_1",-170},{"w_2",-150},{"w_3",-36},{"w_4",-10},
            {"theta",0.0}
        };

        // Holds each constraint as (Label, LHS row, RHS).
        // LHS is an array of length 10 since you have 10 variables {x1, x2, x3, y1, y2, w1, w2, w3, w4, theta}.
        static List<(string Label, double[] LHS, double RHS)> allConstraintData
            = new List<(string Label, double[] LHS, double RHS)>();
        // Global storage for all constraints
        static double[,] AllLHS; // 2D array for LHS coefficients
        static double[] AllRHS;  // 1D array for RHS values



        static void Main(string[] args)
        {
            originalLHS[0, 0] = 1; originalLHS[0, 1] = 1; originalLHS[0, 2] = 1; originalRHS[0] = 500;
            // C2: 2.5x1 + y1 - w1≥200 => -2.5x1 -y1 +w1 ≤ -200
            originalLHS[1, 0] = -2.5; originalLHS[1, 3] = -1; originalLHS[1, 5] = 1; originalRHS[1] = -200;
            // C3: 3x2 + y2 - w2≥240 => -3x2 -y2 + w2 ≤ -240
            originalLHS[2, 1] = -3; originalLHS[2, 4] = -1; originalLHS[2, 6] = 1; originalRHS[2] = -240;
            // C4: w3+w4 ≤20x3 => w3+w4 -20x3≤0
            originalLHS[3, 7] = 1; originalLHS[3, 8] = 1; originalLHS[3, 2] = -20; originalRHS[3] = 0;
            // C5: w3 ≤6000
            originalLHS[4, 7] = 1; originalRHS[4] = 6000;

            // (1) Let user define the number of scenarios
            Console.WriteLine("Enter the number of scenarios:");
            int numScenarios;
            while (!int.TryParse(Console.ReadLine(), out numScenarios) || numScenarios <= 0)
            {
                Console.WriteLine("Invalid input. Please enter a positive integer for the number of scenarios:");
            }

            // (2) Dynamically generate constraints for scenarios
            List<string> allConstraints = GenerateScenarioConstraints(numScenarios);
            // Generate all variables (x_1,x_2,x_3 for first-stage, y_1_s,y_2_s,w_1_s,... for second-stage)
            List<string> allVariables = GenerateAllVariables(numScenarios);
            // Display available constraints for user selection
            Console.WriteLine("\nAvailable constraints:");
            for (int i = 0; i < allConstraints.Count; i++)
            {
                Console.WriteLine($"{allConstraints[i]}");
            }

            // Example: Assuming allConstraintData is already populated
            int numConstraints = allConstraintData.Count;
            int numVars = allConstraintData[0].LHS.Length;

            // Initialize global arrays
            AllLHS = new double[numConstraints, numVars];
            AllRHS = new double[numConstraints];

            // Populate AllLHS and AllRHS
            for (int i = 0; i < numConstraints; i++)
            {
                var (label, lhs, rhs) = allConstraintData[i];

                // Store LHS coefficients in AllLHS
                for (int j = 0; j < numVars; j++)
                {
                    AllLHS[i, j] = lhs[j];
                }

                // Store RHS in AllRHS
                AllRHS[i] = rhs;
            }
            static List<string> GenerateAllVariables(int numScenarios)
            {
                List<string> vars = new List<string>();
                vars.Add("theta");
                vars.Add("x_1");
                vars.Add("x_2");
                vars.Add("x_3");

                for (int s = 1; s <= numScenarios; s++)
                {
                    vars.Add($"y_1_{s}");
                    vars.Add($"y_2_{s}");
                }
                for (int s = 1; s <= numScenarios; s++)
                {
                    vars.Add($"w_1_{s}");
                    vars.Add($"w_2_{s}");
                    vars.Add($"w_3_{s}");
                    vars.Add($"w_4_{s}");
                }
                return vars;
            }

            // (3) Allow user to select constraints
            Console.WriteLine("\nEnter the constraints to include (e.g., C1,C2_S1,C3_S2):");
            string inputConstraints = Console.ReadLine();
            string[] selectedConstraintLabels = inputConstraints.Split(',')
                .Select(s => s.Trim().ToUpper())
                .Where(s => s.StartsWith("C"))
                .ToArray();
            // --------------------------------------------------------------
            // Force "C5" if user picks any constraint referencing "W_3"
            // --------------------------------------------------------------
            bool referencesW3 = selectedConstraintLabels.Any(lbl => lbl.Contains("C4_"));
            if (referencesW3)
            {
                // If "C5" is not already included, forcibly add it
                if (!selectedConstraintLabels.Contains("C5_"))
                {
                    selectedConstraintLabels = selectedConstraintLabels
                        .Concat(new[] { "C5" })
                        .Distinct() // in case user typed duplicates
                        .ToArray();
                    Console.WriteLine(
                        "\nDetected constraints referencing w_3; automatically adding C5 (w_3 ≤ 6000)!"
                    );
                }
            }


            HashSet<int> selectedConstraintIndices = new HashSet<int>();
            foreach (var label in selectedConstraintLabels)
            {
                if (consIndexMap.ContainsKey(label))
                {
                    selectedConstraintIndices.Add(consIndexMap[label]);
                }
                else
                {
                    Console.WriteLine($"Warning: Constraint {label} not found in consIndexMap.");
                }
            }

            // Validate selection
            if (selectedConstraintIndices.Count == 0)
            {
                Console.WriteLine("No valid constraints selected. Exiting program.");
                return;
            }

            // Build master problem with selected constraints
            (double[,] masterLHS, double[] masterRHS, double[] masterCoeffs, List<string> masterVarNames) =
                BuildMasterProblemData(selectedConstraintIndices);

            // master dimension:
            int masterCons = masterRHS.Length;
            int masterVarsCount = masterCoeffs.Length;
            int subCons = AllRHS.Length - masterRHS.Length;
            int subVarsCount = numScenarios * 6 + 3 - masterVarsCount + 2;

            // Subproblem constraints = the constraints not chosen
            List<int> subConsIndices = new List<int>();
            for (int i = 0; i < AllRHS.Length; i++)
                if (!selectedConstraintIndices.Contains(i)) subConsIndices.Add(i);

            // Subproblem variables = the original minus the master var set (and "theta")
            // We find the name "theta" or slack not in subproblem
            HashSet<string> masterVarSet = new HashSet<string>(masterVarNames);
            List<int> subVarIndices = new List<int>();
            for (int idx = 0; idx < allVariables.Count; idx++)
            {
                string v = allVariables[idx];
                if (!masterVarSet.Contains(v) && !v.StartsWith("slack_") && v != "theta")
                {
                    subVarIndices.Add(idx);
                }
            }

            // Check complete recourse: If user picks any constraint that has y_ or w_ => not complete recourse
            bool isCompleteRecourse = true;
            foreach (int ci in selectedConstraintIndices)
            {
                // Check if that constraint includes any y_ or w_
                for (int j = 0; j < 10; j++)
                {
                    if (Math.Abs(originalLHS[ci, j]) > 1e-15)
                    {
                        string varName = allVariables[j];
                        if (varName.StartsWith("y_") || varName.StartsWith("w_"))
                            isCompleteRecourse = false;
                    }
                }
            }
            // Compute scenario multipliers
            double[] scenarioMultipliers = ComputeScenarioMultipliers(numScenarios);

            // L-shaped decomposition setup
            double UB = double.PositiveInfinity;
            double LB = double.NegativeInfinity;
            int maxIter = 500;
            bool optimalFound = false;

            double[] finalMasterSolution = null; // Store final solution
            double finalObjectiveValue = 0.0;    // Store final objective value

            for (int iteration = 1; iteration <= maxIter; iteration++)
            {
                Console.WriteLine($"\n--- L-Shaped Iteration {iteration} ---");

                // 1) Solve the master problem
                var masterResult = SolveWithDualAndPrimal(masterLHS, masterRHS, masterCoeffs);
                if (!masterResult.feasible)
                {
                    Console.WriteLine("Master problem infeasible. Stopping iteration.");
                    break;
                }

                // Extract solution and objective from master
                double[] masterSolution = masterResult.solution;
                double masterObjective = masterResult.objVal;

                // 2) Update upper bound (best known feasible solution)
                UB = Math.Min(UB, masterObjective);

                // 3) Solve subproblems for each scenario
                double totalScenarioCost = 100.0; // some base cost or expression
                bool infeasibleScenario = false;

                for (int s = 0; s < numScenarios; s++)
                {
                    // Solve subproblem for scenario s
                    var (subResult, y, omega, pi, feasibility, subLHS, subRHS, reducedLHS) =
                        SolveSubproblemScenario(
                            masterSolution,
                            masterVarNames,
                            selectedConstraintIndices.ToList(),
                            new List<int>(),
                            scenarioMultipliers[s] * 2.5,
                            scenarioMultipliers[s] * 3.0,
                            scenarioMultipliers[s] * 20.0,
                            isCompleteRecourse,
                            s,
                            numScenarios,
                            allVariables
                        );

                    // -------------------------
                    //    CHECK FEASIBILITY
                    // -------------------------
                    // The CheckFeasibility function can do your dual calculations or
                    // anything else that determines if the subproblem is infeasible
                    // or needs a feasibility cut. For demonstration, let's assume it
                    // returns a boolean plus the cut arrays, e.g.:
                    // (isFeasible, D_r, d_r)
                    var (isFeasible, D_r_Matrix, d_r) = CheckFeasibility(
                        masterSolution,
                        masterVarNames,
                        selectedConstraintIndices.ToList(),
                        new List<int>(),
                        scenarioMultipliers,
                        true,       // isCompleteRecourse
                        numScenarios,
                        masterVarNames.Count,
                        new HashSet<string>(masterVarNames),
                        allVariables
                    );

                    if (!isFeasible)
                    {
                        Console.WriteLine($"Scenario {s} infeasible. Adding feasibility cut.");

                        // -------------------------
                        //   ADD FEASIBILITY CUT
                        // -------------------------
                        // The AddFeasibilityCut function takes the existing master problem
                        // matrices and appends the feasibility cut:
                        (masterLHS, masterRHS, masterCoeffs) = AddFeasibilityCut(
                            masterLHS,
                            masterRHS,
                            masterCoeffs,
                            D_r_Matrix,   // from CheckFeasibility
                            d_r,          // from CheckFeasibility
                            masterVarNames.Count,
                            new HashSet<string>(masterVarNames)
                        );

                        // Mark that we found an infeasible scenario.
                        infeasibleScenario = true;
                        break;
                    }
                    else
                    {
                        // If scenario is feasible, keep track of scenario cost
                        // "omega" is typically subproblem's cost or dual result
                        totalScenarioCost += omega / numScenarios;
                    }
                }

                // If at least one scenario was infeasible, go to the next iteration
                if (infeasibleScenario)
                    continue;

                // 4) Update lower bound
                double thetaValue = masterSolution[masterVarNames.IndexOf("theta")];
                LB = Math.Max(LB, thetaValue);

                // 5) Check for convergence
                if (Math.Abs(UB - LB) < 1e-6)
                {
                    Console.WriteLine("Converged: UB ≈ LB");
                    optimalFound = true;
                    finalMasterSolution = masterSolution;
                    finalObjectiveValue = masterObjective;
                    break;
                }

                // 6) Add optimality cut if needed
                if (thetaValue < totalScenarioCost - 1e-6)
                {
                    Console.WriteLine($"Adding optimality cut: theta ≥ {totalScenarioCost:F2}");

                    (masterLHS, masterRHS, masterCoeffs) = AddOptimalityCut(
                        ref masterLHS,
                        ref masterRHS,
                        ref masterCoeffs,
                        masterVarNames.IndexOf("theta"),
                        totalScenarioCost,
                        masterVarNames.Count,
                        masterSolution,
                        new HashSet<string>(masterVarNames),
                        selectedConstraintIndices.ToList(),
                        new List<int>(),
                        numScenarios,
                        true,
                        scenarioMultipliers,
                        masterVarNames,
                        allVariables
                    );
                }
            }


            // Final solution or termination
            if (optimalFound)
            {
                Console.WriteLine("\nOptimal Solution Found:");
                for (int i = 0; i < masterVarNames.Count; i++)
                {
                    Console.WriteLine($"{masterVarNames[i]} = {finalMasterSolution[i]:F2}");
                }
                Console.WriteLine($"Objective Value: {finalObjectiveValue:F2}");
            }
            else
            {
                Console.WriteLine("Maximum iterations reached or problem did not converge.");
            }
        }


        // Helper: Generate constraints dynamically based on scenarios
        static List<string> GenerateScenarioConstraints(int numScenarios)
        {
            List<string> constraints = new List<string>();
            consIndexMap.Clear(); // Clear existing mapping

            // (A) Add the base constraint (C1: total land constraint) as before:
            // -------------------------------------------------------------
            string baseConstraint = $"{consLabels[0]}: ";
            double[] baseLHSRow = new double[10]; // 10 variables
            for (int j = 0; j < allVars.Count; j++)
            {
                double coeff = originalLHS[0, j];
                baseLHSRow[j] = coeff;
                if (Math.Abs(coeff) > 1e-9)
                {
                    baseConstraint += $"{coeff:F1}*{allVars[j]} + ";
                }
            }
            baseConstraint = baseConstraint.TrimEnd(' ', '+') + $" ≤ {originalRHS[0]:F1}";
            constraints.Add(baseConstraint);

            // Save to allConstraintData 
            allConstraintData.Add((
                Label: consLabels[0],
                LHS: baseLHSRow,
                RHS: originalRHS[0]
            ));

            consIndexMap["C1"] = 0; // Add to mapping

            // (B) Generate scenario-specific constraints for C2–C5
            // -------------------------------------------------------------
            int constraintIndex = 1; // We already used 0 for C1
            for (int scenario = 1; scenario <= numScenarios; scenario++)
            {
                double scenarioMultiplier = 0.8 + (1.2 - 0.8) * (scenario - 1) / (numScenarios - 1);

                // For each of the rows i in the original constraints (C2..C5 => i=1..4)
                // For each original constraint row i in {C2..C5}
                for (int i = 1; i < originalLHS.GetLength(0); i++)
                {
                    string uniqueLabel = $"{consLabels[i]}_S{scenario}";
                    string constraintText = $"{uniqueLabel}: ";

                    // numeric LHS row of length 10
                    double[] numericLHSRow = new double[10];

                    for (int j = 0; j < 10; j++)
                    {
                        double coeff = originalLHS[i, j];
                        // Apply scenario multiplier if needed:
                        if (i == 1 && j == 0) coeff *= scenarioMultiplier; // C2, x1
                        if (i == 2 && j == 1) coeff *= scenarioMultiplier; // C3, x2
                        if (i == 3 && j == 2) coeff *= scenarioMultiplier; // C4, x3

                        numericLHSRow[j] = coeff;

                        if (Math.Abs(coeff) > 1e-9)
                        {
                            // Renaming logic:
                            //  - If the variable is x_1, x_2, or x_3 => keep original name
                            //  - Otherwise => variableName_S{scenario}
                            string baseVarName = allVars[j];

                            // If it's not x_1, x_2, or x_3, rename
                            if (!baseVarName.StartsWith("x_"))
                            {
                                // e.g. if baseVarName="w_3", scenario=2 => "w_3_2"
                                baseVarName += $"_{scenario}";
                            }

                            // Now build the text term
                            constraintText += $"{coeff:F1}*{baseVarName} + ";
                        }
                    }

                    // Build final text “ ≤ originalRHS[i]”
                    double cRHS = originalRHS[i];
                    constraintText = constraintText.TrimEnd(' ', '+') + $" ≤ {cRHS:F1}";

                    // (B1) Save the textual constraint
                    constraints.Add(constraintText);

                    // (B2) Save the numeric form
                    allConstraintData.Add((
                        Label: uniqueLabel,
                        LHS: numericLHSRow,
                        RHS: cRHS
                    ));

                    // Update the index map
                    consIndexMap[uniqueLabel] = allConstraintData.Count - 1; // or `constraintIndex`
                    constraintIndex++;
                }
            }
            return constraints;
        }

        public static (bool isFeasible, double[,] D_r_Matrics, double[] d_r)
    CheckFeasibility(
    double[] mSol,
    List<string> masterVarNames,
    List<int> subConsIndices,
    List<int> subVarIndices,
    double[] scenarioMultipliers,
    bool isCompleteRecourse,
    int numScenarios,
    int totalVars,
    HashSet<string> masterVarSet,
    List<string> allVariables
)
        {
            // Determine the number of subproblem variables (e.g., y1, y2, w1, etc.)
            int subVarCount = 6 * numScenarios + 3 - masterVarSet.Count + 2;

            // Initialize D_r_Matrics and d_r with appropriate dimensions
            // Rows = number of subproblem variables, Columns = number of master variables
            double[,] D_r_Matrics = new double[numScenarios, masterVarSet.Count - 2];
            double[] d_r = new double[subVarCount];

            // Initialize as feasible
            bool isFeasible = true;

            // Loop through each scenario
            for (int s = 0; s < numScenarios; s++)
            {
                double Multi1 = scenarioMultipliers[s] * 2.5;
                double Multi2 = scenarioMultipliers[s] * 3.0;
                double Multi3 = scenarioMultipliers[s] * 20.0;

                var (result, y, omega, pi, feasibility, subLHS, subRHS, reducedLHS) =
                    SolveSubproblemScenario(
                        mSol,
                        masterVarNames,
                        subConsIndices,
                        subVarIndices,
                        Multi1,
                        Multi2,
                        Multi3,
                        isCompleteRecourse,
                        s,
                        numScenarios,
                        allVariables
                    );

                if (!feasibility)
                {
                    // If any scenario is infeasible, set flag and break
                    isFeasible = false;
                    Console.WriteLine($"Scenario {s + 1} is infeasible.");
                    break;
                }

                // Extract dual variables (pi) from the subproblem
                // Assuming pi corresponds to the dual variables for the constraints in subproblem
                // and their order matches the master variables
                for (int j = 0; j < pi.Length; j++)
                {
                    D_r_Matrics[s, j] = pi[j];
                }

                // Calculate d_r[s] = pi^T * b
                // Here, 'b' is the RHS of the subproblem constraints
                // Ensure that 'pi' and 'subRHS' have compatible dimensions
                if (pi.Length != subRHS.Length)
                {
                    throw new InvalidOperationException($"Dual variables length ({pi.Length}) does not match subRHS length ({subRHS.Length}) for scenario {s + 1}.");
                }

                double piT_b = 0.0;
                for (int i = 0; i < pi.Length; i++)
                {
                    piT_b += pi[i] * subRHS[i];
                }

                d_r[s] = piT_b;
            }

            return (isFeasible, D_r_Matrics, d_r);
        }



        public static (double[,] updatedLHS, double[] updatedRHS, double[] updatedObj)
AddFeasibilityCut(
    double[,] masterLHS,
    double[] masterRHS,
    double[] masterCoeffs,
    double[,] D_r_Matrics,
    double[] d_r,
    int totalVars,
    HashSet<string> masterVarSet
)
        {
            // Number of existing constraints and variables
            int oldCons = masterRHS.Length;
            int oldVars = masterLHS.GetLength(1);

            // Number of new feasibility cuts
            int newCuts = D_r_Matrics.GetLength(0); // One cut per scenario

            // Create updated LHS and RHS with new cuts
            double[,] updatedLHS = new double[oldCons + newCuts, oldVars];
            double[] updatedRHS = new double[oldCons + newCuts];
            double[] updatedObj = new double[masterCoeffs.Length]; // Objective remains the same initially

            // Copy existing constraints
            for (int i = 0; i < oldCons; i++)
            {
                for (int j = 0; j < oldVars; j++)
                {
                    updatedLHS[i, j] = masterLHS[i, j];
                }
                updatedRHS[i] = masterRHS[i];
            }

            // Append new feasibility cuts
            for (int k = 0; k < newCuts; k++)
            {
                for (int j = 0; j < oldVars; j++)
                {
                    updatedLHS[oldCons + k, j] = D_r_Matrics[k, j];
                }
                updatedRHS[oldCons + k] = d_r[k];
            }

            // Copy existing objective coefficients
            Array.Copy(masterCoeffs, updatedObj, masterCoeffs.Length);

            // Optionally, update the objective if necessary (e.g., including theta)
            // This depends on how theta is incorporated in the master problem

            return (updatedLHS, updatedRHS, updatedObj);
        }


        static (double[,] lhs, double[] rhs, double[] obj) AddOptimalityCut(ref double[,] masterLHS
            , ref double[] masterRHS, ref double[] masterCoeffs, int idxTheta, double w_v, int totalVars, double[] mSol, HashSet<string> masterVarSet
            , List<int> subConsIndices, List<int> subVarIndices, int numScenarios, bool isCompleteRecourse, double[] scenarioMultipliers, List<string> masterVarNames, List<string> allVariables)
        {

            // Step 3: Calculate E_s and e_s
            double[] E_s_vector = new double[totalVars - 2];
            double e_s = 0.0;
            int variableCount = totalVars - 2;
            double[,] piT = new double[variableCount, numScenarios];
            double P_k = 1.0 / numScenarios;
            for (int s = 0; s < numScenarios; s++)
            {
                double Multi1 = scenarioMultipliers[s] * 2.5;
                double Multi2 = scenarioMultipliers[s] * 3;
                double Multi3 = scenarioMultipliers[s] * 20;


                var (result, y, o, pi, feasibility, subLHS, subRHS, reducedLHS) = SolveSubproblemScenario(mSol, masterVarNames, subConsIndices, subVarIndices,
                 Multi1, Multi2, Multi3, isCompleteRecourse, s, numScenarios, allVariables);

                // Add row: theta≥ w_v => row: theta - w_v≥0
                int oldRows = masterLHS.GetLength(0);
                int oldCols = masterLHS.GetLength(1);

                double[] pi_k = pi; //pi is dual solution 

                int m = subRHS.Length;

                for (int j = 0; j < reducedLHS.GetLength(1); j++)
                {
                    double sum = 0.0;
                    for (int i = 0; i < m; i++)
                    {
                        sum += pi_k[i] * reducedLHS[i, j];
                    }
                    piT[j, s] = sum;
                }

                double pi_h = 0.0;
                for (int i = 0; i < m; i++)
                {
                    pi_h += pi_k[i] * subRHS[i];
                }

                e_s += P_k * pi_h;
            }


            for (int j = 0; j < variableCount; j++)
            {
                for (int s = 0; s <= numScenarios - 1; s++)
                {
                    E_s_vector[j] += P_k * piT[j, s];
                }
            }

            // Print e_s and E_s_vector
            Console.WriteLine("\n=== e_s and E_s_vector ===");
            Console.WriteLine($"e_s: {e_s:F2}");
            Console.Write("E_s_vector: ");
            foreach (var item in E_s_vector)
            {
                Console.Write($"{item:F2} ");
            }
            Console.WriteLine();

            double optimality = e_s;
            int n = totalVars - 2;
            for (int j = 0; j < n; j++)
            {
                optimality -= E_s_vector[j] * mSol[j];
            }

            int oldCons = masterRHS.Length;
            double[,] updatedLHS = new double[oldCons + 1, totalVars + oldCons];
            double[] updatedRHS = new double[oldCons + 1];
            for (int i = 0; i < oldCons; i++)
            {
                for (int j = 0; j < totalVars; j++)
                    updatedLHS[i, j] = masterLHS[i, j];
                updatedRHS[i] = masterRHS[i];
            }
            for (int j = 0; j < totalVars; j++)
            {
                updatedLHS[oldCons, j] = 0.0;
            }
            updatedRHS[oldCons] = e_s;
            Console.WriteLine("\n=== Updated RHS after setting e_s ===");
            PrintArray(updatedRHS);

            for (int k = 0; k < idxTheta; k++)
            {
                updatedLHS[oldCons, k] = E_s_vector[k];

            }

            // Add slack variables to the LHS of the new constraint row
            for (int k = 0; k < oldCons; k++)
            {
                updatedLHS[k + 1, totalVars + k] = 1.0;  // Slack variable coefficients are 1
            }


            for (int i = 1; i < updatedLHS.GetLength(0); i++)
            {
                updatedLHS[i, idxTheta] = 1.0;
            }

            Console.WriteLine("\n=== Updated LHS after copying last column to 1.0 ===\"");
            PrintMatrix(updatedLHS, oldCons + 1, totalVars + oldCons);
            Console.WriteLine($"\nAdding Optimality Cut: theta[{idxTheta}] ≥ {e_s:F2}");

            int updatedObjLength = updatedLHS.GetLength(1);
            double[] updatedObj = new double[updatedObjLength];

            Array.Copy(masterCoeffs, updatedObj, masterCoeffs.Length);

            for (int k = masterCoeffs.Length; k < updatedObj.Length; k++)
            {
                updatedObj[k] = 0.0;
            }

            updatedObj[idxTheta] = -1.0;

            masterLHS = updatedLHS; masterRHS = updatedRHS; masterCoeffs = updatedObj;


            return (masterLHS, masterRHS, masterCoeffs);
        }


        // Build Master Problem Data:
        static (double[,] lhs, double[] rhs, double[] coeffs, List<string> varNames)
        BuildMasterProblemData(HashSet<int> chosenConsIdx)
        {
            // 1) Validate chosen constraints against the allConstraintData list
            if (chosenConsIdx.Any(ci => ci < 0 || ci >= allConstraintData.Count))
            {
                Console.WriteLine("Invalid constraint index in chosenConsIdx.");
                return default;
            }


            // 2) Identify the index of 'theta'
            int thetaIndex = allVars.IndexOf("theta");
            if (thetaIndex == -1)
            {
                Console.WriteLine("'theta' variable not found in allVars.");
                return default;
            }

            // 3) Determine which variables appear in the chosen constraints
            //    We'll also ensure 'theta' is included.
            HashSet<int> masterVarIndices = new HashSet<int>();
            foreach (int ci in chosenConsIdx)
            {
                // For each selected constraint
                var (constraintLabel, lhsArray, rhsVal) = allConstraintData[ci];
                for (int j = 0; j < lhsArray.Length; j++)
                {
                    if (Math.Abs(lhsArray[j]) > 1e-15)
                    {
                        masterVarIndices.Add(j);
                    }
                }
            }

            // Ensure we include 'theta'
            masterVarIndices.Add(thetaIndex);

            // 4) Sort the variable indices for consistent ordering
            List<int> mList = masterVarIndices.ToList();
            mList.Sort();

            // 5) The master problem has one slack variable for each chosen constraint
            int masterCons = chosenConsIdx.Count;
            int slack = masterCons;
            int totalVars = mList.Count + slack;

            // 6) Prepare arrays for the master LHS, RHS, and variable names
            double[,] masterLHS = new double[masterCons, totalVars];
            double[] masterRHS = new double[masterCons];
            List<string> masterVarNames = new List<string>();

            // 7) Structural variable names
            foreach (int idx in mList)
            {
                masterVarNames.Add(allVars[idx]);
            }
            // 8) Slack variable names
            for (int i = 0; i < masterCons; i++)
            {
                masterVarNames.Add($"slack_{i + 1}");
            }

            // 9) Fill master LHS and RHS
            int rowIndex = 0;
            foreach (int ci in chosenConsIdx)
            {
                var (lbl, rowLHSArray, rowRHSval) = allConstraintData[ci];

                // For each structural variable in mList
                for (int j = 0; j < mList.Count; j++)
                {
                    int globalVarIdx = mList[j];
                    masterLHS[rowIndex, j] = rowLHSArray[globalVarIdx];
                }

                // Slack variable for this row
                masterLHS[rowIndex, mList.Count + rowIndex] = 1.0;

                // RHS
                masterRHS[rowIndex] = rowRHSval;
                rowIndex++;
            }

            // 10) Build the objective coefficients
            double[] masterCoeffs = new double[totalVars];
            for (int j = 0; j < mList.Count; j++)
            {
                string v = allVars[mList[j]];
                masterCoeffs[j] = varObj[v];
            }
            // Slack variables have zero cost

            // Return the data needed for simplex
            return (masterLHS, masterRHS, masterCoeffs, masterVarNames);
        }



        static void PrintMatrix(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    Console.Write($"{matrix[i, j]:F2} "); // Format to 2 decimal places
                }
                Console.WriteLine(); // New line after each row
            }
        }

        static void PrintArray(double[] array, string label)
        {
            Console.WriteLine($"{label}: {string.Join(", ", array.Select(x => x.ToString("F2")))}");
        }


        // Build initial tableau for primal simplex
        static (double[,] lhs, double[] rhs, double[] obj) BuildMasterTableau(double[,] masterLHS, double[] masterRHS, double[] masterCoeffs)
        {
            int m = masterRHS.Length;  // #constraints
            int tvars = masterCoeffs.Length; //#vars
                                             // We'll transform Ax ≤ b into a tableau form
                                             // The last m columns will be the slack variables, presumably an identity basis.
            double[,] lhs = new double[m, tvars];
            double[] rhs = new double[m];
            double[] obj = new double[tvars];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < tvars; j++)
                {
                    lhs[i, j] = masterLHS[i, j];
                }
                rhs[i] = masterRHS[i];
            }
            for (int j = 0; j < tvars; j++)
                obj[j] = masterCoeffs[j];
            return (lhs, rhs, obj);
        }

        // Solve subproblem scenario with internal approach
        // For demonstration, we do a direct approach:
        // subproblem constraints = subConsIndices with yield variations substituted
        // subproblem objective from varObj for subVarIndices
        public struct SubproblemResult { public bool feasible; public double cost; }
        static (SubproblemResult result, double[] y, double omega, double[] pi, bool feasibility, double[,] subLHS, double[] subRHS, double[,] reducedLHS)
            SolveSubproblemScenario(double[] masterSolVal, List<string> masterVarNames,
            List<int> subConsIndices, List<int> subVarIndices, double wheatYield, double cornYield, double sugarYield,
            bool isCompleteRecourse, int s, int numScenarios, List<string> allVariables)
        {
            // Build subproblem tableau
            // subproblem constraints => subConsIndices
            // subproblem variables => subVarIndices
            // objective => sum(varObj[vName]*subVar)
            // For demonstration: we'll treat constraints as ≤ form, add slack => BFS
            // Then do primal simplex.

            // Build LHS,RHS for subproblem
            int subM = subConsIndices.Count;
            int subVarCount = subVarIndices.Count;
            double[,] subLHS = new double[subM, subVarCount + subM];
            double[] subRHS = new double[subM];
            double[] subObj = new double[subVarCount + subM];

            for (int j = 0; j < subVarCount; j++)
            {
                string v = allVariables[subVarIndices[j]];
                subObj[j] = varObj[v];
            }
            for (int i = 0; i < subM; i++) subObj[subVarCount + i] = 0.0;

            for (int i = 0; i < subM; i++)
            {
                int ci = subConsIndices[i];
                double[] cRow = new double[10];
                for (int kk = 0; kk < 10; kk++)
                {
                    cRow[kk] = AllLHS[ci, kk];
                }
                double cRHS = AllRHS[ci];

                // Adjust yields:
                if (ci == 1)
                {
                    cRow[0] = -wheatYield; //C2
                }
                if (ci == 2)
                {
                    cRow[1] = -cornYield; //C3
                }
                if (ci == 3)
                {
                    cRow[2] = -sugarYield; //C4
                }

                // Sub out master variables
                double hVal = cRHS;

                for (int mv = 0; mv < masterVarNames.Count; mv++)
                {
                    string varM = masterVarNames[mv];
                    int globalIdx = allVars.IndexOf(varM);
                    if (globalIdx >= 0 && Math.Abs(cRow[globalIdx]) > 1e-15)
                    {
                        hVal -= cRow[globalIdx] * masterSolVal[mv];
                        cRow[globalIdx] = 0.0;
                    }
                }

                // Fill subLHS row
                for (int j = 0; j < subVarCount; j++)
                {
                    int globalID = subVarIndices[j];
                    subLHS[i, j] = cRow[globalID];
                }
                subLHS[i, subVarCount + i] = 1.0; // slack
                subRHS[i] = hVal;
            }


            // Get the indices of variables from subLHS that match the master variables
            List<int> selectedVarIndices = new List<int>();
            int masterLength = allVariables.Count - subVarCount - 1;

            foreach (string masterVar in masterVarNames)
            {
                int globalIndex = allVariables.IndexOf(masterVar);

                if (masterVar == "theta")
                {
                    break;
                }

                if (globalIndex >= 0)
                {
                    selectedVarIndices.Add(globalIndex);
                }

            }

            // Initialize a smaller matrix for storing the reduced LHS
            int numRows = subLHS.GetLength(0);
            int numCols = selectedVarIndices.Count;
            double[,] reducedLHS = new double[numRows, numCols];

            // Fill the reduced LHS with the relevant coefficients
            for (int i = 0; i < numRows; i++)
            {
                for (int j = 0; j < numCols; j++)
                {
                    reducedLHS[i, j] = 0;
                }
            }

            for (int i = 0; i < numRows; i++)
            {
                for (int j = 0; j < numCols; j++) // 遍歷 selectedVarIndices
                {

                    int index = selectedVarIndices[j]; // 對應的變數索引
                    reducedLHS[i, j] = AllLHS[i + 1, index];

                    if (index == 0 && i == 0) // 變數 x1 的倍數/一定是第一個
                    {
                        reducedLHS[index, index] = -wheatYield;
                    }
                    else if (index == 1 && selectedVarIndices[0] == 0 && i == 1) //變數為x2時前面有x1(x2就會在[1,1])
                    {
                        reducedLHS[index, index] = (-cornYield);
                    }
                    else if (index == 1 && selectedVarIndices[0] != 0 && i == 0)  //變數為x2時前面沒有x1(x2有倍數的地方就會在[0,1])
                    {
                        reducedLHS[j, index] = (-cornYield);
                    }
                    else if (index == 2 && selectedVarIndices[0] == 0 && selectedVarIndices[1] == 1 && i == 2) // 變數 x3 的倍數前面有x1和x2
                    {
                        reducedLHS[index, index] = -sugarYield;
                    }
                    else if (index == 2 && selectedVarIndices[0] == 0 && selectedVarIndices[1] != 1 && i == 1)// 變數 x3 的前面有x1沒x2
                    {
                        reducedLHS[j, index] = -sugarYield;
                    }
                    else if (index == 2 && selectedVarIndices[0] == 1 && i == 1)// 變數 x3 的前面有x2沒x1
                    {
                        reducedLHS[j, index] = -sugarYield;
                    }
                    else if (index == 2 && selectedVarIndices[0] == 2 && i == 0)  // 變數 x3 前面甚麼都沒有
                    {
                        reducedLHS[j, index] = -sugarYield;
                    }
                    else // 其他變數
                    {
                        reducedLHS[i, j] = AllLHS[i + 1, index];
                    }
                }
            }

            // Now primal simplex
            var spTableau = BuildSubproblemTableau(subLHS, subRHS, subObj, subM, subVarCount);

            var spRes = PrimalSimplexSolveFull(spTableau.lhs, spTableau.rhs, spTableau.obj);

            var y = spRes.solution;
            var omega = spRes.objVal;
            var feasibility = spRes.feasible;

            var spResDual = DualSimplexSolveFull(subLHS, spTableau.rhs, spTableau.obj);
            var pi = spRes.solution;

            SubproblemResult result = new SubproblemResult();
            result.feasible = spRes.feasible;
            if (!spRes.feasible)
            {
                result.cost = Double.PositiveInfinity;

                return (result, y, omega, pi, feasibility, subLHS, subRHS, reducedLHS);
            }
            result.cost = spRes.objVal;

            return (result, y, omega, pi, feasibility, subLHS, subRHS, reducedLHS);
        }

        static (double[,] lhs, double[] rhs, double[] obj) BuildSubproblemTableau(double[,] subLHS, double[] subRHS, double[] subObj, int subM, int subVarCount)
        {
            // subproblem total vars = subVarCount + subM
            // We'll transform them into a tableau
            int totalVars = subVarCount + subM;
            double[,] lhs = new double[subM, totalVars];
            double[] rhs = new double[subM];
            double[] obj = new double[totalVars];
            for (int i = 0; i < subM; i++)
            {
                for (int j = 0; j < totalVars; j++)
                    lhs[i, j] = subLHS[i, j];
                rhs[i] = subRHS[i];
            }
            for (int j = 0; j < totalVars; j++)
                obj[j] = subObj[j];

            return (lhs, rhs, obj);
        }

        // Primal Simplex result
        public struct SimplexResult
        {
            public bool feasible;
            public double[] solution;
            public double objVal;
        }

        // Full primal simplex solver 
        static SimplexResult DualSimplexSolveFull(double[,] lhs, double[] rhs, double[] obj)
        {
            int m = rhs.Length; // Number of constraints
            int n = obj.Length; // Number of variables
            int iteration = 0;

            // Make local copies of inputs
            double[,] A = (double[,])lhs.Clone();
            double[] B = (double[])rhs.Clone();
            double[] c = (double[])obj.Clone();

            int[] basicVar = new int[m];
            for (int i = 0; i < m; i++)
                basicVar[i] = n - m + i; // Initialize slack variables as basic variables

            const double EPS = 1e-9; // Precision threshold
            HashSet<string> pastStates = new HashSet<string>();
            HashSet<(int, int)> pivotHistory = new HashSet<(int, int)>(); // Track pivots

            while (true)
            {
                iteration++;

                if (iteration >= 100)
                {
                    var blandResult = SolveUsingBlandsRule(lhs, rhs, obj);
                    return blandResult;
                }

                // Step 1: Check dual feasibility (find the leaving variable)
                int leaving = -1;
                double minB = double.MaxValue;
                for (int i = 0; i < m; i++)
                {
                    if (B[i] < -EPS && (B[i] < minB || (Math.Abs(B[i] - minB) < EPS && i < leaving)))
                    {
                        minB = B[i];
                        leaving = i;
                    }
                }

                // If no leaving variable, the dual is feasible (primal optimality)
                if (leaving == -1) break;

                // Step 2: Find the entering variable (minimum reduced cost ratio)
                int entering = -1;
                double minRatio = double.MaxValue;
                for (int j = 0; j < n; j++)
                {
                    if (A[leaving, j] > EPS) // Positive coefficient required for pivot
                    {
                        double ratio = c[j] / A[leaving, j];
                        if (ratio < minRatio || (Math.Abs(ratio - minRatio) < EPS && j < entering))
                        {
                            minRatio = ratio;
                            entering = j;
                        }
                    }
                }

                // If no entering variable, the problem is unbounded
                if (entering == -1)
                {
                    Console.WriteLine("Dual problem is unbounded.");
                    return new SimplexResult
                    {
                        feasible = false,
                        solution = null,
                        objVal = double.PositiveInfinity
                    };
                }

                if (IsPivotRepeated(entering, leaving, pivotHistory))
                {
                    Console.WriteLine("Repeated pivot detected in dual simplex. Returning failure.");
                    return new SimplexResult { feasible = false, solution = null, objVal = double.PositiveInfinity };
                }

                // Check for degeneracy
                if (CheckDegeneracy(B, pastStates))
                {
                    Console.WriteLine("Degeneracy detected. Perturbing tableau.");
                    PerturbTableau(ref A, ref B, ref c);
                    continue;
                }

                // Step 3: Perform the pivot operation
                double pivot = A[leaving, entering];
                if (Math.Abs(pivot) < EPS)
                {
                    throw new InvalidOperationException("Pivot element is too small, leading to numerical instability.");
                }

                // Normalize the pivot row
                for (int j = 0; j < n; j++)
                    A[leaving, j] /= pivot;
                B[leaving] /= pivot;

                // Update other rows
                for (int i = 0; i < m; i++)
                {
                    if (i != leaving)
                    {
                        double factor = A[i, entering];
                        for (int j = 0; j < n; j++)
                            A[i, j] -= factor * A[leaving, j];
                        B[i] -= factor * B[leaving];
                    }
                }

                // Update cost coefficients
                double objFactor = c[entering];
                for (int j = 0; j < n; j++)
                    c[j] -= objFactor * A[leaving, j];

                // Update basic variables
                basicVar[leaving] = entering;
            }

            // Build solution
            double[] solution = new double[n];
            for (int i = 0; i < m; i++)
                solution[basicVar[i]] = B[i];

            double objVal = 0.0;
            for (int j = 0; j < n; j++)
                objVal += solution[j] * obj[j];

            return new SimplexResult
            {
                feasible = true,
                solution = solution,
                objVal = objVal
            };

        }

        static bool CheckDegeneracy(double[] rhs, HashSet<string> pastStates)
        {
            string currentState = string.Join(",", rhs.Select(x => x.ToString("F6")));
            if (pastStates.Contains(currentState))
            {
                Console.WriteLine("Degeneracy detected. Forcing perturbation.");
                return true;
            }
            pastStates.Add(currentState);
            return false;
        }

        static void PerturbTableau(ref double[,] A, ref double[] B, ref double[] c)
        {
            Random rand = new Random();
            for (int i = 0; i < B.Length; i++)
            {
                B[i] += rand.NextDouble() * 1e-6; // Small perturbation to RHS
            }

            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    A[i, j] += rand.NextDouble() * 1e-6; // Small perturbation to LHS
                }
            }

            for (int j = 0; j < c.Length; j++)
            {
                c[j] += rand.NextDouble() * 1e-6; // Small perturbation to objective
            }
        }

        static bool IsPivotRepeated(int entering, int leaving, HashSet<(int, int)> pivotHistory)
        {
            if (pivotHistory.Contains((entering, leaving)))
            {
                return true;
            }
            pivotHistory.Add((entering, leaving));
            return false;
        }

        static SimplexResult PrimalSimplexSolveFull(double[,] lhs, double[] rhs, double[] obj)
        {
            int m = rhs.Length; // Number of constraints
            int n = obj.Length; // Number of variables

            // Make local copies
            double[,] A = (double[,])lhs.Clone();
            double[] B = (double[])rhs.Clone();
            double[] c = (double[])obj.Clone();

            int[] basicVar = new int[m];
            for (int i = 0; i < m; i++)
                basicVar[i] = n - m + i; // Initialize slack variables as basic variables

            const double EPS = 1e-9; // Precision threshold
            int iteration = 0;

            while (true)
            {
                iteration++;

                // Step 1: Compute reduced costs
                double[] reducedCosts = new double[n];
                for (int j = 0; j < n; j++)
                {
                    reducedCosts[j] = c[j];
                    for (int i = 0; i < m; i++)
                    {
                        reducedCosts[j] -= c[basicVar[i]] * A[i, j];
                    }
                }

                // Step 2: Check for optimality (all reduced costs >= 0)
                int entering = -1;
                for (int j = 0; j < n; j++)
                {
                    if (reducedCosts[j] < -EPS)
                    {
                        entering = j;
                        break;
                    }
                }
                if (entering == -1) break; // Optimal solution found

                // Step 3: Perform ratio test to find the leaving variable
                int leaving = -1;
                double minRatio = double.MaxValue;
                for (int i = 0; i < m; i++)
                {
                    if (A[i, entering] > EPS)
                    {
                        double ratio = B[i] / A[i, entering];
                        if (ratio < minRatio)
                        {
                            minRatio = ratio;
                            leaving = i;
                        }
                    }
                }

                if (leaving == -1)
                {
                    Console.WriteLine("Primal problem is unbounded.");
                    return new SimplexResult { feasible = false, solution = null, objVal = double.PositiveInfinity };
                }

                // Step 4: Perform the pivot operation
                double pivot = A[leaving, entering];
                for (int j = 0; j < n; j++)
                    A[leaving, j] /= pivot;
                B[leaving] /= pivot;

                for (int i = 0; i < m; i++)
                {
                    if (i != leaving)
                    {
                        double factor = A[i, entering];
                        for (int j = 0; j < n; j++)
                            A[i, j] -= factor * A[leaving, j];
                        B[i] -= factor * B[leaving];
                    }
                }

                double objFactor = c[entering];
                for (int j = 0; j < n; j++)
                    c[j] -= objFactor * A[leaving, j];

                // Update basic variables
                basicVar[leaving] = entering;
            }

            // Build the solution
            double[] solution = new double[n];
            for (int i = 0; i < m; i++)
                solution[basicVar[i]] = B[i];

            double objVal = 0.0;
            for (int j = 0; j < n; j++)
                objVal += solution[j] * obj[j];

            return new SimplexResult
            {
                feasible = true,
                solution = solution,
                objVal = objVal
            };
        }

        static SimplexResult SolveUsingBlandsRule(double[,] lhs, double[] rhs, double[] obj)
        {
            Console.WriteLine("\n[SolveUsingBlandsRule] Start");

            // Dimensions
            int numConstraints = rhs.Length; // Number of constraints
            int totalVars = obj.Length;      // Number of variables

            // Local copies of inputs
            double[,] A = (double[,])lhs.Clone();
            double[] B = (double[])rhs.Clone();
            double[] objCoeffs = (double[])obj.Clone();

            // Initialize basic variables
            int[] basicVars = new int[numConstraints];
            int originalVars = totalVars - numConstraints; // Slack variables start after original vars
            for (int i = 0; i < numConstraints; i++)
                basicVars[i] = originalVars + i;

            // Iterative solving
            bool proceed = true;
            while (proceed)
            {
                // Compute dual variables
                double[] pi = new double[numConstraints];
                for (int i = 0; i < numConstraints; i++)
                    pi[i] = objCoeffs[basicVars[i]];

                // Compute reduced costs
                double[] redCosts = new double[totalVars];
                for (int j = 0; j < totalVars; j++)
                {
                    redCosts[j] = objCoeffs[j];
                    for (int i = 0; i < numConstraints; i++)
                        redCosts[j] -= pi[i] * A[i, j];
                }

                // Find entering variable (Bland's Rule)
                int entering = -1;
                for (int j = 0; j < totalVars; j++)
                {
                    if (redCosts[j] > 1e-9)
                    {
                        entering = j;
                        break;
                    }
                }
                if (entering < 0) break; // Optimal solution found

                // Find leaving variable
                int leaving = -1;
                double minRatio = double.MaxValue;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (A[i, entering] > 1e-9)
                    {
                        double ratio = B[i] / A[i, entering];
                        if (ratio < minRatio)
                        {
                            minRatio = ratio;
                            leaving = i;
                        }
                        else if (Math.Abs(ratio - minRatio) < 1e-9)
                        {
                            if (basicVars[i] < basicVars[leaving]) leaving = i; // Bland's tie-breaking
                        }
                    }
                }
                if (leaving < 0)
                {
                    Console.WriteLine("Problem is unbounded.");
                    return new SimplexResult { feasible = false, solution = null, objVal = double.PositiveInfinity };
                }

                // Perform pivot operation
                double pivot = A[leaving, entering];
                for (int j = 0; j < totalVars; j++)
                    A[leaving, j] /= pivot;
                B[leaving] /= pivot;

                for (int i = 0; i < numConstraints; i++)
                {
                    if (i != leaving)
                    {
                        double fac = A[i, entering];
                        for (int j = 0; j < totalVars; j++)
                            A[i, j] -= fac * A[leaving, j];
                        B[i] -= fac * B[leaving];
                    }
                }

                // Update objective coefficients
                double objFac = objCoeffs[entering];
                for (int j = 0; j < totalVars; j++)
                    objCoeffs[j] -= objFac * A[leaving, j];

                // Update basic variables
                basicVars[leaving] = entering;
            }

            // Build solution
            double[] solution = new double[totalVars];
            for (int i = 0; i < numConstraints; i++)
                solution[basicVars[i]] = B[i];

            double objVal = 0.0;
            for (int j = 0; j < totalVars; j++)
                objVal += solution[j] * obj[j];

            Console.WriteLine("[SolveUsingBlandsRule] Done");
            return new SimplexResult
            {
                feasible = true,
                solution = solution,
                objVal = objVal
            };
        }


        static SimplexResult SolveWithDualAndPrimal(double[,] lhs, double[] rhs, double[] obj)
        {
            Console.WriteLine("Attempting to solve with Dual Simplex...");
            var dualResult = DualSimplexSolveFull(lhs, rhs, obj);
            if (dualResult.feasible)
            {
                Console.WriteLine("Dual Simplex succeeded.");
                return dualResult;
            }

            Console.WriteLine("Dual Simplex failed. Switching to Primal Simplex...");
            var primalResult = PrimalSimplexSolveFull(lhs, rhs, obj);
            if (primalResult.feasible)
            {
                Console.WriteLine("Primal Simplex succeeded.");
                return primalResult;
            }

            Console.WriteLine("Both Dual and Primal Simplex methods failed.");
            return new SimplexResult { feasible = false, solution = null, objVal = double.PositiveInfinity };
        }

        static string DetermineSimplexMethod(double[] rhs, double[] objCoeffs)
        {
            bool all_bi_nonnegative = true;
            bool all_cN_nonpositive = true;
            bool some_bi_zero = false;
            bool some_bi_negative = false;

            foreach (double bi in rhs)
            {
                if (bi < -1e-9) { all_bi_nonnegative = false; some_bi_negative = true; }
                if (Math.Abs(bi) < 1e-9) some_bi_zero = true;
            }
            foreach (double cN in objCoeffs)
            {
                if (cN > 1e-9) all_cN_nonpositive = false;
            }

            if (all_bi_nonnegative && all_cN_nonpositive)
                return "Optimal";
            else if (some_bi_zero)
                return "BlandsRule";
            else if (all_bi_nonnegative && !all_cN_nonpositive)
                return "Primal";
            else if (some_bi_negative && all_cN_nonpositive)
                return "Dual";
            return "Primal";
        }

        static void SolveUsingPrimalSimplex(double[,] lhs, double[] rhs,
                double[] objectiveCoeffs, int totalVars, int numConstraints)
        {
            Console.WriteLine("\n[SolveUsingPrimalSimplex] Start");
            int originalVars = totalVars - numConstraints;
            bool proceed = true;
            double[] objCoeffs = new double[totalVars];
            Array.Copy(objectiveCoeffs, objCoeffs, totalVars);

            int[] basicVars = new int[numConstraints];
            for (int i = 0; i < numConstraints; i++)
            {
                basicVars[i] = originalVars + i;
            }
            while (proceed)
            {
                // compute dual
                double[] pi = new double[numConstraints];
                for (int i = 0; i < numConstraints; i++)
                    pi[i] = objCoeffs[basicVars[i]];
                double[] redCosts = new double[totalVars];
                for (int j = 0; j < totalVars; j++)
                {
                    redCosts[j] = objCoeffs[j];
                    for (int i = 0; i < numConstraints; i++)
                    {
                        redCosts[j] -= pi[i] * lhs[i, j];
                    }
                }
                int entering = -1; double maxRC = 0.0;
                for (int j = 0; j < totalVars; j++)
                {
                    if (redCosts[j] > maxRC)
                    {
                        maxRC = redCosts[j];
                        entering = j;
                    }
                }
                if (entering == -1)
                {
                    proceed = false; // optimal
                    break;
                }
                int leaving = -1;
                double minRatio = double.MaxValue;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (lhs[i, entering] > 1e-12)
                    {
                        double ratio = rhs[i] / lhs[i, entering];
                        if (ratio < minRatio)
                        {
                            minRatio = ratio;
                            leaving = i;
                        }
                    }
                }
                if (leaving == -1)
                {
                    proceed = false; // unbounded
                    break;
                }
                double pivot = lhs[leaving, entering];
                for (int j = 0; j < totalVars; j++)
                    lhs[leaving, j] /= pivot;
                rhs[leaving] /= pivot;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (i != leaving)
                    {
                        double fac = lhs[i, entering];
                        for (int j = 0; j < totalVars; j++)
                            lhs[i, j] -= fac * lhs[leaving, j];
                        rhs[i] -= fac * rhs[leaving];
                    }
                }
                double objFac = objCoeffs[entering];
                for (int j = 0; j < totalVars; j++)
                {
                    objCoeffs[j] -= objFac * lhs[leaving, j];
                }
                basicVars[leaving] = entering;
            }
            for (int j = 0; j < totalVars; j++)
                objectiveCoeffs[j] = objCoeffs[j];
            Console.WriteLine("[SolveUsingPrimalSimplex] Done");
        }

        static void SolveUsingDualSimplex(double[,] lhs, double[] rhs,
                double[] objCoeffs, int totalVars, int numConstraints)
        {
            Console.WriteLine("\n[SolveUsingDualSimplex] Start");
            int originalVars = totalVars - numConstraints;
            int[] basicVars = new int[numConstraints];
            for (int i = 0; i < numConstraints; i++)
                basicVars[i] = originalVars + i;
            bool proceed = true;
            while (proceed)
            {
                int leavingRow = -1; double mostNeg = 0;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (rhs[i] < mostNeg)
                    {
                        mostNeg = rhs[i];
                        leavingRow = i;
                    }
                }
                if (leavingRow < 0) break; // optimum
                int enteringCol = -1;
                double minRatio = double.MaxValue;
                for (int j = 0; j < totalVars; j++)
                {
                    if (lhs[leavingRow, j] < -1e-12)
                    {
                        double reducedCost = objCoeffs[j];
                        for (int ii = 0; ii < numConstraints; ii++)
                            reducedCost -= objCoeffs[basicVars[ii]] * lhs[ii, j];
                        double ratio = reducedCost / lhs[leavingRow, j];
                        if (ratio < minRatio)
                        {
                            minRatio = ratio;
                            enteringCol = j;
                        }
                    }
                }
                if (enteringCol < 0)
                {
                    proceed = false; // infeasible
                    break;
                }
                double pivot = lhs[leavingRow, enteringCol];
                for (int j = 0; j < totalVars; j++)
                    lhs[leavingRow, j] /= pivot;
                rhs[leavingRow] /= pivot;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (i != leavingRow)
                    {
                        double fac = lhs[i, enteringCol];
                        for (int j = 0; j < totalVars; j++)
                            lhs[i, j] -= fac * lhs[leavingRow, j];
                        rhs[i] -= fac * rhs[leavingRow];
                    }
                }
                double objFac = objCoeffs[enteringCol];
                for (int j = 0; j < totalVars; j++)
                    objCoeffs[j] -= objFac * lhs[leavingRow, j];
                basicVars[leavingRow] = enteringCol;
            }
            Console.WriteLine("[SolveUsingDualSimplex] Done");
        }

        static void SolveUsingBlandsRule(double[,] lhs, double[] rhs,
                double[] originalObj, int totalVars, int numConstraints)
        {
            Console.WriteLine("\n[SolveUsingBlandsRule] Start");
            double[] objCoeffs = new double[totalVars];
            Array.Copy(originalObj, objCoeffs, totalVars);
            int originalVars = totalVars - numConstraints;
            int[] basicVars = new int[numConstraints];
            for (int i = 0; i < numConstraints; i++)
                basicVars[i] = originalVars + i;
            bool proceed = true;
            while (proceed)
            {
                double[] pi = new double[numConstraints];
                for (int i = 0; i < numConstraints; i++)
                    pi[i] = objCoeffs[basicVars[i]];
                double[] redCosts = new double[totalVars];
                for (int j = 0; j < totalVars; j++)
                {
                    redCosts[j] = objCoeffs[j];
                    for (int i = 0; i < numConstraints; i++)
                        redCosts[j] -= pi[i] * lhs[i, j];
                }
                int entering = -1;
                for (int j = 0; j < totalVars; j++)
                {
                    if (redCosts[j] > 1e-9)
                    {
                        entering = j;
                        break;
                    }
                }
                if (entering < 0) break; // optimum
                int leaving = -1; double minRatio = double.MaxValue;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (lhs[i, entering] > 1e-9)
                    {
                        double ratio = rhs[i] / lhs[i, entering];
                        if (ratio < minRatio)
                        {
                            minRatio = ratio;
                            leaving = i;
                        }
                        else if (Math.Abs(ratio - minRatio) < 1e-9)
                        {
                            if (basicVars[i] < basicVars[leaving]) leaving = i;
                        }
                    }
                }
                if (leaving < 0)
                {
                    proceed = false; // unbounded
                    break;
                }
                double pivot = lhs[leaving, entering];
                for (int j = 0; j < totalVars; j++)
                    lhs[leaving, j] /= pivot;
                rhs[leaving] /= pivot;
                for (int i = 0; i < numConstraints; i++)
                {
                    if (i != leaving)
                    {
                        double fac = lhs[i, entering];
                        for (int j = 0; j < totalVars; j++)
                            lhs[i, j] -= fac * lhs[leaving, j];
                        rhs[i] -= fac * rhs[leaving];
                    }
                }
                double objFac = objCoeffs[entering];
                for (int j = 0; j < totalVars; j++)
                    objCoeffs[j] -= objFac * lhs[leaving, j];
                basicVars[leaving] = entering;
            }
            for (int j = 0; j < totalVars; j++)
                originalObj[j] = objCoeffs[j];
            Console.WriteLine("[SolveUsingBlandsRule] Done");
        }



        // Helper function to print a matrix
        static void PrintMatrix(double[,] matrix, int rows, int cols)
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    Console.Write($"{matrix[i, j]:F2} ");
                }
                Console.WriteLine();
            }
        }

        // Helper function to print an array
        static void PrintArray(double[] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                Console.Write($"{array[i]:F2} ");
            }
            Console.WriteLine();
        }
        static void PrintOriginalProblem()
        {
            Console.WriteLine("--- Original Farmer's Problem (No ILOG)---");
            Console.WriteLine("Objective: min 150x_1 +230x_2 +260x_3 +238y_1 -170w_1 +210y_2 -150w_2 -36w_3 -10w_4");
            Console.WriteLine("Constraints:");
            Console.WriteLine("C1: x_1 + x_2 + x_3 ≤ 500");
            Console.WriteLine("C2: 2.5x_1 + y_1 - w_1 ≥ 200 => -2.5x_1 - y_1 + w_1 ≤ -200");
            Console.WriteLine("C3: 3x_2 + y_2 - w_2 ≥ 240 => -3x_2 - y_2 + w_2 ≤ -240");
            Console.WriteLine("C4: w_3 + w_4 ≤ 20x_3 => w_3 + w_4 -20x_3 ≤ 0");
            Console.WriteLine("C5: w_3 ≤ 6000");
        }

        static void PrintConstraint(double[,] LHS, int row, List<string> varNames, double rhs)
        {
            bool first = true;
            for (int j = 0; j < varNames.Count; j++)
            {
                double coeff = LHS[row, j];
                if (Math.Abs(coeff) > 1e-15)
                {
                    if (!first && coeff > 0) Console.Write("+ ");
                    Console.Write($"{coeff:F2}*{varNames[j]} ");
                    first = false;
                }
            }
            Console.WriteLine($" ≤ {rhs:F2}");
        }

        static double[] ComputeScenarioMultipliers(int numScenarios)
        {
            double[] arr = new double[numScenarios];
            Random rand = new Random();

            for (int i = 0; i < numScenarios; i++)
            {
                // Generate a random double uniformly in the range [0.8, 1.2]
                arr[i] = 0.8 + (1.2 - 0.8) * i / (numScenarios - 1);
            }

            return arr;
        }

    }

    public class FeasibilityHelper
    {
        // Get the Left-Hand Side matrix (LHS) for the Simplex method
        public static double[,] GetLHS(double[,] W)
        {
            int rows = W.GetLength(0); // Number of constraints (m)
            int cols = W.GetLength(1); // Number of y variables (n)

            // New matrix dimensions: m x (n + m)
            // Typically, for dual variables, you might not need to append additional variables
            double[,] extendedMatrix = new double[rows, cols];

            // Step 1: Copy W into the extended matrix
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    extendedMatrix[i, j] = W[i, j];
                }
            }

            return extendedMatrix;
        }

        // Get the Right-Hand Side (RHS) vector for the Simplex method
        public static double[] GetRHS(double[] h_k, double[,] T_k, double[] x)
        {
            int rows = h_k.Length;  // Number of constraints (rows of h_k)
            int cols = x.Length;    // Number of variables (columns of T_k)

            // Validate dimensions
            if (T_k.GetLength(0) != rows || T_k.GetLength(1) != cols)
            {
                throw new ArgumentException("Dimensions of T_k and x do not match with h_k.");
            }

            // Compute T_k * x using the helper method
            double[] TkX = MultiplyMatrixVector(T_k, x);

            // Compute h_k - T_k * x
            double[] rhs = new double[rows];
            for (int i = 0; i < rows; i++)
            {
                rhs[i] = h_k[i] - TkX[i];
            }

            return rhs;
        }

        // Helper method for matrix-vector multiplication
        public static double[] MultiplyMatrixVector(double[,] matrix, double[] vector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            double[] result = new double[rows];

            for (int i = 0; i < rows; i++)
            {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++)
                {
                    result[i] += matrix[i, j] * vector[j];
                }
            }

            return result;
        }

        // Get the Objective Coefficients vector for the Simplex method
        public static double[] GetObjectiveCoefficients(int numConstraints)
        {
            // The length of the coefficient vector is numConstraints
            double[] objCoeffs = new double[numConstraints];

            // Initialize the coefficient vector with all values set to -1 (for minimization)
            for (int i = 0; i < numConstraints; i++)
            {
                objCoeffs[i] = -1.0;
            }

            return objCoeffs;
        }
    }

}

