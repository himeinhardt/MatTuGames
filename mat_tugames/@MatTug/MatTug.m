classdef MatTug < TuGame
% MatTuGames: A Matlab Game Theory Toolbox 
% Version 1.0 (R2015b) 24-Dec-2018
%
% aux: Some auxiliary files
%--------------------------
% FrameToImage                        - Converts a frame to an image.
% PlayCoreMovie                       - Plays a movie from a collection of frames.
% SaveFrames                          - Saves converted image files to the hard disk.
% ToSimplex3d                         - Projects data from 4d to 3d.
%
% 
% bin: Scirpt File
%----------------
% corevert                            - External bash script to call the cdd library. 
%
%
%
% doc: Document Files
%-------------------
% getting_started.m                   - Checks the installation
% getting_started.out                 - Reference results of getting_started
% MatTuGames_Version_0.8.m            - Additions and changes in version 0.8
% manual_mat_tugames.pdf              - Manual (PDF)
% graphics_manual.pdf                 - Manual for using the graphical features.
% ReadMe.pdf                          - Installation instruction (PDF)
% ReadMe.ps                           - Installation instruction (PS)
% ReadMe.txt                          - Installation instruction (TXT)
% testcase_graphics                   - Checking basic graphic installation.

%
%
%
%
% graphics: Graphic Example Files
%-------------------------------
% core_exp_all.pdf                    - Core plot example 3d.
% core_exp.pdf                        - Core plot example 3d.
% core_exp_prk.pdf                    - Core plot example 3d.
% core_exp_prn.pdf                    - Core plot example 3d.
% core_exp_shap.pdf                   - Core plot example 3d.
% core_exp_sol_all.pdf                - Core plot example 3d.
% core_exp_sol_none.pdf               - Core plot example 3d.
% core_exp_sol_prk.pdf                - Core plot example 3d.
% core_exp_sol_prn.pdf                - Core plot example 3d.
% core_exp_sol_shap.pdf               - Core plot example 3d.
% manual_exp2_core01.pdf              - Core plot example 3d.
%
%
%
% mama: Mathematica Symbolic Toolbox Functions to call the Mathematica Package TuGames
%------------------------------------------------------------------------------------
% tug_AdjustedWorthVectors            - Computes the adjusted worth vectors of k-convex games.   
% tug_AllAntiSurpluses                - Computes the minimum surpluses.
% tug_AllMaxSurpluses                 - Computes the maximum surpluses.
% tug_AntiPreKernel                   - Computes an anti-pre-kernel point.
% tug_AntiPreKernelQ                  - Checks if an imputation is an anti-pre-kernel point.
% tug_AvConvexQ                       - Checks on average convexity.
% tug_AverageConvexQ                  - Checks on average convexity.
% tug_BalancedKSelectionQ             - Checks if an imputation induces a k-balanced selection.
% tug_BalancedSelectionQ              - Checks if an imputation induces a balanced selection.
% tug_Bankruptcy                      - Creates a modest bankruptcy game.
% tug_BelongToCoreQ                   - Checks if an imputation belongs to the core.
% tug_BestCoalToMatrix                - Computes an equivalence matrix. 
% tug_Bsc                             - Returns the set of most effective coalitions.
% tug_CharacteristicValues            - Computes the characteristic values.
% tug_Coal2Dec                        - List of proper coalitions in Mathematica order.
% tug_CollectionBalancedQ             - Checks if a collection is balanced.
% tug_CollectionOfDecreasingExcess    - Creates the collection of decreasing excesses.
% tug_Concession                      - Computes the concession vector.
% tug_ContestedGarment                - Computes the contested garment. 
% tug_ConvexQ                         - Checks convexity.
% tug_ConvexUnanConditionQ            - Checks convexity while relying on the unanimity coordinates.
% tug_CoreElementsQ                   - Checks if an imputation belongs to the core.
% tug_CoreQ                           - Checks if the core is non-empty. 
% tug_CostSavings                     - Creates the cost savings game.
% tug_CriticalVal                     - Computes some critical epsilon values.
% tug_DetQuasiAvConvex                - Determines a quasi average convex game.
% tug_DetRandCoord                    - Returns random unanimity coordinates.
% tug_DetUCoord                       - Determines the missing unanimity coordinates of size greater than 2. 
% tug_Disagreement                    - Computes the disagreement vector.
% tug_DualGame                        - Creates the dual of a Tu-game.
% tug_EpsCore                         - Computes the least core.
% tug_EqClass                         - Determines the equivalence classes from the set of most effective coalitions.
% tug_EvalSumMinCoord                 - Calculates at most (n-1) inequalities of the unanimity coordinates constraints of nonnegative sums.
% tug_ExcessValues                    - Determines the excesses.
% tug_FindPreKernel                   - Computes a pre-kernel element.
% tug_GameMonotoneQ                   - Checks on monotonicity.
% tug_Gap                             - Computes the gap function.
% tug_GrandCoalitionLargestValueQ     - Checks if the grand coalition has largest value.
% tug_GreedyBankruptcy                - Creates the greedy bankruptcy game.
% tug_HarsanyiDividends               - Creates the unanimity coordinates.
% tug_ImpToVec                        - Converts an imputation to a set of vectors.
% tug_ImputationQ                     - Checks if a payoff vector is an imputation.
% tug_IntersectionOfMaxExcessSets     - Determines if the set of proper coalitions having largest excesses has an empty intersection.
% tug_IntersectionUpperLowerSetQ      - Checks if the intersection of the lower and upper set is non-empty. 
% tug_kCover                          - Determines from the Tu-game the corresponding k-game.
% tug_KernelCalculation               - Computes a or some kernel element(s).
% tug_KernelImputationQ               - Checks if an imputation is a kernel point.
% tug_Kernel                          - Computes a kernel point.
% tug_KernelVertices                  - Computes a kernel segment.
% tug_LargestAmount                   - Computes the largest amount.
% tug_LeastCore                       - Determine the least core.
% tug_LexiCenter                      - Computes the lexi center.
% tug_LowerSetIncImputationQ          - Checks if the lower set is included in the imputation set.
% tug_LowerSetQ                       - Checks if an imputation belongs to the lower set.
% tug_MargValue                       - Determines the marginal contribution vector.
% tug_MaxExcessBalanced               - Checks if the maximum surpluses are balanced.
% tug_MaxExcessSets                   - Computes the set of proper coalitions having largest excesses.
% tug_MinExcessBalanced               - Determines if the minimum surpluses are balanced.
% tug_MinUnanimityCoordinates         - Returns the minimum unanimity coordinates. 
% tug_MKernel                         - Determines a kernel point.
% tug_MLExtension                     - Computes the multi-linear extension.
% tug_Mnuc                            - Determines the nucleolus.
% tug_MonotoneQ                       - Checks on monotonicity.
% tug_Nuc                             - Computes the nucleolus.
% tug_OneNormalization                - Creates a one normalized game.
% tug_PreKernelEl                     - Computes a pre-kernel element.
% tug_PreKernelEqualsKernelQ          - Checks if the pre-kernel coincides with the kernel.
% tug_PreKernel                       - Computes a pre-kernel element.
% tug_PreKernelQ                      - Checks if an imputation is a pre-kernel element.
% tug_PreNuc                          - Computes the pre-nucleolus.
% tug_ProperAmount                    - Computes the proper amount.
% tug_Quota                           - Computes the quotas.
% tug_ReasonableOutcome               - Computes the reasonable outcome.
% tug_ReasonableSet                   - Computes the reasonable set.
% tug_ScrbSolution                    - Determines the Scrb solution.
% tug_SetsToVec                       - Converts the set of most effective coalitions to a set of vectors.
% tug_ShapleyValue                    - Determines the Shapley value.
% tug_ShapleyValueML                  - Determines the Shapley value using multi-linear extension.
% tug_SmallestContribution            - Determines the smallest contribution vector.
% tug_StrictlyConvexUnanConditionQ    - Examines the sufficient condition of convexity in terms of unanimity coordinates.
% tug_SuperAdditiveQ                  - Checks on super-additivity.
% tug_SymGameSizeK                    - Returns a special type of symmetric game.
% tug_SymGameType2                    - Returns a special type of symmetric game.
% tug_SymGameType3                    - Returns a special type of symmetric game.
% tug_SymGameType4                    - Returns a special type of symmetric game.
% tug_TalmudicRule                    - Computes the Talmudic distribution rule.
% tug_TauValue                        - Determines the Tau value.
% tug_UnanAvConvexQ                   - Checks if the coordinates satisfy the sufficient and necessary condition of average convexity.
% tug_UnanConvexQ                     - Checks if the coordinates satisfy the sufficient and necessary condition of convexity.
% tug_UnanimityCoordinates            - Determines all unanimity coordinates of the game
% tug_UpperSetIncImputationQ          - Checks if the upper set is included in the imputation set.
% tug_UpperSetQ                       - Checks if an imputation belongs to the upper set.
% tug_UtopiaVector                    - Computes the utopia payoff. 
% tug_ValueExcess                     - Computes an objective function to compute a pre-kernel element.
% tug_VerticesCore                    - Determines the vertices of the core.
% tug_WeaklySuperAdditiveQ            - Checks if the Tu-game is weakly super-additive.
% tug_WeightedMajority                - Creates the weighted majority game.
% tug_ZeroMonotoneQ                   - Checks on zero-monotonicity.
% tug_ZeroNormalization               - Creates the zero normalized game.
% tug_ZeroOneNormalization            - Creates the zero-one normalized game.
%
%
%
% mat_tugames: Serial Computing
%-----------------------------
% additive_game                       - Creates an additive game.
% AdditiveQ                           - Checks if the game v is additive.
% admissibleGame                      - Computes a symmetric compromise admissible game.
% ADvalue                             - Computes the Aumann-Dreze value.
% AllMarginalContributions            - Computes all marginal contributions of a Tu game.
% AllSubGames                         - Computes all subgames.
% Anti_balancedCollectionQ            - Checks the reversal of Kohlberg's criterion
% AntiCorePlot                        - Plots the anti-core.
% AntiCoreVertices                    - Evaluates the vertices of the anti-core.
% AntiImputationVertices              - Computes all vertices of the anti imputation set.
% anti_coreQ                          - Checks the existence for the anti core of the game.
% Anti_Kernel                         - Computes an anti-kernel point.
% Anti_kernelQ                        - Checks if an imputation is an anti kernel point.
% Anti_PreKernel                      - Computes an anti-prekernel point.
% Anti_PrekernelQ                     - Checks if an imputation is an anti prekernel point.
% Anti_PreNucl                        - Computes the anti pre-nucleolus of game v.
% apex_game                           - Creates an apex game.
% assignment_game                     - Creates an assignment game.
% average_convexQ                     - Checks the Tu-game on average convexity.
% balancedCollectionQ                 - Checking Kohlberg's criterion.
% balancedQ                           - Verifies whether the collection of coalitions is balanced.
% bankruptcy_game                     - Creates a bankruptcy game.
% banzhaf                             - Computes the Banzhaf value.
% basis_coordinates                   - Determines the basis coordinates of a Tu game.
% basis_game                          - Determines bases games.
% belongToAntiCoreQ                   - Checks if a payoff vector belongs to the anti-core.
% belongToCoreQ                       - Checks if a payoff vector belongs to the core.
% belongToImputationSetQ              - Checks if a payoff vector belongs to imputation set.
% belongToLeastCoreQ                  - Checks if a payoff vector belongs to the least core.
% belongToLowerSetQ                   - Checks if a payoff vector belongs to lower set.
% belongToUpperSetQ                   - Checks if a payoff vector belongs to upper set.
% BestCoalitions                      - Computes the set of most effective coalitions.
% cardinality_game                    - Assigns zero to a coalition of size<=k<n, otherwise its cardinality.
% CddAntiCorePlot                     - Plots the anti-core of a game using cddmex. 
% CddAntiCoreQ                        - Checks if the anti-core exists (cddmex).
% CddAntiCoreSimplexPlot              - Plots the anti-core using simplex projection.
% CddAntiCoreSimplexVertices          - Computes all anti-core vertices using simplex projection.
% CddAntiCoreVertices                 - Computes the vertices of the anti-core (cddmex).
% CddAntiImputationSimplexVertices    - Computes all vertices of the anti-imputation set using simplex projection.
% CddAntiImputationVertices           - Computes all vertices of the anti-imputation set.
% CddAntiPrenucl                      - Computes the anti pre-nucleolus of game v (cddmex).
% CddBelongToLeastCoreQ               - Checks if a payoff vector belongs to the least-core.
% CddCoreCoverPlot                    - Plots the core cover of a TU game.
% CddCoreCoverSimplexPlot             - Plots the core cover (simplex projection).
% CddCoreCoverSimplexVertices         - Computes all vertices of the core cover (simplex).
% CddCoreCoverVertices                - Computes all vertices of the core cover of a TU game.
% CddCoreMovie                        - Creates a movie w.r.t. the strong epsilon-cores.
% CddCorePlot                         - Plots the core of a game using cddmex. 
% CddCoreQ                            - Checks if the core exists (cddmex).
% CddCoreSimplexMovie                 - Creates a movie w.r.t. the strong epsilon-cores (simplex projection).
% CddCoreSimplexPlot                  - Plots the core (simplex projection).
% CddCoreSimplexVertices              - Computes the vertices of the core (simplex).
% CddCoreVertices                     - Computes the vertices of the core (cddmex).
% CddImputationSimplexVertices        - Computes the vertices of the imputation set (simplex).
% CddImputationVertices               - Computes the vertices of the imputation set (cddmex).
% CddKernelCatchers                   - Draws some kernel catchers (cddmex).
% CddKernelCatchersSimplex            - Draws some kernel catchers (simplex).
% CddLeastCore                        - Computes the leastcore (cddmex).
% CddLeastCoreVertices                - Computes the leastcore vertices (cddmex). 
% CddLowerSetSimplexVertices          - Computes the vertices of the lower set (simplex). 
% CddLowerSetVertices                 - Computes the vertices of the lower set (cddmex). 
% CddNucl                             - Computes the nucleolus using the CDD solver (cddmex).
% CddPreKernel                        - Computes a pre-kernel element (cddmex).
% CddPrenucl                          - Computes the prenucleolus using the CDD solver (cddmex).
% CddReasonableSetSimplexVertices     - Computes the vertices of the reasonable set (simplex). 
% CddReasonableSetVertices            - Computes the vertices of the reasonable set (cddmex). 
% CddStrongCorePlot                   - Plots a strong epsilon core. 
% CddStrongCoreSimplexPlot            - Plots the strong epsilon core (simplex projection).
% CddUpperSetSimplexVertices          - Computes the vertices of the upper set (simplex). 
% CddUpperSetVertices                 - Computes the vertices of the upper set (cddmex). 
% CddWeberSet                         - Computes the vertices of the Weber Set.
% CddWeberSetPlot                     - Plots the Weber set.
% CddWeberSetSimplex                  - Computes the vertices of the Weber Set (simplex).
% CddWeberSetSimplexPlot              - Plots the Weber set (simplex).
% clp_kernel                          - Computes a kernel point using the CLP solver.
% cls_kernel                          - Computes a kernel point using the CLS solver.
% clToMatlab                          - Computes the unique integer representation of coalitions.
% CoalitionSolidarity                 - Determines the coalition solidarity value.
% coeff_linearbasis                   - Determines the coefficients (dividends) of a linear basis from a TU game.
% compromiseAdmissibleQ               - Checks if the core cover a TU game v is non-empty.
% compromiseStableQ                   - Checks if the game is compromise stable.
% Converse_RGP_Q                      - Checks if an imputation satisfies the CRGP.
% convex_gameQ                        - Checks the convexity of a Tu-game.
% CoreCoverQ                          - Checks if the core cover a TU game v is non-empty.
% CorePlot                            - Plots the core.
% coreQ                               - Checks the non-emptiness of the core.
% CoreVertices                        - Computes the vertices of the core.
% cplex_AntiPreNucl                   - Computes the anti prenucleolus using the CPLEX solver.
% cplex_kernel                        - Computes a kernel point using the CPLEX solver.
% cplex_LeastCore                     - Computes the least core using cplexmex.
% cplex_nucl                          - Computes the nucleolus using the CPLEX solver.
% cplex_prekernel                     - Computes a prekernel point using the CPLEX solver.
% cplex_prenucl                       - Computes the prenucleolus using the CPLEX solver.
% critical_value1                     - Computes the biggest gain of any group of players.
% critical_value2                     - Computes a critical value w.r.t. the strong epsilon-core.
% critical_value_star                 - Computes a critical value which contains the intersection of the imputation and reasonable set
% COV_propertyQ                       - Verifies if the payoff x satisfies COV property.
% cvx_kernel                          - Computes a kernel point using the CVX solver. 
% cvx_prekernel                       - Computes a prekernel point using the CVX solver. 
% DiscShapleyValue                    - Computes the discounted Shapley value.
% DM_Reduced_game                     - Computes all Davis-Maschler reduced games.
% dual_game                           - Creates the dual of a Tu-game.
% excess                              - Determines the excesses w.r.t. a payoff vector. 
% game_basis                          - Computes a game basis of the n-person TU game space.
% game_space                          - Computes the game space which replicates a payoff as a pre-kernel element.
% gameToMama                          - Converts a TU-game into Mathematica representation.
% gameToMatlab                        - Converts a Tu-game into Matlab representation.
% Gap                                 - Determines the gap function.
% genUnionStable                      - Creates a union stable system.
% getgame                             - Creates a Tu-game from the unanimity coordinates.
% glpk_AntiPreNucl                    - Computes the anti pre-nucleolus using the GLPK solver.
% glpk_kernel                         - Computes a kernel point using the GLPK solver. 
% glpk_nucl                           - Computes the nucleolus using the GLPK solver.
% glpk_prekernel                      - Computes a prekernel point using the GLPK solver. 
% glpk_prenucl                        - Computes the prenucleolus using the GLPK solver.
% greedy_bankruptcy                   - Creates the greedy bankruptcy game.
% gurobi_AntiPreNucl                  - Computes the anti prenucleolus using the GUROBI solver.
% gurobi_kernel                       - Computes a kernel point using the GUROBI solver. 
% gurobi_nucl                         - Computes the nucleolus using the GUROBI solver.
% gurobi_prekernel                    - Computes a prekernel point using the GUROBI solver. 
% gurobi_prenucl                      - Computes the prenucleolus using the GUROBI solver.
% harsanyi_dividends                  - Determines the the unanimity coordinates.
% HMS_ImputSavingReducedGame          - Computes from (v,x) all Hart/Mas-Colell ISR games.
% HMS_RedGame                         - Creates a Hart/Mas-Colell reduced games. 
% HMS_Reduced_game                    - Creates all Hart/Mas-Colell reduced games. 
% holler                              - Computes the Holler index.
% homogeneous_representationQ         - Checks if the weighted majority game possesses a homogeneous representation.
% hsl_prekernel                       - Computes a prekernel point using HSL solvers. 
% hypergraphQ                         - Checks whether the system is a hypergraph communication situation.
% ImputationVertices                  - Computes the vertices of the imputation set.
% ImputSavingReducedGame              - Computes from (v,x) all imputation saving reduced games.
% ImpSetEqsLwsQ                       - Checks if the imputation set coincides with the lower set.
% InteractionSets                     - Determines a system of interaction sets.
% interest_game                       - Computes from an interest problem the corresponding game.
% ipopt_kernel                        - Computes a kernel point using the IPOPT solver. 
% ipopt_prekernel                     - Computes a prekernel point using the IPOPT solver.
% ISRG_propertyQ                      - Checks whether an imputation x satisfies the ISR game property.
% KrEqsPrkQ                           - Checks if the kernel is equal to the pre-kernel.
% k_Converse_RGP_Q                    - Checks if an imputation satisfies the k-CRGP.
% k_convexQ                           - Checks k-convexity of the Tu-game.
% k_cover                             - Determines from the Tu-game the corresponding k-game.
% Kernel                              - Computes a kernel point using optimization toolbox.
% kernelQ                             - Checks if an imputation is a kernel point.
% k_Reconfirmation_propertyQ          - Checks the k-RCP.
% k_Reduced_game_propertyQ            - Checks the k-RGP.
% k_StrConverse_RGP_Q                 - Checks the strong k-CRGP. 
% LeastCore                           - Computes the leastcore using optimization toolbox.
% linear_basis                        - Determines the linear basis of the n-person TU game space.
% lin_prekernel                       - Computes a prekernel point using optimization toolbox.
% lowersetQ                           - Checks the existence of the lower set. 
% LS_Nucl                             - Computes the least square nucleolus of a game.
% LS_PreNucl                          - Computes the least square pre-nucleolus of a game.
% MLextension                         - Computes the multi-linear extension.
% MyersonValue                        - Computes the Myerson value of a Tu game.
% market_game                         - Determines from two disjoint sets a market game.
% minimal_winning                     - Computes the minimal winning coalitions.
% monotone_gameQ                      - Checks monotonicity of the TU game.
% monotonic_cover                     - Determines the monotinic cover from a TU game.
% msk_alphaVector                     - Computes recursively an alpha vector from the positive core.
% msk_AntiPreNucl                     - Computes the anti prenucleolus using the MOSEK solver.
% msk_kernel                          - Computes a kernel point using the MOSEK solver.
% msk_nucl                            - Computes the nucleolus using the MOSEK solver.
% msk_prekernel                       - Computes a prekernel point using the MOSEK solver.
% msk_prenucl                         - Computes the prenucleolus using the MOSEK solver.
% nucl                                - Computes the nucleolus using optimization toolbox.
% OwenValue                           - Computes the Owen value.
% oases_kernel                        - Computes a kernel point using the OASES solver.
% oases_prekernel                     - Computes a prekernel point using the OASES solver.
% ols_prekernel                       - Computes a prekernel point using optimization toolbox.
% PartitionSA                         - Computes a partition of S w.r.t. a hypergraph communication situation.
% PartitionSL                         - Computes a partition of S w.r.t. a communication situation.
% PermutationGame                     - Computes from an assignment matrix the permutation game.
% PositionValue                       - Computes the position value.
% Potential                           - Determines the potential of a TU game (recursive).
% PowerSet                            - Computes all subsets from a set representation.
% PreKernel                           - Computes a prekernel element.
% PrekernelQ                          - Checks if an imputation is a pre-kernel point.
% PreNucl                             - Computes the prenucleolus using optimization toolbox.
% PreNucl2                            - Computes the prenucleolus using optimization toolbox.
% PrenuclQ                            - Checks if an imputation is the pre-nucleolus using Kohlberg's criterion.
% production_game                     - Creates an affine production game.
% production_game_sq                  - Creates a quadratic production game. 
% profit_matrix                       - Creates the profit matrix of an assignment game. 
% proper_amount                       - Computes the largest amount players contribute to a proper coalition.
% PropNucl                            - Computes the proportional nucleolus.
% PropPreNucl                         - Computes the proportional pre-nucleolus.
% pure_overhead                       - Creates the matrix of pure overhead games.
% qpBB_kernel                         - Computes a kernel point using the QPBB solver.
% qpc_kernel                          - Computes a kernel point using the QPC solver.
% qpc_prekernel                       - Computes a prekernel point using the QPC solver.
% qrg_prekernel                       - Computes a prekernel point using qrginv instead of pinv.
% quotas                              - Determines the quotas of a game.
% reasonable_outcome                  - Determines the reasonable outcome.
% ReasSetEqsUpsQ                      - Checks if the reasonable set coincides with the upper set.
% Reconfirmation_propertyQ            - Checks the RCP.
% RedGame                             - Creates a Davis-Maschler reduced game.
% Reduced_game_propertyQ              - Checks the RGP.
% replicate_prk                       - Replicates a pre-kernel solution as a pre-kernel of a game space. 
% replicate_Shapley                   - Replicates the Shapley value for a game space. 
% savings_game                        - Creates a saving game from a cost game.
% scrb_solution                       - Computes separable costs-remaining benefits allocation. 
% select_starting_pt                  - Selects a starting point for the pre-kernel computation.
% semi_convexQ                        - Checks semi-convexity.
% separable_cost_allocation           - Computes the separable cost allocation.
% ShapleyQ                            - Checks if the imputation x is a Shapley value of game v.
% ShapleyValue                        - Computes the Shapley value (potential).
% ShapleyValueLB                      - Computes the Shapley value from the linear basis.
% ShapleyValueM                       - Computes the Shapley value based on all marginal contributions.
% ShapleyValueML                      - Computes the Shapley value using multi-linear extension. 
% simple_game                         - Creates a simple game.
% smallest_amount                     - Computes the smallest amount vector.
% SolidarityShapleyValue              - Determines the solidarity Shapley value. 
% SolidarityValue                     - Determines the solidarity value. 
% SortMg                              - Sorts a sub/power set w.r.t. its cardinality.
% sortsets                            - Sorts a sub/power set w.r.t. its cardinality.
% StandardSolution                    - Determines the standard solution.
% StrConverse_RGP_Q                   - Checks the strong RGP.
% streps_value                        - Determines the strong epsilon-game.
% SubCoalitions                       - Computes the power set (subsets) from an array.
% SubDual                             - Determines the dual of a subgame.
% SubGame                             - Creates a subgame.
% SubSets                             - Creates all subsets of super set.
% Sum_Marg_Contributions              - Returns 1 whenever for a coalition the sum of marginal contributions is positive.
% super_additiveQ                     - Checks the Tu-game on super additivity.
% Talmudic_Rule                       - Computes the Talmudic rule.
% TauValue                            - Computes the Tau value.
% unanimity_games                     - Computes the unanimity coordinates.
% UnionStableBasis                    - Determines a basis of a union stable system.
% union_stableQ                       - Checks whether a system is union stable.
% uppersetQ                           - Checks the existence of the upper set.
% UtopiaPayoff                        - Computes the utopia and minimum claim vector of game v.
% vclToMatlab                         - Computes a Tu-game and the corresponding unique integer representation of coalitions
% veto_players                        - Determines the veto players of a simple game. 
% weakly_super_additiveQ              - Checks the Tu-game on weakly super additivity.
% zero_monotonicQ                     - Checks zero monotonicity.
% zero_normalization                  - Creates a zero normalized game.
% ZeroOne_Normalization               - Creates a zero-one normalized game.
%
%
% Class Objects
% -------------
% TuGame                              - to perform several computations for retrieving and modifying game data.
% TuProp                              - subclass object of TuGame (game properties).
% TuSol                               - subclass object of TuGame (game solutions).
% TuVal                               - subclass object of TuGame (fairness and related values).
% TuCore                              - subclass object of TuSol (core plot).
% TuACore                             - subclass object of TuSol (anti-core plot).
% TuAVert                             - subclass object of TuSol (anti-core vertices).
% TuVert                              - subclass object of TuSol (core vertices).
% TuRep                               - subclass object of TuSol (prk replication).
% TuShRep                             - subclass object of TuSol (Shapley value replication).
% TuCons                              - subclass object of TuSol (consistency).
% TuKcons                             - subclass object of TuSol (generalized consistency).
%
%
% pct_tugames: Parallel Computing
%-------------------------------
% p_ADvalue                           - Computes the Aumann-Dreze value.
% p_AllMarginalContributions          - Computes all marginal contributions of a Tu-game.
% p_Anti_PreKernel                    - Computes an anti-pre-kernel element.
% p_Anti_PrekernelQ                   - Checks if an imputation is an anti prekernel point.
% p_assignment_game                   - Creates  an assignment game.
% p_average_convexQ                   - Checks on average convexity.
% p_banzhaf                           - Computes  the Banzhaf value.
% p_basis_coordinates                 - Determines the basis coordinates of a Tu game.
% p_basis_game                        - Determines bases games.
% p_BestCoalitions                    - Computes  the set of most effective coalitions.
% p_clp_kernel                        - Computes a kernel point using the CLP solver.
% p_cls_kernel                        - Computes a kernel point using the CLS solver.
% p_CoalitionSolidarity               - Determines the coalition solidarity value.
% p_coeff_linearbasis                 - Determines the coefficients (dividends) of a linear basis from a TU game.
% p_Converse_RGP_Q                    - Checks if an imputation satisfies the CRGP.
% p_convex_gameQ                      - Checks on convexity.
% p_coreQ                             - Checks the non-emptiness of the core.
% p_cplex_kernel                      - Computes a kernel point using the CPLEX solver.
% p_cplex_prekernel                   - Computes a prekernel point using the CPLEX solver.
% p_COV_propertyQ                     - Verifies if the payoff x satisfies COV property.
% p_cvx_kernel                        - Computes a kernel point using the CVX solver. 
% p_cvx_prekernel                     - Computes a prekernel point using the CVX solver. 
% p_DM_Reduced_game                   - Computes all Davis-Maschler reduced games.
% p_excess                            - Computes the excesses.
% p_game_basis                        - Computes a game basis of the n-person TU-game space.
% p_game_space                        - Computes the game space which replicates a payoff as a pre-kernel element.
% p_game_space_red                    - Computes the game space which replicates a payoff as a pre-kernel element.
% p_Gap                               - Determines the gap function.
% p_genUnionStable                    - Creates a union stable system.
% p_getgame                           - Creates a Tu-game from the unanimity coordinates. 
% p_glpk_kernel                       - Computes a kernel point using the GLPK solver. 
% p_glpk_prekernel                    - Computes a prekernel point using the GLPK solver. 
% p_gurobi_kernel                     - Computes a kernel point using the GUROBI solver. 
% p_gurobi_prekernel                  - Computes a prekernel point using the GUROBI solver. 
% p_harsanyi_dividends                - Determines the the unanimity coordinates.
% p_HMS_ImputSavingReducedGame        - Computes from (v,x) all Hart/Mas-Colell ISR games.
% p_HMS_Reduced_game                  - Creates all Hart/Mas-Colell reduced games. 
% p_holler                            - Computes the Holler index.
% p_homogeneous_representationQ       - Checks if the weighted majority game possesses a homogeneous representation.
% p_hsl_prekernel                     - Computes a prekernel point using HSL solvers. 
% p_ImputSavingReducedGame            - Computes from (v,x) all imputation saving reduced games.
% p_ipopt_kernel                      - Computes a kernel point using the IPOPT solver.
% p_ipopt_prekernel                   - Computes a prekernel point using the IPOPT solver.
% p_ISRG_propertyQ                    - Checks whether an imputation x satisfies the ISR game property.
% p_k_Converse_RGP_Q                  - Checks if an imputation satisfies the k-CRGP.
% p_k_convexQ                         - Checks k-convexity of the Tu-game.
% p_k_cover                           - Determines from the Tu-game the corresponding k-game.
% p_Kernel                            - Computes a kernel point using the optimization toolbox.
% p_k_Reconfirmation_propertyQ        - Checks the k-RCP.
% p_k_Reduced_game_propertyQ          - Checks the k-RGP.
% p_k_StrConverse_RGP_Q               - Checks the strong k-CRGP. 
% p_linear_basis                      - Determines the linear basis of the n-person TU game space.
% p_lin_prekernel                     - Computes a prekernel point using optimization toolbox.
% p_LS_Nucl                           - Computes the least square nucleolus of a game.
% p_LS_PreNucl                        - Computes the least square pre-nucleolus of a game.
% p_minimal_winning                   - Computes the minimal winning coalitions.
% p_monotone_gameQ                    - Checks monotonicity of the Tu-game.
% p_msk_kernel                        - Computes a kernel point using the MOSEK solver.
% p_msk_prekernel                     - Computes a prekernel point using the MOSEK solver.
% p_MyersonValue                      - Computes the Myerson value of a Tu game.
% p_oases_kernel                      - Computes a kernel point using the OASES solver.
% p_oases_prekernel                   - Computes a prekernel point using the OASES solver.
% p_ols_prekernel                     - Computes a prekernel point using optimization toolbox.
% p_OwenValue                         - Computes the Owen value.
% p_PermutationGame                   - Computes from an assignment matrix the permutation game.
% p_PreKernel                         - Computes a pre-kernel element.
% p_PrekernelQ                        - Checks if an imputation is a pre-kernel point.
% p_pure_overhead                     - Creates the matrix of pure overhead games.
% p_qpBB_kernel                       - Computes a kernel point using the QPBB solver.
% p_qpc_kernel                        - Computes a kernel point using the QPC solver.
% p_qpc_prekernel                     - Computes a prekernel point using the QPC solver.
% p_qrg_prekernel                     - Computes a prekernel point using qrginv instead of pinv.
% p_reasonable_outcome                - Determines the reasonable outcome.
% p_Reconfirmation_propertyQ          - Checks the RCP.
% p_RedGame                           - Creates a Davis-Maschler reduced game.
% p_Reduced_game_propertyQ            - Checks the RGP.
% p_replicate_prk                     - Replicates a pre-kernel solution as a pre-kernel of a game space. 
% p_replicate_Shapley                 - Replicates the Shapley value for a game space. 
% p_select_starting_pt                - Selects a starting point for the pre-kernel computation.
% p_semi_convexQ                      - Checks semi-convexity.
% p_ShapleyValue                      - Computes the Shapley value (potential).
% p_ShapleyValueLB                    - Computes the Shapley value from the linear basis.
% p_ShapleyValueM                     - Computes the Shapley value while relying on all marginal contributions.
% p_SolidarityShapleyValue            - Determines the solidarity Shapley value. 
% p_SolidarityValue                   - Determines the solidarity value. 
% p_StrConverse_RGP_Q                 - Checks the strong RGP.
% p_SubSets                           - Creates all subsets of super set.
% p_super_additiveQ                   - Checks the Tu-game on super additivity.
% p_TauValue                          - Computes the Tau value.
% p_unanimity_games                   - Computes the unanimity coordinates.
% p_union_stableQ                     - Checks whether a system is union stable.
% p_zero_monotonicQ                   - Checks zero monotonicity.
%
%
% Class Objects
% -------------
% p_TuProp                            - subclass object of TuGame (game properties).
% p_TuSol                             - subclass object of TuGame (game solutions).
% p_TuVal                             - subclass object of TuGame (fairness and related values).
% p_TuRep                             - subclass object of p_TuSol (prk replication).
% p_TuShRep                           - subclass object of p_TuSol (Shapley value replication).
% p_TuCons                            - subclass object of p_TuSol (consistency).
% p_TuKcons                           - subclass object of p_TuSol (generalized consistency).
%
%
% tools: Sed File
%---------------
% sed_core                            - Converts cdd file format into Matlab format.
%
%
%  -----------------------------------------------------
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/11/2012        0.3             hme
%   09/25/2013        0.4             hme
%   11/03/2014        0.5             hme
%   03/06/2015        0.6             hme
%   04/06/2015        0.7             hme    
%   21/12/2015        0.8             hme        
%   24/04/2018        1.0             hme
%
 
end
