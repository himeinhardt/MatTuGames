% ChangeLog file for the MatTuGames Toolbox
%
%
%  Author:        Holger I. Meinhardt 
%  E-Mail:        holger.meinhardt@partner.kit.edu
%  Institution:   Karlsruhe Institute of Technology (KIT)
%
%
%
% Requirements
% ------------
% This release of *MatTuGames* was developed and tested using *Matlab
% R2024b* and earlier releases. A set of functions use the *Optimization Toolbox*
% and the *cdd-library* by *Komei Fukuda*, which can be found at the URL:
%
% [CDD](http://www.inf.ethz.ch/personal/fukudak/cdd_home/)
%
% as well as the Matlab interface to the cdd solver *CDDMEX*: 
%
% [CDDMEX](http://control.ee.ethz.ch/~hybrid/cdd.php.)
%
% Alternatively, in order to get even full scope of operation of the graphical features, one can also install the *MPT3* toolbox that can be downloaded from 
%
% [MPT3](http://control.ee.ethz.ch/~mpt/3/Main/Installation)
%
% which ships with *CDDMEX*.  We strongly recommend the user to apply
% the *MPT3 toolbox*, in particular of using the graphical features of
% our toolbox.
% 
% For the computation of the pre-kernel and related solutions the *SuiteSparse* for *Matlab* is 
% recommended that can be got from the URL
%
% [SuitSpare](https://github.com/DrTimothyAldenDavis/SuiteSparse)
%
% If you do not want to use *SuiteSparse*, then replace the function `qr_dec` by `pinv` in all functions
% for the pre-kernel and related solution. The same argument applies for the function `qrginv`.
%
% To run the toolbox even in parallel mode, *Matlab's Parallel Computing Toolbox* is needed.
%
% For connecting the *Mathematica* Package *TuGames*, the *Mathematica Symbolic Toolbox* is required, which can be found under
% the URL:
%
% [Mathematica Symbolic Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/6044-mathematica-symbolic-toolbox-for-matlab-version-2-0)
%
% whereas *TuGames* Version 3.1.4 can be downloaded from the URL:
%
% [TuGames](https://github.com/himeinhardt/TuGames)
%
% We recommend a custom installation with Paclet, which can be found at
%
% [Paclet](https://github.com/himeinhardt/TuGames/releases)
%
% The *MatTuGames* toolbox should work with all platforms.
%
% Moreover, the toolbox works also with the game theory toolbox written by *Jean Derks*, which can be downloaded from:
% 
% [Derks](http://www.personeel.unimaas.nl/Jean-Derks/downlds/CGinMatlab20050731.zip)
%
% This toolbox can be used to compute the pre-nucleolus up to 10-persons, if one has no license of Matlab's optimization toolbox. 
% 
% Finally, the toolbox offers interfaces to access the solvers of CVX, CPLEX, GLPK, GUROBI, HSL, IPOPT, MOSEK, and OASES. 
%
% To summarize, apart from the mentioned software, the toolbox requires the following MATLAB toolboxes:
%
% *MATLAB Parallel Server*,
% *Optimization Toolbox*,
% *Parallel Computing Toolbox*,
% *Signal Processing Toolbox*,
% *Statistics and Machine Learning Toolbox*,
% *Symbolic Math Toolbox*
%
% to get full functionality in serial as well as in parallel. 
%
%
% Installation
% ------------
% To install the MatTuGames Toolbox, unzip the zip-file mat_tugV1d0.zip, 
% and place the folder containing the functions on a local hard drive or 
% a network drive accessible to your computer. In the next step rename 
% the folder mat_tugV1d0 to mat_tug before including the folder location 
% in the MATLAB path. To set the MATLAB path, start MATLAB and then
% select the File/Set Path menu item. Then select Add Folder. Use the 
% navigation window to select the folder containing the functions. Click 
% OK and then click Save. The functions will then be ready for use within
% MATLAB. For additional installation instructions see the ReadMe.pdf
% file in the doc folder.
%
%
% Copyright
% ---------
% Copyright (c) 2012-25, Holger I. Meinhardt
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the author's institution nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%                                                 
% Disclaimer of Warranty
% ----------------------
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%
% The MatTuGames Toolbox was not created by The MathWorks Inc., nor is it 
% supported by the MathWorks.
%
% About Version 1.9.2
% -----------------
% Release Date: 11-Feb-2025
%
% 
% SUMMARY
% New in Version 1.9.2 is the implementation of a set of functions to get from an airport cost allocation 
% problem some capital cost distributions based on well-established game solutions. Moreover, a set of 
% functions is conceived to study Lorenz and Anti-Lorenz domination. Furthermore, some routines are added 
% to compute derived Banzhaf indices. In addition, some tools have been introduced to handle an MC-nets 
% representation. We have also rewritten some m-functions to cope with the depreciation of the MATLAB 
% interface for CPLEX. Finally, some minor bugs were fixed.
%
%
% NEW FUNCTIONS
%
% AUX: Some auxiliary files
%   
%   * FrameToImage2                              - Converts a frame to an image using MYAA.
%   * SaveFrames                                 - Saves the frames of a movie to different file formats.
%   * ToSimplex                                  - Projects data from 3d to 2d.
%   * toSymbols                                  - Converts a vector (vec, n) into digits. 
%
%
%  SERIAL COMPUTING
%
%   * ANucAirportProb                           - Computes the anti-nucleolus from an airport capital cost problem.
%   * AllMarginalContributionsResToS            - Computes all marginal worth vectors of a TU-game v restricted to coalition S.
%   * Anti_CPCore                               - Computes the closest point of the anti-core to x.
%   * Anti_LorenzDom                            - Checks if x anti-Lorenz dominates y in game v.
%   * Anti_LorenzMinACoreQ                      - Checks if x is Anti_Lorenz minimal in the anti-core of game v, i.e., x is in the Anti_Lorenz set.
%   * Anti_LorenzSol                            - Determines the anti-Lorenz solution of game v.
%   * BanzhafColeman                            - Computes the Banzhaf/Coleman index of a simple game sv (normalized Banzhaf value by 2^n-1).
%   * BanzhafOwenValue                          - Computes the Banzhaf-Owen value w.r.t. a priori unions cs.
%   * BanzhafPenrose                            - Computes the Banzhaf/Penrose and Banzhaf/Coleman index of a simple game sv.
%   * CanonicalOrder                            - Orders a set of coalitions by a canonical order.
%   * CddAnti_WeberSetPlot                      - Plots the anti-Weber set of game v.
%   * CddAnti_WeberSetSimplexPlot               - Plots the anti-Weber set of game v projected to the simplex.
%   * ChebyCenterCore                           - Computes the Cheby Center of the core of game v using MPT3.
%   * ColemanOwenValue                          - Computes the Coleman-Owen value w.r.t. a priori unions cs.
%   * ComposeMarkets                            - Composite from a parameter set defining number of players involved in each separated market situation a market game (merged markets).
%   * Composition                               - Composite at least 2 TU games up to 15 to a new extended game.
%   * Cost_Nucl                                 - Computes the nucleolus of a cost game v using the optimization toolbox.
%   * Cost_PreNucl                              - Computes the pre-nucleolus of a cost game v using the optimization toolbox.
%   * DecomposeInPositiveGames                  - Decomposes a TU game v into the difference of two positive games (convex games).
%   * EqDistDividends                           - Computes the equally distributed dividends of coalitions to which players belong, i.e., Shapley-value, of a TU-game v.
%   * GetAirPortProb                            - Computes a pseudo-random airport cost allocation problem of size m+1.
%   * GetMCNetRules                             - Transforms a cell array of the MC-nets representation of a TU game into a structure array.
%   * GetMarketGame                             - Determines from a random generated game the corresponding market game.
%   * GetPlayersCharacter                       - Determines from the set of players of the weighted majority game the characters Step, Sum and Null Player.
%   * GetProbDist                               - Tries to find for the game v the corresponding probability distribution over the players' orderings w.r.t. pre-imputation x.
%   * HarsanyiValue                             - Computes a Harsanyi-value of a TU-game v.
%   * LexMaxMinWin                              - Determines the lexicographical maximal coalition from the clm represented as a cell, matrix, or array of integers.
%   * LexOrder                                  - Sorts a set of coalitions represented as a cell, matrix, or array of integers into the corresponding lexicographical order.
%   * LorenzMaxCoreQ                            - Checks if x is Lorenz maximal in the core of game v, i.e., x is in the Lorenz set.
%   * MarketGameQ                               - Checks whether the game v is a market game.
%   * NucAirportProb                            - Computes the nucleolus from an airport capital cost problem.
%   * PartitionGameQ                            - Verifies if (th,w_vec) induces a partition game.
%   * PreNucl_mod                               - Computes the pre-nucleolus of game v using the optimization toolbox and a spanning method.
%   * PreNucl_mod2                              - Computes the pre-nucleolus of game v using the optimization toolbox and an alternative spanning method.
%   * PreNucl_mod3                              - Computes the pre-nucleolus of game v using the optimization toolbox and an alternative spanning method.
%   * PreNucl_mod4                              - Computes the pre-nucleolus of game v using the optimization toolbox and an alternative spanning method.
%   * PrkEqPnQ                                  - Checks whether the pre-kernel of a game v is a singleton.
%   * ProportionalGame                          - Computes from (q,co) a proportional game.
%   * ReverseMCNetsRep                          - Reverse the MC-nets (concise) representation of a cooperative game into a full representation.
%   * ShapAirPortMod2                           - Computes the Shapley value from an airport capital cost problem (modified).
%   * ShapleyAirportProb                        - Computes the Shapley value from an airport capital cost problem.
%   * SplitSimpleGame                           - Splits a simple game with winning players into sub-games with and without winning players.
%   * TauValAirportProb                         - Computes the tau value from an airport capital cost problem.
%   * WedgeProdGame                             - Computes from an MC-nets (concise) representation a cooperative production game (wedge production game).
%   * airport_costgame                          - Computes from an airport problem the associated airport cost game.
%   * balancedCoverQ                            - Checks whether the characteristic function v is equal to a balanced cover.
%   * belongToAllSubCores                       - Checks whether all projections of an imputation x are a member of an associated core of a subgame.
%   * belongToWeberSetQ                         - Checks if the imputation x belongs to the Weber set using MPT3.
%   * bidding_collusion_game                    - Determines a bidding collusion game from a bid vector.
%   * bs_PreKernel                              - Computes from (v,x) a pre-kernel element using MATLAB's backslash instead of pinv.
%   * clp_weightedKernel                        - Computes a weighted kernel point using the CLP solver.
%   * cls_weightedKernel                        - Computes a weighted kernel point using the CLS solver.
%   * complementary_basis                       - Computes the complementary basis of the n-person TU game space.
%   * complementary_dividends                   - Computes the complementary dividends.
%   * cp_kernel                                 - Computes a kernel element using the cone programming solver of MATLAB's Optimization Toolbox.
%   * cp_prekernel                              - Computes a pre-kernel element using the cone programming solver of MATLAB's Optimization Toolbox.
%   * cs_Banzhaf                                - Computes the Banzhaf value w.r.t. a communication situation.
%   * cvx_weightedKernel                        - Computes a weighted kernel point using the CVX solver.
%   * cvx_weightedPreKernel                     - Computes a weighted prekernel point using the CVX solver.
%   * diffOperator                              - Computes the difference operator of the game w.r.t. coalition T.
%   * essentialSet                              - Computes the set of essential coalitions.
%   * intersection_basis                        - Computes the intersection basis of the n-person TU game space.
%   * landlord                                  - Computes a production game arising from l-landlords and t-tenants.
%   * linprog_game                              - Computes from a production matrix (A,H,mB,mD,p) a linear programming game.
%   * ma57_prekernel                            - Computes from (v,x) a pre-kernel element using HSL MA57.
%   * ma86_prekernel                            - Computes from (v,x) a pre-kernel element using HSL MA86.
%   * ma87_prekernel                            - Computes from (v,x) a pre-kernel element using HSL MA87.
%   * ma97_prekernel                            - Computes from (v,x) a pre-kernel element using HSL MA97.
%   * msk_prenucl_mod                           - Computes the pre-nucleolus of game v using mosekmex and a spanning method.
%   * msk_prenucl_mod2                          - Computes the pre-nucleolus of game v using mosekmex and an alternative spanning method.
%   * msk_prenucl_mod3                          - Computes the pre-nucleolus of game v using mosekmex and an alternative spanning method.
%   * msk_prenucl_mod4                          - Computes the pre-nucleolus of game v using mosekmex and an alternative spanning method.
%   * nucl_formula                              - Computes the nucleolus of a three-person super-additive game v from a formula.
%   * one_concaveQ                              - Checks whether the game v is 1-concave.
%   * one_convexQ                               - Checks whether the game v is 1-convex.
%   * positive_gameQ                            - Returns true (1) if all Harsanyi dividends are non-negative, otherwise false (zero).
%   * probability_game                          - Computes from a positive vector x the corresponding probability game.
%   * rootedQ                                   - Checks whether the game v is rooted.
%   * secDiffOperatorQ                          - Checks whether the second order differences are all positive which indicates convexity.
%   * strictconvex_gameQ                        - Returns 1 whenever the game v is strictly convex.
%   * totallyBalancedCoverQ                     - Checks whether the characteristic function v is equal to a totally balanced cover.
%
%
%  PARALLEL COMPUTING
%
%   * p_AntiBalancedSetQ                         - Verifies whether the set of induced coalitions is an anti balanced collection. Checking Kohlberg's criterion.
%   * p_BanzhafColeman                           - Computes the Banzhaf/Coleman index of a simple game sv using MATLAB's PCT.
%   * p_BanzhafOwenValue                         - Computes the Banzhaf-Owen value w.r.t. a priori unions cs using Matlab's PCT.
%   * p_BanzhafPenrose                           - Computes the Banzhaf/Penrose and Banzhaf/Coleman index of a simple game sv using MATLAB's PCT.
%   * p_ColemanOwenValue                         - Computes the Coleman-Owen value w.r.t. a priori unions cs using Matlab's PCT. Is not normalized to one.
%   * p_DM_TwoReduced_game                       - Computes from (v,x) all single and two-person reduced games on S at x of game v using MATLAB's PCT.
%   * p_GetMarketGame                            - Determines from a random generated game the corresponding market game.
%   * p_HMS_TwoReduced_game                      - Computes from (v,x) all Hart/Mas-Colell singleton and two-person reduced games on S at x of game v.
%   * p_MarketGameQ                              - Checks whether the game v is a market game.
%   * p_ShapleyValueML                           - Computes the Shapley-value of a TU-game v using the multilinear extension.
%   * p_WeakReduced_game_propertyQ               - Checks whether an imputation x satisfies the weak reduced game property (consistency)
%   * p_bankruptcy_game                          - Computes for a bankruptcy situation (E,d_vec) the corresponding bankruptcy.
%   * p_belongToAllSubCores                      - Checks whether all projections of an imputation x are a member of an associated core of a subgame.
%   * p_bs_PreKernel                             - Computes from (v,x) a pre-kernel element using Matlab's PCT and backslash operation.
%   * p_concave_gameQ                            - Returns 1 whenever the game v is concave.
%   * p_cp_kernel                                - Computes from (v,x) a kernel element using PCT and the cone programming solver of MATLAB's Optimization Toolbox.
%   * p_cp_prekernel                             - Computes from (v,x) a pre-kernel element using PCT and the cone programming solver of MATLAB's Optimization Toolbox.
%   * p_landlord                                 - Computes a production game arising from l-landlords and t-tenants.
%   * p_one_concaveQ                             - Checks whether the game v is 1-concave.
%   * p_one_convexQ                              - Checks whether the game v is 1-convex.
%   * p_secDiffOperatorQ                         - Checks whether the second order differences are all positive which indicates convexity.
%   * p_totallyBalancedCoverQ                    - Checks whether the characteristic function v is equal to a totally balanced cover.
%   * p_weighted_EssentialSet                    - Computes the set of weighted essential coalitions using Matlab's PCT.
%   * p_weighted_majority                        - Computes from the threshold th and the weights w_vec a weighted majority game using Matlab's PCT.
%   * p_winning_coalitions                       - Computes from a pre-defined set of winning coalitions (e.g. minimal winning coalitions) the whole set of winning coalitions. 
%
%
% BUG FIXES
%  Incorrect computation of PartitionPlySet, center_solution, and minimal_representation. The latter is replaced by the new m-function min_homogrep. The former was overhauled. However, it is recommended to replace the former by GetPlayersCharacter.
%
%
% MODIFIED FUNCTIONS
% SERIAL
%  airport_costgame,airport_game,AllMarginalContributions,AllMarginalContributionsResToS,AlmostAverageConcaveQ,AlmostAverageConvexQ,
%  AlmostConcave_gameQ,AlmostConvex_gameQ,alphaVector,Anti_B0_balancedCollectionQ,Anti_balancedCollectionQ,anti_coreQ,Anti_CPCore,Anti_Kernel,
%  Anti_LorenzDom,Anti_LorenzMinACoreQ,Anti_LorenzSol,Anti_Modiclus,Anti_modiclusQ,Anti_ModPreKernel,Anti_ModPreKernel2,Anti_Nucl,Anti_Nucl_llp,
%  Anti_PModPreKernel,Anti_PreKernel,Anti_PreNucl,Anti_PreNucl_llp,Anti_Weak_balancedCollectionQ,ANucAirportProb,apu_SolidarityValue,
%  B0_balancedCollectionQ,B0_balancedQ,balancedCollectionQ,balancedCoverQ,balancedQ,belongToCoreQ,BestCoalitions,bs_PreKernel,
%  CddCorePlot,CddCoreQ,CddCoreSimplexPlot,CddNucl,CddNucl_llp,CddPreKernel,CddPrenucl,CddPrenucl_llp,CddWeightedPreKernel,center_solution,
%  cfr_Anti_balancedCollectionQ,cfr_Anti_Kernel,cfr_Anti_Nucl,cfr_Anti_PreKernel,cfr_Anti_PreNucl,cfr_Anti_Weak_balancedCollectionQ,cfr_balancedCollectionQ,
%  cfr_Kernel,cfr_nucl,cfr_PreKernel,cfr_PreNucl,cfr_weak_balancedCollectionQ,clp_kernel,clp_weightedKernel,cls_kernel,cls_weightedKernel,
%  Converse_RGP_Q,coreQ,Cost_Nucl,Cost_PreNucl,COV_propertyQ,CPCore,cp_kernel,cplex_kernel,cplex_linear_production,cplex_linprog_game,cplex_nucl,cplex_nucl_llp,
%  cplex_prekernel,cplex_prenucl,cplex_prenucl_llp,cplex_prenucl_mod,cplex_prenucl_mod2,cplex_prenucl_mod3,cplex_prenucl_mod4,cplex_weightedKernel,cplex_weightedPreKernel,
%  cp_prekernel,cs_Anti_balancedCollectionQ,cs_Anti_Kernel,cs_Anti_Nucl,cs_Anti_PreKernel,cs_Anti_PreNucl,cs_Anti_Weak_balancedCollectionQ,
%  cs_balancedCollectionQ,cs_Banzhaf,cs_Kernel,cs_nucl,cs_PreKernel,cs_PreNucl,cs_Weak_balancedCollectionQ,cvx_kernel,cvx_prekernel,cvx_weightedKernel,
%  cvx_weightedPreKernel,dc_PreKernel,Derived_game_propertyQ,diffOperator,Dual_Cover_propertyQ,Dual_Floor_propertyQ,DuallyEssentialSet,
%  EC_DGP_Q,EC_propertyQ,EC_RGP_Q,EqClasses,essentialSet,exact_game,FindPartition,flow_game,Gap,glpk_kernel,glpk_nucl,glpk_nucl_llp,glpk_prekernel,glpk_prenucl,
%  glpk_prenucl_llp,glpk_weightedKernel,glpk_weightedPreKernel,gurobi_kernel,gurobi_nucl,gurobi_nucl_llp,gurobi_prekernel,gurobi_prenucl,gurobi_prenucl_llp,
%  gurobi_weightedKernel,gurobi_weightedPreKernel,HMS_AntiReduced_game,HMS_ImputSavingReducedGame,HMS_RedGame,HMS_Reduced_game,homogeneous_representationQ,
%  hsl_prekernel,hsl_weightedPreKernel,intersection_basis,ipopt_kernel,ipopt_prekernel,ipopt_weightedKernel,ipopt_weightedPreKernel,ItChebyCenter,
%  k_anticover,Kernel,landlord,LeastCore,Ledcons_propertyQ,LEDcons_propertyQ,linear_production,lin_prekernel,linprog_game,lin_weightedPreKernel,market_game,
%  MarketGameQ,MaxConsistencyQ,min_homogrep,minNoBlockPayoff,Modiclus,modiclusQ,ModLeastCore,ModPreKernel,ModPreKernel2,msk_AntiPreKernel,msk_kernel,
%  msk_nucl,msk_nucl_llp,msk_prekernel,msk_PreKernel,msk_prenucl,msk_prenucl_llp, msk_weightedKernel,msk_weightedPreKernel,NetworkCenterSol,NetworkMinimalRep,
%  nucl_formula,nucl_llp,oases_kernel,oases_prekernel,oases_weightedKernel,oases_weightedPreKernel,ols_prekernel,ols_weightedPreKernel,one_concaveQ,one_convexQ,
%  OwenValue,PartitionPlySet,PermutationGame,PlayersCharacter,PlyEqFirstStep,PModPreKernel,PModPreKernel2,positive_gameQ,positiveKernel,positivePreKernel,
%  Potential,PreKernel,PrekernelQ,PreNucl,PreNucl2,PreNucl_llp,PreNucl_mod,PreNucl_mod2,PreNucl_mod3,PreNucl_mod4,PrkEqPnQ,product_game,PropNucl,ProportionalGame,
%  PropPreNucl,qpBB_kernel,qpBB_weightedKernel,qpc_kernel,qpc_prekernel,qpc_weightedKernel,qpc_weightedPreKernel,qrg_prekernel,qrg_weightedPreKernel,
%  Reconfirmation_propertyQ,RedGame,Reduced_game_propertyQ,satisfaction,separating_collectionQ,ShapleyValue,ShapleyValueM,
%  ShapleyValueML,sortsets,SplitSimpleGame,StandardSolution,TauValue,totallyBalancedCoverQ,totallyBalancedQ,Weak_balancedCollectionQ,WeakReduced_game_propertyQ,
%  WedgeProdGame,weightedAnti_B0_balancedCollectionQ,weightedAntiBalancedCollectionQ,weightedAnti_Kernel,weightedAnti_Nucl,weightedAnti_Nucl_llp,
%  weightedAnti_PreKernel,weightedAnti_PreNucl,weightedAnti_PreNucl_llp,weightedB0_balancedCollectionQ,weightedBalancedCollectionQ,weighted_cfr_balancedCollectionQ,
%  weighted_cfr_Kernel,weighted_cfr_nucl,weighted_cfr_PreKernel,weighted_cfr_PreNucl,weighted_cfr_weak_balancedCollectionQ,WeightedComposition,
%  weightedCsB0_balancedCollectionQ,weightedCsBalancedCollectionQ,weightedCsKernel,weightedCsNucl,weightedCsPreKernel,weightedCsPreNucl,weightedGraphGame,
%  weightedKernel,weightedNucl,weightedNucl_llp,weightedPreKernel,weightedPreNucl,weightedPreNucl_llp,weighted_truncated,WSysKernelQ,WSysNuclQ,WSysPreKernelQ,WSysPreNuclQ
%  
%
%
% PARALLEL
%  p_AllMarginalContributions,p_Anti_ModPreKernel,p_Anti_PModPreKernel,p_balancedSetQ,p_clp_kernel,p_clp_weightedKernel,p_cls_kernel,p_cls_weightedKernel,p_Converse_RGP_Q,
%  p_COV_propertyQ,p_cplex_kernel,p_cplex_prekernel,p_cplex_weightedKernel,p_cplex_weightedPreKernel,p_cs_CoreQ,p_cvx_kernel,p_cvx_prekernel,p_cvx_weightedKernel,
%  p_cvx_weightedPreKernel,p_dc_PreKernel,p_Derived_game_propertyQ,p_DualEssentialSet,p_DuallyEssentialSet,p_EC_DGP_Q,p_EC_RGP_Q,p_essentialSet,p_exact_game,p_flow_game,
%  p_GenGap,p_glpk_kernel,p_glpk_prekernel,p_glpk_weightedKernel,p_glpk_weightedPreKernel,p_gurobi_kernel,p_gurobi_prekernel,p_gurobi_weightedKernel,p_gurobi_weightedPreKernel,
%  p_HMS_AntiReduced_game,p_HMS_ImputSavingReducedGame,p_HMS_Reduced_game,p_hsl_prekernel,p_hsl_weightedPreKernel,p_ipopt_kernel,p_ipopt_prekernel,p_ipopt_weightedKernel,
%  p_ipopt_weightedPreKernel,p_Kernel,p_Ledcons_propertyQ,p_LED_RGPvsDGP,p_lin_prekernel,p_lin_weightedPreKernel,p_ModPreKernel,p_msk_kernel,p_msk_prekernel,p_msk_weightedKernel,
%  p_msk_weightedPreKernel,p_MyersonValue,p_oases_kernel,p_oases_prekernel,p_oases_weightedKernel,p_oases_weightedPreKernel,p_ols_prekernel,p_ols_weightedPreKernel,p_PermutationGame,
%  p_PModPreKernel,p_positiveKernel,p_positivePreKernel,p_PreKernel,p_qpBB_kernel,p_qpBB_weightedKernel,p_qpc_kernel,p_qpc_prekernel,p_qpc_weightedKernel,p_qpc_weightedPreKernel,
%  p_qrg_prekernel,p_qrg_weightedPreKernel,p_Reconfirmation_propertyQ,p_Reduced_game_propertyQ,p_replicate_prk,p_ShapleyValue,p_ShapleyValueM,p_weightedAnti_PreKernel,
%  p_weighted_DualEssentialSet,p_weighted_DuallyEssentialSet,p_weightedKernel,p_weightedPreKernel
%
%
%
%
% About Version 1.9
% -----------------
% Release Date: 25-Jun-2021
%
% 
% SUMMARY
% New in Version 1.9 is the implementation of a set of functions to compute 
% the public good index w.r.t. a network. Moreover, functions to compute 
% chi-, Gately-, or Equal Collective Gains value have been added. In addition, 
% a set of routines to check whether a game is almost-convex almost-concave etc. 
% were implemented. Due to some miscalculations, the function to compute a minimum cost 
% spanning game has been completely overhauled, and some new routines to study 
% this game class were enclosed. Finally, some minor bugs were fixed.  
%
%
% NEW FUNCTIONS
%
% AUX: Some auxiliary files
%   * NONE
%
%  SERIAL COMPUTING
%
%   * AP_DummyPlayer_propertyQ                   - Checks if the solution x satisfies the AP-Dummy player property.
%   * AP_DummyPlayers                            - Returns the player who are AP-Dummy players.
%   * AP_NullPlayer_propertyQ                    - Checks if the solution x satisfies the AP-Null player property. 
%   * AP_NullPlayers                             - Returns the players who are AP-Null players. 
%   * A_DummyPlayer_propertyQ                    - Checks if the solution x satisfies the A-Dummy player property. 
%   * A_NullPlayer_propertyQ                     - Checks if the solution x satisfies the A-Null player property. 
%   * A_NullPlayers                              - Returns the players who are A-Null players. 
%   * AlmostConcave_gameQ                        - Returns true whenever the game v is almost concave.
%   * AlmostConvex_gameQ                         - Returns true whenever the game v is almost convex.
%   * Anti_ChiValue                              - Computes the anti-chi-value of a TU-game v. 
%   * Anti_Gap                                   - Computes the anti-gap function from game v. 
%   * Anti_GenGap                                - Computes the anti-generalized gap function from game v. 
%   * Anti_Monotonic_Cover                       - Computes the anti-monotonic cover of game v. 
%   * Anti_TauValue                              - Computes the anti-tau-value of a TU-game v. 
%   * Anti_Weak_balancedCollectionQ              - Checking reverse weak Kohlberg's criterion.             
%   * BaryCenter                                 - Computes the barycenter of the core.
%   * CddLinProg_game                            - Computes from a production matrix (A,H,mB,mD,p) a linear programming game using cddmex. 
%   * CddLinear_Production                       - Computes from a production problem (A,mB,p) a linear production game using cddmex. 
%   * ChiValue                                   - Computes the chi-value of a TU-game v. This is a generalized Tau value.
%   * CmpConsistencyQ                            - Checks whether an imputation x satisfies the complement consistency.
%   * CmpRedGame                                 - Computes from (v,x,S) a complement reduced game vS on S at x for game v.
%   * Complement_Reduced_game                    - Computes from (v,x) all complement reduced games on S at x of game v.
%   * Converse_CmpConsistencyQ                   - Checks whether an imputation x satisfies the converse complement consistency property.
%   * DecomposableQ                              - Checks whether the game v is decomposable w.r.t. the coalition structure cs. 
%   * DeeganPackel                               - Computes the Deegan-Packel index from the set of minimal winning coalitions. 
%   * DeeganPackel_SV                            - Computes the Deegan-Packel index from a simple game to construct the set of minimal winning coalitions. 
%   * DerivedCostMatrix                          - Computes from a cost matrix and a partition of the player set N the corresponding derived cost matrix. 
%   * DummyPlayer_propertyQ                      - Checks the dummy player property. 
%   * DummyPlayers                               - Returns the list of dummy players of game v. 
%   * EANSCValue                                 - Computes the Equal Allocation of Non-Separable Contribution/Cost Value.
%   * ECGValue                                   - Computes the Equal Collective Gains value of a TU-game v. 
%   * EPSDValue                                  - Computes the egalitarian proportional surplus division value of an individually positive TU-game.
%   * ExtShapleyValue                            - Computes the extended Shapley-value.
%   * GatelyValue				 - Computes the Gately point of an essential game v.
%   * GenGap                                     - Computes the generalized gap function from game v.
%   * InessGame                                  - Computes the inessential game from a payoff vector.
%   * IrredAntiCore                              - Computes from a cost matrix the corresponding extreme points of the irreducible anti-core of the associated m.c.s.t. game. 
%   * IrredCostMatrix                            - Computes from a cost matrix and a solution tree the irreducible cost matrix. 
%   * Johnston                                   - Computes the Johnston power index from the set of winning coalitions. 
%   * MIMC                                       - Computes the vector of minimum increase in players marginal contribution when they leave the grand coalition.
%   * MTRCostMatrix                              - Computes from a cost matrix and a solution tree the cost matrix of a minimal spanning tree. 
%   * MaxConsistencyQ                            - Checks whether an imputation x satisfies maximal consistency.
%   * ModDeeganPackel                            - Computes the modified Deegan-Packel index from the set of winning coalitions. 
%   * ModDeeganPackel_SV                         - Computes the Deegan-Packel index from a simple game to construct the set of minimal winning coalitions. 
%   * ModPGI                                     - Computes the modified public good index from the set of minimal winning coalitions.
%   * ModPGI_SV                                  - Computes the modified public good index from a simple game to determine the set of minimal winning coalitions. 
%   * NetworkBanzhaf                             - Computes the network Banzhaf power index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * NetworkDeeganPackel                        - Computes the network Deegan-Packel index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * NetworkJohnston                            - Computes the network Johnston power index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * NetworkMajorityGame                        - Computes a network majority TU game (simple game).
%   * NetworkMinimalRep                          - Computes from the set of edges, threshold th and the weights w_vec the minimal homogeneous representation of a homogeneous network weighted majority game. 
%   * NetworkModDeeganPackel                     - Computes the network modified Deegan-Packel index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * NetworkModPGI                              - Computes the modified network public good index from the set of minimal winning coalitions of a network.
%   * NetworkPGI                                 - Computes the network public good index from the set of minimal winning coalitions of a network.
%   * NetworkShapleyShubik                       - Computes the network Shapley-Shubik power index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * NullPlayers                                - Returns the list of null players of game v. 
%   * PDValue                                    - Computes the proportional division value of an individually positive TU-game.
%   * PGI                                        - Computes the public good index from the set of minimal winning coalitions.
%   * PGI_SV                                     - Computes the public good index from a simple game to determine the set of minimal winning coalitions. 
%   * PS_GameBasis                               - Computes the basis for the class of PS games.
%   * PartitionPlySet                            - Partitions the set of players of the weighted majority game into character Sum, Step, and Null-Player. 
%   * PlayersCharacter                           - Partitions the set of players of the weighted majority game into the character Sum, Step, and Null-Player.
%   * PlotCostGraph                              - Plots from a cost matrix the associated cost spanning graph. 
%   * SolidarityPGI                              - Computes the solidarity Holler index w.r.t. a priori unions cs. 
%   * SubCoalitions                              - Computes the power set (subsets) from an array.
%   * SuperAddSolQ                               - Checks if the vector x is an element of a super additive solution of the game v.
%   * SuperSets                                  - Computes the super-sets of set S. 
%   * Tk_lin_basis                               - Computes the (T,k) intermediate linear basis of the n-person TU game space.
%   * anti_partition                             - Computes from a partition its anti partition.
%   * apu_PGI                                    - Computes the Holler index w.r.t. a priori unions cs. 
%   * average_concaveQ                           - Returns true whenever the game v is average-concave.
%   * bint_AssignmentGame                        - Creates an assignment game (bintprog).
%   * contentment                                - Computes the contentment vector of game v w.r.t. x.
%   * complementary_game                         - Generates a producer and buyer game.
%   * compromiseAntiAdmissibleQ                  - Checks if the anti-core cover a TU game v is non-empty.
%   * contentment                                - Computes the contentment vector of game v w.r.t. x.
%   * disagreement                               - Computes the disagreement vector of game v.
%   * flow_game                                  - Computes from a flow problem a TU flow game using the optimization toolbox.
%   * flow_probMinCut                            - Computes from a flow problem a minimal cut.
%   * getCOV                                     - Computes from a sample of observations obs and for n-assets the covariance matrix V of a portfolio with indefinite risk-return relationship. 
%   * getCOV2                                    - Computes from a sample of observations obs and for n-assets the covariance matrix V of a portfolio with negative risk-return relationship.          
%   * getCOV3                                    - Computes from a sample of observations obs and for n-assets the covariance matrix V of a portfolio with positive risk-return relationship. 
%   * getMinimalWinning                          - Computes from a simple game the minimal winning coalitions.
%   * getPSgame                                  - Computes a PS game from the PS game basis.
%   * getSymCostMatrix                           - Computes a symmetric cost matrix from the cardinality of the player set and an upper bound value to specify the range from which the random number are drawn. 
%   * linear_production                          - Computes from a production problem (A,mB,p) a linear production game. 
%   * linprog_game                               - Computes from a production matrix (A,H,mB,mD,p) a linear programming game. 
%   * market2_game                               - Determines from two disjoint sets a market game.
%   * min_aspiration                             - Computes the minimum aspiration level of players of game v.
%   * min_epsshift                               - Computes for an almost-convex game the min epsilon shift to construct a convex game. 
%   * min_game                                   - Generates a minimum game.
%   * min_homogrep                               - Computes from the threshold th and the weights w_vec the minimal homogeneous representation of a homogeneous weighted majority game. 
%   * nullShapley                                - Determines a basis of the null space for the Shapley-value for n-persons.
%   * nullShapleyLB                              - Determines a counting basis of the null space for the Shapley-value for n-persons.
%   * oddeven_game                               - Assigns |S|-1 if S is odd and |S|+1 if S is even.
%   * product_game                               - Computes from a vector x the corresponding product game.
%   * production_game2                           - Creates an affine production game.
%   * ps_gameQ                                   - Checks whether a game is a PS game.
%   * psstar_gameQ                               - Checks whether a game is a PS* game.
%   * root_game                                  - Computes from game v its associated root game.
%   * semi_concaveQ                              - Checks semi-concavity.
%   * separating_collectionQ                     - Verifies if a collection is separating.
%   * superadditive_cover                        - Computes from game v its superadditive cover.
%   * surplus_game                               - Computes from a cost game c the corresponding surplus game v. 
%   * totallyAntiBalancedQ                       - Checks whether the anti-core of all subgames is non-empty. 
%   * tricameral_assembly                        - Computes from a set of parameters a simple game. 
%   * veto_rich_players                          - Returns a list of veto players for the TU-game v.
%   * winning_players                            - Computes from a pre-defined set of winning coalitions (e.g. minimal winning coalitions) the set of winning players. 
%
%
%  PARALLEL COMPUTING
%
%   * p_AP_DummyPlayer_propertyQ                 - Checks if the solution x satisfies the AP-Dummy player property.
%   * p_AP_DummyPlayers                          - Returns the player who are AP-Dummy players. 
%   * p_AP_NullPlayer_propertyQ                  - Checks if the solution x satisfies the AP-Null player property. 
%   * p_AP_NullPlayers                           - Returns the players who are AP-Null players. 
%   * p_A_DummyPlayer_propertyQ                  - Checks if the solution x satisfies the A-Dummy player property. 
%   * p_A_NullPlayer_propertyQ                   - Checks if the solution x satisfies the A-Null player property.
%   * p_A_NullPlayers                            - Returns the players who are A-Null players. 
%   * p_AlmostConcave_gameQ			 - Returns true whenever the game v is almost concave.
%   * p_AlmostConvex_gameQ			 - Returns true whenever the game v is almost convex.
%   * p_Anti_ChiValue                            - Computes the anti-chi-value of a TU-game v. 
%   * p_Anti_Gap                                 - Computes the anti-gap function from game v. 
%   * p_Anti_GenGap                              - Computes the anti-generalized gap function from game v. 
%   * p_Anti_TauValue                            - Computes the anti-tau-value of a TU-game v. 
%   * p_CddCoreMovie                             - Creates a movie w.r.t. the strong epsilon-cores.
%   * p_ChiValue				 - Computes the chi-value of a TU-game v.
%   * p_CmpConsistencyQ                          - Checks whether an imputation x satisfies the complement consistency.
%   * p_CmpRedGame				 - Computes from (v,x,S) a complement reduced game vS on S at x for game v.
%   * p_Complement_Reduced_game                  - Computes from (v,x) all complement reduced games on S at x of game v.
%   * p_Converse_CmpConsistencyQ		 - Checks whether an imputation x satisfies the converse complement consistency property.
%   * p_DeeganPackel                             - Computes the Deegan-Packel index from the set of minimal winning coalitions. 
%   * p_DeeganPackel_SV                          - Computes the Deegan-Packel index from a simple game to construct the set of minimal winning coalitions. 
%   * p_DummyPlayer_propertyQ                    - Verifies if x satisfies the dummy player property. 
%   * p_DummyPlayers                             - Returns the list of dummy players of game v. 
%   * p_ECGValue                                 - Computes the Equal Collective Gains value of a TU-game v. 
%   * p_GenGap                                   - Computes the generalized gap function from game v.
%   * p_GenGap                                   - Computes the generalized gap function from game v. 
%   * p_InessGame                                - Computes the inessential game from a payoff vector.
%   * p_Johnston                                 - Computes the Johnston power index from the set of winning coalitions. 
%   * p_MaxConsistencyQ                          - Checks whether an imputation x satisfies maximal consistency.
%   * p_ModDeeganPackel                          - Computes the modified Deegan-Packel index from the set of winning coalitions. 
%   * p_ModDeeganPackel_SV                       - Computes the Deegan-Packel index from a simple game to construct the set of minimal winning coalitions. 
%   * p_ModPGI                                   - Computes the modified public good index from the set of winning coalitions.
%   * p_ModPGI_SV                                - Computes the modified public good index from a simple game to determine the set of minimal winning coalitions. 
%   * p_NetworkBanzhaf                           - Computes the network Banzhaf power index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * p_NetworkDeeganPackel                      - Computes the network Deegan-Packel index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * p_NetworkJohnston                          - Computes the network Johnston power index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * p_NetworkMajorityGame			 - Computes from a network problem (E,c,th) a network majority TU game (simple game).
%   * p_NetworkModDeeganPackel                   - Computes the network modified Deegan-Packel index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * p_NetworkModPGI                            - Computes the network modified public good index from the set of winning coalitions of a network E while imposing a threshold.
%   * p_NetworkPGI				 - Computes the network public good index from the set of minimal winning coalitions of a network E while imposing a threshold.
%   * p_NetworkShapleyShubik                     - Computes the network Shapley-Shubik power index from the set of winning coalitions of a network E while imposing a threshold of th. 
%   * p_NullPlayer_propertyQ                     - Verifies if x satisfies the null player property. 
%   * p_NullPlayers                              - Returns the list of null players of game v. 
%   * p_PGI					 - Computes the public good index from the set of minimal winning coalitions.
%   * p_PGI_SV                                   - Computes the public good index from a simple game to determine the set of minimal winning coalitions. 
%   * p_UpperPayoff				 - Computes the utopia and minimum claim vector of game v.
%   * p_average_concaveQ			 - Returns true whenever the game v is average-concave.
%   * p_bint_AssignmentGame                      - Creates an assignment game (bintprog).
%   * p_disagreement                             - Computes the disagreement vector of game v.
%   * p_flow_game                                - Computes from a flow problem a TU flow game using the optimization toolbox.
%   * p_getMinimalWinning                        - Computes from a simple game the minimal winning coalitions.
%   * p_grMaxFlowGame                            - Computes from a flow problem a TU flow game.
%   * p_min_aspiration                           - Computes the minimum aspiration level of players of game v.
%   * p_nullShapley                              - Determines a basis of the null space for the Shapley-value for n-persons.
%   * p_proper_amount                            - Computes the largest amount players can contribute to a proper coalition.
%   * p_semi_concaveQ                            - Checks semi-concavity.
%   * p_smallest_amount                          - Computes the smallest amount vector of the game.
%   * p_superadditive_cover                      - Computes from game v its superadditive cover.
%   * p_tricameral_assembly                      - Computes from a set of parameters a simple game. 
%   * p_veto_rich_players			 - Returns a list of veto players for the TU-game v.
%
%
%
% BUG FIXES
%  Incorrect computation of LorenzDom due to unsorted vectors has been fixed. Fixing incorrect computation of mcst_game, and of the anti-nucleolus/anti-kernel functions.
%
%
% MODIFIED FUNCTIONS
% SERIAL
%  AllMarginalContributions, AntiCoreVertices, AntiImputationVertices, AntiUtopiaPayoff, Anti_B0_balancedCollectionQ, Anti_Converse_DGP_Q, 
%  Anti_Kernel, Anti_ModPreKernel, Anti_Nucl, Anti_Nucl_llp, Anti_PModPreKernel, Anti_PreKernel, Anti_PreNucl, Anti_PreNucl_llp, 
%  Anti_PropModPreKernel, Anti_kernelQ, B0_balancedQ, COV_propertyQ, CddAntiCorePlot, CddAntiCoreSimplexPlot, CddAntiCoreVertices, 
%  CddAntiImputationVertices, CddAntiNucl, CddAntiNucl_llp, CddAntiPrenucl, CddAntiPrenucl_llp, CddCoreQ, CddCoreVertices, CddLeastCore, 
%  CddModCoreVertices, CddModLeastCoreVertices, CddModiclus, CddNucl, CddNucl_llp, CddPreKernel, CddPrenucl, CddTotallyBalancedQ, 
%  ConstantSumQ, Converse_RGP_Q, CoreVertices, DCP_propertyQ, DFP_propertyQ, DM_AntiReduced_game, DM_Anti_Derived_game, DM_Derived_game, 
%  DM_TwoReduced_game, DRP_propertyQ, DerivedGame, Derived_game_propertyQ, DualCover, DualFloor, Dual_Cover_propertyQ, DuttaRay, 
%  ECCoverGame, ECFloorGame, EC_DGP_Q, EC_RGP_Q, EC_propertyQ, HMS_AntiReduced_game, HMS_Anti_Derived_game, HMS_Derived_game, HMS_RedGame, 
%  HMS_Reduced_game, HMS_TwoReduced_game, ImputationVertices, Kernel, LED, LED_propertyQ, LEDcons_propertyQ, LeastCore, LedcoconsQ, 
%  Ledcons_propertyQ, LorenzDom, LorenzSol, MMExcess, ModLeastCore, ModPreKernel, ModPrekernelQ, Modiclus, PModPreKernel, PModPrekernelQ, 
%  PRP_propertyQ, PreKernel, PreNucl, REAS_LED_DCGame, REAS_propertyQ, REC_propertyQ, RE_RGP, Reconfirmation_propertyQ, Reduced_game_propertyQ, 
%  SDCP_propertyQ, SDFP_propertyQ, SD_ShapleyValue, SED, SED_propertyQ, SedcoconsQ, Sedcons_propertyQ, StrConverse_DGP_Q, StrConverse_RGP_Q, 
%  StrLedcoconsQ, StrSedcoconsQ, SuperSets, TauValue, UtopiaPayoff, airport_game, anti_coreQ, balancedCollectionQ, balancedQ, belongToGameSpan, 
%  clToMatlab, clp_kernel, cls_kernel, cplex_AntiNucl, cplex_AntiNucl_llp, cplex_AntiPreNucl, cplex_AntiPreNucl_llp, cplex_LeastCore, cplex_exact_game, 
%  cplex_kernel, cplex_modiclus, cplex_nucl, cplex_nucl_llp, cplex_nucl_mod, cplex_prekernel, cplex_prenucl, cplex_prenucl_llp, cvx_kernel, 
%  cvx_prekernel, equal_treatmentQ, game_basis, game_space, getgame, glpk_AntiNucl, glpk_AntiNucl_llp, glpk_AntiPreNucl, glpk_AntiPreNucl_llp, 
%  glpk_kernel, glpk_nucl, glpk_nucl_llp, glpk_prekernel, glpk_prenucl, gurobi_AntiNucl, gurobi_AntiNucl_llp, gurobi_AntiPreNucl, gurobi_AntiPreNucl_llp, 
%  gurobi_AssignmentGame, gurobi_kernel, gurobi_nucl, gurobi_nucl_llp, gurobi_prekernel, gurobi_prenucl, gurobi_prenucl_llp, homogeneous_representationQ, 
%  hsl_prekernel, ipopt_kernel, ipopt_prekernel, lin_prekernel, mcst_game, modiclusQ, msk_AntiNucl, msk_AntiNucl_llp, msk_AntiPreNucl, msk_AntiPreNucl_llp, 
%  msk_exact_game, msk_kernel, msk_nucl, msk_nucl_llp, msk_prekernel, msk_prenucl, msk_prenucl_llp, oases_kernel, oases_prekernel, ols_prekernel, 
%  ps_gameQ, psstar_gameQ, qpBB_kernel, qpc_kernel, qpc_prekernel, qrg_prekernel, reasonable_outcome, replicate_Shapley, replicate_prk, 
%  select_starting_pt, sm_Kernel, sm_PreKernel, sm_PreNucl, sm_nucl,
%
%  
%
% PARALLEL
%  p_AllMarginalContributions, p_Anti_Converse_DGP_Q, p_Anti_Derived_game_propertyQ, p_Anti_ModPreKernel, p_Anti_ModPrekernelQ, p_Anti_PModPrekernelQ, 
%  p_Anti_PropModPreKernel, p_COV_propertyQ, p_CddTotallyBalancedQ, p_Converse_RGP_Q, p_DM_AntiReduced_game, p_DM_Anti_Derived_game, p_DM_Reduced_game,
%  p_DualCover, p_DualFloor, p_ECCoverGame, p_ECFloorGame, p_EC_DGP_Q, p_EC_RGP_Q, p_EC_propertyQ, p_HMS_AntiReduced_game, p_HMS_Anti_Derived_game, 
%  p_HMS_Derived_game, p_HMS_Reduced_game, p_Kernel, p_LedcoconsQ, p_Ledcons_propertyQ, p_ModHoller, p_ModPreKernel, p_ModPrekernelQ, p_PModPreKernel, 
%  p_PModPrekernelQ, p_PreKernel, p_PrekernelQ, p_PropModPreKernel, p_REAS_propertyQ, p_REC_propertyQ, p_Reconfirmation_propertyQ, p_Reduced_game_propertyQ, 
%  p_SD_ShapleyValue, p_SedcoconsQ, p_Sedcons_propertyQ, p_StrConverse_DGP_Q, p_StrConverse_RGP_Q, p_StrLedcoconsQ, p_StrSedcoconsQ, p_airport_profit, 
%  p_clp_kernel, p_cls_kernel, p_concave_gameQ, p_cplex_exact_game, p_cplex_kernel, p_cplex_prekernel, p_cvx_kernel, p_cvx_prekernel, p_glpk_kernel, 
%  p_glpk_prekernel, p_gurobi_AssignmentGame, p_gurobi_exact_game, p_gurobi_kernel, p_gurobi_prekernel, p_homogeneous_representationQ, p_hsl_prekernel, 
%  p_ipopt_kernel, p_ipopt_prekernel, p_lin_prekernel, p_msk_AssignmentGame, p_msk_exact_game, p_msk_kernel, p_msk_prekernel, p_oases_kernel, 
%  p_oases_prekernel, p_ols_prekernel, p_qpBB_kernel, p_qpc_kernel, p_qpc_prekernel, p_qrg_prekernel, p_semi_concaveQ, p_totallyBalancedQ
%
%
%
%
%
% About Version 1.8
% -----------------
% Release Date: 28-Jan-2020
%
% 
% SUMMARY
% A major highlight of Version 1.8 is the implementation of a set
% of functions to compute the modified and proper modified
% pre-kernel of a game. In this respect, the axiomatization of
% these solutions and of the modiclus has been implemented.
% In addition, some functions have modified. Finally, some minor 
% bugs were fixed. 
%
%
% NEW FUNCTIONS
%
% AUX: Some auxiliary files
%   * NONE
%
%  SERIAL COMPUTING
%
%   * AntiReduced_game_propertyQ.m		- Checks whether an imputation x satisfies the anti-reduced game property. 
%   * Anti_B0_balancedCollectionQ.m		- Verifies whether the set of induced coalitions is a anti B0_balanced collection.
%   * Anti_BestCoalitions.m			- Computes the set of less effective coalitions.
%   * Anti_Converse_DGP_Q.m			- Checks whether an imputation x satisfies the anti-converse derived game property.
%   * Anti_DerivedGame.m			- Computes from (v,x,S) a modified Davis-Maschler anti-derived game vS on S at x for game v.
%   * Anti_Derived_game_propertyQ.m		- Checks whether an imputation x satisfies a modified anti-derived game property.
%   * Anti_ModPreKernel.m			- Computes from (v,x) an anti-modified pre-kernel element.
%   * Anti_ModPrekernelQ.m			- Checks whether the imputation x is a modified anti-pre-kernel element of the TU-game v.
%   * Anti_PModPreKernel.m			- Computes from (v,x) an anti-proper-modified pre-kernel element.
%   * Anti_PModPrekernelQ.m			- Checks whether the imputation x is a proper modified anti-pre-kernel element of the TU-game v.
%   * Anti_PropModPreKernel.m                   - Checks whether the imputation x is a proper modified anti-pre-kernel element of the TU-game v.	    
%   * Anti_balancedCollectionQ.m		- Verifies whether the set of induced coalitions is an anti balanced collection.
%   * CddModLeastCoreVertices.m                 - Computes the vertices of the least core of game v using CDDMEX.
%   * CddTotallyBalancedQ.m			- Checks whether the core of all subgames is non-empty (cddmex).
%   * Converse_DGP_Q.m                          - Checks whether an imputation x satisfies the converse derived game property.
%   * DCP_propertyQ.m                           - Checks wheter the solution x satisfies the dual cover property.	    
%   * DFP_propertyQ.m                           - Checks wheter the solution x satisfies the dual floor property.	    
%   * DM_AntiReduced_game.m			- Computes from (v,x) all anti-reduced games on S at x of game v.
%   * DM_Anti_Derived_game.m                    - Computes from (v,x) a modified Davis-Maschler anti-reduced game vS on S at x for game v.
%   * DM_Derived_game.m                         - Computes from (v,x) a modified Davis-Maschler reduced game vS on S at x for game v.
%   * DM_TwoReduced_game.m			- Computes from (v,x) all single and two-person reduced games on S at x of game v.
%   * DRP_propertyQ.m                           - Checks wheter the solution x satisfies the dual replication property.	    
%   * DerivedGame.m				- Computes from (v,x,S) a modified Davis-Maschler derived game vS on S at x for game v.
%   * Derived_game_propertyQ.m                  - Checks whether an imputation x satisfies a modified derived game property.
%   * DualCover.m				- Computes the maximum characteristic values from the primal or dual game.
%   * DualFloor.m				- Computes the minimum characteristic values from the primal or dual game.
%   * Dual_Cover_game.m                         - Computes from (v,x) a modified Davis-Maschler reduced game vS on S at x for game v.
%   * Dual_Cover_propertyQ.m                    - Checks whether an imputation x satisfies a modified reduced game property	    
%   * DuttaRay.m				- Computes the Dutta-Ray solution for convex games.
%   * ECCoverGame.m				- Computes from (v,x) an excess comparability cover of game v.
%   * ECFloorGame.m				- Computes from (v,x) an excess comparability floor of game v.
%   * EC_DGP_Q.m				- Checks whether the solution x satisfies excess comparability for each derived game.
%   * EC_RGP_Q.m				- Checks whether the solution x satisfies excess comparability for each reduced game.
%   * EC_propertyQ.m                            - Checks whether the solution x satisfies excess comparability.	    
%   * HMS_AntiReduced_game.m                    - Computes from (v,x) all Hart/Mas-Colell anti-reduced games on S at x of game v.	    
%   * HMS_Anti_Derived_game.m                   - Computes from (v,x,S) a modified Hart-Mas-Colell anti-reduced game vS on S at x for game v.	    
%   * HMS_Derived_game.m			- Computes from (v,x,S) a modified Hart-Mas-Colell reduced game vS on S at x for game v.
%   * HMS_DervGame.m                            - Computes from (v,x,S) a modified Hart-Mas-Colell derived game vS on S at x for game v.	    
%   * HMS_RedGame.m				- Computes from (v,x,S) a Hart-MasColell reduced game vS on S at x for game v.
%   * HMS_TwoReduced_game.m			- Computes from (v,x) all Hart/Mas-Colell singleton and two-person reduced games on S at x of game v.
%   * LED.m					- Computes the large excess difference w.r.t. the payoff x.
%   * LED_propertyQ.m                           - Checks wheter the solution x satisfies large excess difference property.	    
%   * LEDcons_propertyQ.m			- Checks whether an imputation x satisfies large excess difference converse consistency
%   * LedcoconsQ.m				- Checks whether an imputation x satisfies large excess difference converse consistency.
%   * Ledcons_propertyQ.m			- Checks whether an imputation x satisfies the ledcons property 
%   * LorenzSet.m                               - Determines the Lorenz set of game v.
%   * LorenzSol.m				- Determines the Lorenz solution of game v.
%   * MMExcess.m				- Computes the minimal and maximal excess vector of game v and its dual.
%   * ModHoller.m				- Computes a modified Holler index from the set of winning coalitions.
%   * ModPreKernel.m                            - Computes from (v,x) a modified pre-kernel element.	    
%   * ModPrekernelQ.m                           - Checks whether the imputation x is a modified pre-kernel element of the TU-game v.	    
%   * PModPreKernel.m                           - Computes from (v,x) a proper modified pre-kernel element.	    
%   * PModPrekernelQ.m                          - Checks whether the imputation x is a proper modified pre-kernel element of the TU-game v.
%   * PRP_propertyQ.m                           - Checks wheter the solution x satisfies the primal replication property.	    
%   * PrkEqsModPrkQ.m                           - Checks whether a pre-kernel element is also an element of the modified as well as proper modified pre-kernel	    
%   * PropModPreKernel.m			- Computes from (v,x) a proper modified pre-kernel element from the dual cover game.
%   * REAS_LED_DCGame.m                         - Verifies that x is a reasonable vector of game v, then the shifted ducal cover game statisfies LED w.r.t. the replicated vector (x,x).
%   * REAS_propertyQ.m                          - Checks if the vector x satisfies the reasonableness on both sides
%   * REC_propertyQ.m                           - Checks wheter the solution x satisfies reverse excess comparability.	    
%   * RE_RGP.m                                  - Checks whether an imputation x is reasonable from both sides for all reduced games.
%   * SDCP_propertyQ.m                          - Checks wheter the solution x satisfies a strong dual cover property.
%   * SDFP_propertyQ.m                          - Checks wheter the solution x satisfies a strong dual floor property.
%   * SD_ShapleyValue.m                         - Computes the surplus division Shapley value. 
%   * SED.m					- Computes the small excess difference w.r.t. the payoff x.
%   * SED_propertyQ.m                           - Checks wheter the solution x satisfies small excess difference property.	    
%   * SedcoconsQ.m				- Checks whether an imputation x satisfies small excess difference converse consistency.
%   * Sedcons_propertyQ.m			- Checks whether an imputation x satisfies the sedcons property.
%   * StrConverse_DGP_Q.m			- Checks whether an imputation x satisfies the strong converse derived game property.
%   * StrLedcoconsQ.m                           - Checks whether an imputation x satisfies strong large excess difference converse consistency.	    
%   * StrSedcoconsQ.m                           - Checks whether an imputation x satisfies strong small excess difference converse consistency (SEDCOCONS).	    
%   * anti_coreQ.m				- Checks the existence of the anti-core of game v.
%   * anti_game_space.m                         - Computes the game space which replicates x as an anti-pre-kernel element.
%   * belongToModLeastCoreQ.m                   - Checks whether the imputation x belongs to the core of game v.	    
%   * cplex_cs_modiclus.m			- Computes the modiclus of game v w.r.t. a coaliton structure cs using cplexmex.
%   * exact_gameQ.m				- Checks whether game v is an exact game using Matlab's Optimization toolbox.
%   * formatPowerSet.m                          - Formats the Matlab cell output that contains the representation of coalitions into matrix form.
%   * replicate_aprk.m                          - Replicates the anti-pre-kernel solution as an anti-pre-kernel of the game space v_sp.
%   * shiftGame.m				- Computes from the game v the t-shift game of v.
%   * sm_PreKernel.m                            - Computes an element of the simplified pre-kernel of game v.	    
%   * totallyBalancedQ.m			- Checks whether the core of all subgames is non-empty.
%   * value_matrix.m                            - Computes from an assignment matrix the corresponding value matrix for a permutation game.	    
%
%
%  PARALLEL COMPUTING
%
%   * p_AntiReduced_game_propertyQ.m	        - Checks whether an imputation x satisfies the anti-reduced game property.
%   * p_Anti_Converse_DGP_Q.m		        - Checks whether an imputation x satisfies the anti-converse derived game property.
%   * p_Anti_Derived_game_propertyQ.m	        - Checks whether an imputation x satisfies a modified anti-derived game property.
%   * p_Anti_ModPreKernel.m			- Computes from (v,x) a modified pre-kernel element.
%   * p_Anti_ModPrekernelQ.m		        - Checks whether the imputation x is a modified anti-pre-kernel element of the TU-game v.
%   * p_Anti_PModPreKernel.m		        - Computes from (v,x) a proper modified anti-pre-kernel element.
%   * p_Anti_PModPrekernelQ.m		        - Checks whether the imputation x is a proper modified anti-pre-kernel element of the TU-game v.
%   * p_Anti_PropModPreKernel.m		        - Computes from (v,x) a proper modified anti-pre-kernel element.
%   * p_COV_propertyQ.m			        - Verifies if the payoff x satisfies the covariance with strategic equivalence property w.r.t. (v,m,t).
%   * p_CddTotallyBalancedQ.m		        - Checks whether the core of all subgames is non-empty (cddmex).
%   * p_Converse_DGP_Q.m			- Checks whether an imputation x satisfies the converse derived game property.
%   * p_DM_AntiReduced_game.m		        - Computes from (v,x) all anti-reduced games on S at x of game v.
%   * p_DM_Anti_Derived_game.m		        - Computes from (v,x) a modified Davis-Maschler anti-reduced game vS on S at x for game v.
%   * p_DM_Derived_game.m			- Computes from (v,x) a modified Davis-Maschler reduced game vS on S at x for game v.
%   * p_Derived_game_propertyQ.m		- Checks whether an imputation x satisfies a modifed derived game property.
%   * p_DualCover.m				- The maximum characteristic values from the primal or dual game.
%   * p_DualFloor.m				- The minimum characteristic values from the primal or dual game.
%   * p_ECCoverGame.m			        - Computes from (v,x) an excess comparability cover of game v.
%   * p_ECFloorGame.m			        - Computes from (v,x) an excess comparability floor of game v.
%   * p_EC_DGP_Q.m				- Checks whether the solution x satisfies excess comparability for each derived game.
%   * p_EC_RGP_Q.m				- Checks whether the solution x satisfies excess comparability for each reduced game.
%   * p_EC_propertyQ.m			        - Checks wheter the solution x satisfies excess comparability.
%   * p_HMS_AntiReduced_game.m		        - Computes from (v,x) all Hart/Mas-Colell anti-reduced games on S at x of game v.
%   * p_HMS_Anti_Derived_game.m		        - Computes from (v,x,S) a modified Hart-Mas-Colell anti-reduced game vS on S at x for game v.
%   * p_HMS_Derived_game.m			- Computes from (v,x,S) a modified Hart-Mas-Colell reduced game vS on S at x for game v.
%   * p_LedcoconsQ.m			        - Checks whether an imputation x satisfies large excess difference converse consistency. 
%   * p_Ledcons_propertyQ.m			- Checks whether an imputation x satisfies the ledcons property.
%   * p_ModHoller.m				- Computes the modified Holler index from the set of winning coalitions.
%   * p_ModPreKernel.m			        - Computes from (v,x) a modified pre-kernel element.
%   * p_ModPrekernelQ.m			        - Checks whether the imputation x is a modified pre-kernel element.
%   * p_PModPreKernel.m			        - Computes from (v,x) a proper modified pre-kernel element.
%   * p_PModPrekernelQ.m			- Checks whether the imputation x is a proper modfied pre-kernel element.
%   * p_PropModPreKernel.m			- Computes from (v,x) a proper modified pre-kernel element.
%   * p_REAS_propertyQ.m			- Checks if the vector x satisfies the reasonableness on both sides.
%   * p_REC_propertyQ.m			        - Checks wheter the solution x satisfies reverse excess comparability.
%   * p_SD_ShapleyValue.m			- Computes the surplus division Shapley value.
%   * p_SedcoconsQ.m			        - Checks whether an imputation x satisfies small excess difference converse consistency.
%   * p_Sedcons_propertyQ.m			- Checks whether an imputation x satisfies the sedcons property.
%   * p_StrConverse_DGP_Q.m			- Checks whether an imputation x satisfies the strong converse derived game property.
%   * p_StrConverse_RGP_Q.m			- Checks whether an imputation x satisfies the strong CRGP.
%   * p_StrLedcoconsQ.m			        - Checks whether an imputation x satisfies strong large excess difference converse consistency.
%   * p_StrSedcoconsQ.m			        - Checks whether an imputation x satisfies strong small excess difference converse consistency.
%   * p_totallyBalancedQ.m			- checks whether the core of all subgames is non-empty.
%
%
%
%
% BUG FIXES
%  Incorrect computation of LorenzDom due to unsorted vectors has been fixed. 
%
%
% MODIFIED FUNCTIONS
% SERIAL
%   B0_balancedCollectionQ,B0_balancedQ,COV_propertyQ,CddCoreQ,CddCoreVertices,CddModiclus,HMS_Reduced_game,
%   PreKernel,PreNucl,RedGame,Reduced_game_propertyQ,Weak_balancedCollectionQ,average_convexQ,balancedCollectionQ,
%   balancedQ,cls_kernel,convex_gameQ,cplex_AntiNucl,cplex_AntiNucl_llp,cplex_AntiPreNucl,cplex_AntiPreNucl_llp,
%   cplex_LeastCore,cplex_exact_game, cplex_modiclus,cplex_nucl,cplex_nucl_llp,cplex_prenucl,cplex_prenucl_llp,
%   cplex_weightedPreNucl,cplex_weightedPreNucl_llp,glpk_modiclus,gurobi_prenucl,holler,ipopt_kernel,modiclusQ,
%   msk_AntiPreNucl,oases_prekernel,vclToMatlab
%
%
% PARALLEL
%
%   p_Anti_PreKernel,p_HMS_Reduced_game,p_PreKernel,p_Reduced_game_propertyQ,p_average_convexQ,
%   p_cplex_exact_game,p_cplex_flow_game,p_holler,p_ipopt_kernel,p_msk_prekernel,p_oases_prekernel,p_winning_coalitions  
%
%
% About Version 1.7.5
% -----------------
% Release Date: 22-Jan-2019
%
% 
% SUMMARY
% A major highlight of Version 1.7.5 is the implementation of a set
% of functions to compute the modiclus as well as the simplified
% (pre-)nucleolus and simplified (pre-)kernel of a game.  Furthermore, 
% a set of functions to compute a minimum cost spanning has been added. 
% Some minor bugs were fixed. 
%
%
% NEW FUNCTIONS
%
%
%
%  SERIAL COMPUTING
%
%   * Anti_B0_balancedCollectionQ               - Checks the reversal of weak Kohlberg's criterion.
%   * Anti_Modiclus                             - Computes the anti modiclus of a game.
%   * Anti_Nucl                                 - Computes the anti nucleolus of a game.
%   * Anti_Nucl_llp                             - Computes the anti nucleolus of a game (fast).
%   * Anti_PreNucl_llp                          - Computes the anti pre-nucleolus of game v.
%   * Anti_balancedCollectionQ                  - Checks the reversal of Kohlberg's criterion.
%   * Anti_modiclusQ                            - Verifies whether the set of induced coalitions is a bi-balanced collection.
%   * B0_balancedCollectionQ                    - Checking weak Kohlberg's criterion.
%   * B0_balancedQ                              - Verifies whether the collection of coalitions is weakly balanced.
%   * CddAntiLeastCore                          - Computes the least core of game v using (cddmex).
%   * CddAntiLeastCoreVertices                  - Computes the vertices of the anti least core of game v (cddmex).
%   * CddAntiNucl                               - Computes the anti nucleolus of game v (cddmex).
%   * CddAntiNucl_llp                           - Computes the anti nucleolus of game v (cddmex).
%   * CddAntiPrenucl_llp                        - Computes the anti pre-nucleolus of game v (cddmex).
%   * CddExactGame                              - Computes the exact game from v (cddmex).
%   * CddModiclus                               - Computes the modiclus of game v using cddmex.
%   * Modiclus                                  - Computes the modiclus of a game.
%   * One_Normalization                         - Computes from the game v the corresponding one-normalized game.
%   * Weak_balancedCollectionQ                  - Checking weak Kohlberg's criterion.
%   * average_excess                            - Computes the average excess of game v.
%   * concave_gameQ                             - Checks the concavity of a Tu-game.
%   * cplex_AntiNucl                            - Computes the anti nucleolus of game v using the CPLEX solver.
%   * cplex_AntiNucl_llp                        - Computes the anti nucleolus of game v using the CPLEX solver.
%   * cplex_AntiPreNucl_llp                     - Computes the anti prenucleolus using the CPLEX solver.
%   * cplex_exact_game                          - Computes the exact game from v using the CPLEX solver.
%   * cplex_modiclus                            - Computes the modiclus of game v using cplexmex.
%   * exact_game                                - Computes the exact game from v using Matlab's Optimization toolbox.
%   * glpk_AntiNucl                             - Computes the anti nucleolus of game v using the GLPK solver.
%   * glpk_AntiNucl_llp                         - Computes the anti nucleolus of game v using the GLPK solver.
%   * glpk_AntiPreNucl_llp                      - Computes the anti pre-nucleolus using the GLPK solver.
%   * glpk_exact_game                           - Computes the exact game from v using the GLPK solver.
%   * glpk_modiclus                             - Computes the modiclus of game v using glpkmex.
%   * gurobi_AntiNucl                           - Computes the anti nucleolus of game v using the GUROBI solver.
%   * gurobi_AntiNucl_llp                       - Computes the anti nucleolus of game v using the GUROBI solver.
%   * gurobi_AntiPreNucl_llp                    - Computes the anti prenucleolus using the GUROBI solver.
%   * gurobi_exact_game                         - Computes the exact game from v using the GUROBI solver.
%   * gurobi_modiclus                           - Computes the modiclus of game v using the GUROBI.
%   * ireffq                                    - Checks if a payoff satisfies IR as well as the Eff.
%   * k_anticover                               - Determines from the Tu-game the corresponding anti k-game.
%   * k_concaveQ                                - Checks k-concavity of the Tu-game.
%   * mcst_game                                 - Computes from a cost matrix the corresponding mcst game.
%   * minNoBlockPayoff                          - Computes the minimum no blocking payoff from game v.
%   * modiclusQ                                 - Verifies whether the set of induced coalitions is a bi-balanced collection.
%   * msk_AntiNucl                              - Computes the anti nucleolus of game v using the MOSEK solver.
%   * msk_AntiNucl_llp                          - Computes the anti nucleolus of game v using the MoSEK solver.
%   * msk_AntiPreNucl_llp                       - Computes the anti prenucleolus using the MOSEK solver.
%   * msk_exact_game                            - Computes the exact game from v using the MOSEK solver.
%   * msk_modiclus                              - Computes the modiclus of game v using the MOSEK solver.
%   * sm_Kernel.m                               - Computes an element of the simplified Kernel of a game.
%   * sm_PreKernel.m                            - Computes an element of the simplified Pre-kernel of a game.
%   * sm_PreNucl                                - Computes the simplified pre-nucleolus of a game.
%   * sm_nucl                                   - Computes the simplified nucleolus of a game.
%
%
%
%  PARALLEL COMPUTING
%
%   * p_AntiB0_balancedCollectionQ              - Checks the reversal of weighted Kohlberg's criterion.
%   * p_B0_balancedCollectionQ                  - Checking weak Kohlberg's criterion.
%   * p_CddExactGame                            - Computes the exact game from v (cddmex).
%   * p_balancedSetQ                            - Verifies whether the set of induced coalitions is a balanced collection.
%   * p_cplex_exact_game                        - Computes the exact game from v using the CPLEX solver.
%   * p_exact_game                              - Computes the exact game from v using Matlab's Optimization toolbox.
%   * p_glpk_exact_game                         - Computes the exact game from v using the GLPK solver.
%   * p_gurobi_exact_game                       - Computes the exact game from v using the GUROBI solver.
%   * p_mcst_game                               - Computes from a cost matrix the corresponding mcst game.
%   * p_msk_exact_game                          - Computes the exact game from v using the MOSEK solver.
%
%
%
%
% Class Objects
% -------------
%  SERIAL COMPUTING
%   * TuAPrn                                    - subclass object of TuSol (anti pre-nucleolus from various solvers).
%   * TuASol                                    - subclass object of TuGame (game solutions).
%   * TuKrn                                     - subclass object of TuSol (kernel solutions from various solvers).
%   * TuNuc                                     - subclass object of TuSol (nucleolus from various solvers).
%
%  PARALLEL COMPUTING
%   * p_TuKrn                                   - subclass object of p_TuSol (kernel solutions from various solvers).
%
%
%
% BUG FIXES
%  Incorrect computations due to a wrong payoff space specifications related to the Anti_Kernel
%  function have been fixed. 
%
% MODIFIED FUNCTIONS
% SERIAL
%   balancedCollectionQ,clToMatlab,cvx_prekernel,ipopt_kernel,market_game,msk_kernel,msk_prekernel,
%   msk_weightedNucl_llp,msk_weightedNucl,oases_kernel,oases_prekernel,nucl,nucl_llp,cplex_nucl_llp,
%   cplex_nucl,gurobi_nucl_llp,gurobi_nucl,glpk_nucl,glpk_nucl_llp,msk_nucl,msk_nucl_llp,CddNucl,
%   CddNucl_llp,msk_prenucl,msk_prenucl_llp,glpk_prenucl,glpk_prenucl_llp,gurobi_prenucl_llp,gurobi_prenucl,
%   cplex_prenucl,cplex_prenucl_llp,PreNucl_llp,CddPrenucl_llp,CddPrenucl,PreNucl,weightedBalancedCollectionQ,
%   CddAntiNucl,glpk_AntiNucl,gurobi_AntiNucl,msk_AntiNucl,cplex_AntiNucl,Anti_Nucl,Anti_Nucl_llp,
%   msk_AntiNucl_llp,glpk_AntiNucl_llp, CddAntiNucl_llp,gurobi_AntiNucl_llp,cplex_AntiNucl_llp,
%   cplex_AntiPreNucl_llp,cplex_AntiPreNucl,zero_monotonicQ.
%
% PARALLEL
%   p_ipopt_kernel,p_oases_kernel,p_oases_prekernel,p_ISRG_propertyQ,p_msk_kernel,p_msk_prekernel.
%
%
%
%
% About Version 1.7 (outdated version numbering: 1.0)
% -----------------
% Release Date: 24-Apr-2018
%
% 
% SUMMARY
% The functions to compute the (pre-)nucleolus/kernel and the
% associated function set of third party solvers like CPLEX, GUROBI
% and MOSEK have been modified. Some minor bug fixes. 
%
% About Version 0.8
% -----------------
% Release Date: 21-Dec-2015
%
% 
% SUMMARY
% The set of graphical functions has been extended, for instance,
% the core-cover, and some kernel catchers can now be plotted. The
% anti-core functions have been revised, and are more robust. A
% simplex projection method from 4-d into 3-d has been added.
% Some minor bug fixes. 
%
% NEW FUNCTIONS
%
%
% AUX: Some auxiliary files
%
%   * FrameToImage                        - Converts a frame to an image.
%   * PlayCoreMovie                       - Plays a movie from a collection of frames.
%   * SaveFrames                          - Saves converted image files to the hard disk.
%   * ToSimplex3d                         - Projects data from 4d to 3d.
%   * testcase_graphics                   - Checking basic graphic installation.
%
%  SERIAL COMPUTING
%   * Anti_balancedCollectionQ            - Checks the reversal of Kohlberg's criterion
%   * AntiImputationVertices              - Computes all vertices of the anti imputation set.
%   * balancedQ                           - Verifies whether the collection of coalitions is balanced.
%   * CddAntiCoreSimplexPlot              - Plots the anti-core using simplex projection.
%   * CddAntiCoreSimplexVertices          - Computes all anti-core vertices using simplex projection.
%   * CddAntiImputationSimplexVertices    - Computes all vertices of the anti-imputation set using simplex projection.
%   * CddAntiImputationVertices           - Computes all vertices of the anti-imputation set.
%   * CddAntiPrenucl                      - Computes the anti pre-nucleolus of game v (cddmex).
%   * CddCoreCoverSimplexPlot             - Plots the core cover (simplex projection).
%   * CddCoreCoverSimplexVertices         - Computes all vertices of the core cover (simplex).
%   * CddCoreSimplexMovie                 - Creates a movie w.r.t. the strong epsilon-cores (simplex projection).
%   * CddCoreSimplexPlot                  - Plots the core (simplex projection).  
%   * CddCoreSimplexVertices              - Computes the vertices of the core (simplex).
%   * CddImputationSimplexVertices        - Computes the vertices of the imputation set (simplex).  
%   * CddKernelCatchersSimplex            - Draws some kernel catchers (simplex).  
%   * CddLowerSetSimplexVertices          - Computes the vertices of the lower set (simplex).   
%   * CddReasonableSetSimplexVertices     - Computes the vertices of the reasonable set (simplex).   
%   * CddStrongCoreSimplexPlot            - Plots the strong epsilon core (simplex projection).  
%   * CddUpperSetSimplexVertices          - Computes the vertices of the upper set (simplex).   
%   * CddWeberSetSimplex                  - Computes the vertices of the Weber Set (simplex).  
%   * CddWeberSetSimplexPlot              - Plots the Weber set (simplex).  
%   * COV_propertyQ                       - Verifies if the payoff x satisfies COV property.  
%   * cplex_AntiPreNucl                   - Computes the anti prenucleolus using the CPLEX solver.  
%   * cplex_LeastCore                     - Computes the least core using cplexmex.  
%   * DecomposeGame                       - Computes the unique decomposition of a TU-game.
%   * DiscShapleyValue                    - Computes the discounted Shapley value.  
%   * glpk_AntiPreNucl                    - Computes the anti pre-nucleolus using the GLPK solver.  
%   * gurobi_AntiPreNucl                  - Computes the anti prenucleolus using the GUROBI solver.  
%   * HMS_ImputSavingReducedGame          - Computes from (v,x) all Hart/Mas-Colell ISR games.
%   * ImpSetEqsLwsQ                       - Checks if the imputation set coincides with the lower set. 
%   * ImputSavingReducedGame              - Computes from (v,x) all imputation saving reduced games.  
%   * interest_game                       - Computes from an interest problem the corresponding game.  
%   * ISRG_propertyQ                      - Checks whether an imputation x satisfies the ISR game property.  
%   * KrEqsPrkQ                           - Checks if the kernel is equal to the pre-kernel.  
%   * msk_AntiPreNucl                     - Computes the anti prenucleolus using the MOSEK solver.  
%   * PropNucl                            - Computes the proportional nucleolus.  
%   * PropPreNucl                         - Computes the proportional pre-nucleolus.  
%   * ReasSetEqsUpsQ                      - Checks if the reasonable set coincides with the upper set.  
%
%  PARALLEL COMPUTING
%   * p_COV_propertyQ                     - Verifies if the payoff x satisfies COV property.
%   * p_HMS_ImputSavingReducedGame        - Computes from (v,x) all Hart/Mas-Colell ISR games.
%   * p_ImputSavingReducedGame            - Computes from (v,x) all imputation saving reduced games.
%   * p_ISRG_propertyQ                    - Checks whether an imputation x satisfies the ISR game property.
%
%
% MODIFIED FUNCTIONS
% SERIAL
%   The set of graphical functions has been revised. 
%   AntiCoreVertices, AntiImputationVertices, coreQ,
%   CddAntiCoreVertices, CddAntiImputationVertices,
%   CddImputationVertices, CddCoreQ,
%   LeastCore, UtopiaPayoff. 
%
%
%
% These functions have now a modified behavior:
%
% About Version 0.7
% -----------------
% Release Date: 06-Apr-2015
%
% 
% SUMMARY
% All functions that are using the linprog command from
% Matlab's Optimization Toolbox apply now a dual-simplex method.
% As a consequence, these functions are not any more backward compatible.
% Some minor bug fixes. 
%
% NEW FUNCTIONS
%
%
%
% MODIFIED FUNCTIONS
% SERIAL
%   balancedCollectionQ, LeastCore, nucl, PreNucl, PrenuclQ.
%
%
%
% These functions have now a modified behavior:
%   * balancedCollectionQ           - Dual-Simplex; only three input arguments are allowed.
%   * LeastCore                     - Dual-Simplex.
%   * nucl                          - Dual-Simplex.
%   * PreNucl                       - Dual-Simplex.
%   * PrenuclQ                      - Dual-Simplex; only three input arguments are allowed.
%
%
%
% About Version 0.6
% -----------------
% Release Date: 06-Mar-2015
%
% 
% SUMMARY
% Extended the graphical capabilities of the toolbox by adding
% some new function.
%
%
%
%  SERIAL COMPUTING
%   * proper_amount                       - Computes the max amount players contribute to a proper coalition.
%   * scrb_solution                       - Computes separable costs-remaining benefits allocation. 
%   * smallest_amount                     - Computes the smallest amount vector.
%
%
%
%
% BUG FIXES
%   * separable_cost_allocation                          - Fix an incorrect computation for method 'SCRB'.
%
% MODIFIED FUNCTIONS
% SERIAL
%   * balancedCollectionQ
%
%
%
% About Version 0.5
% -----------------
% Release Date: 03-Nov-2014
%
% SUMMARY
% Complete code revision of the function ShapleyValue. It is
% now in average ten times faster, but needs more memory. Added new
% graphic functions and some code revision for the old ones.
% They require now the Multi-Parametric Toolbox 3, and their behavior have
% changed. Finally, a set of new functions has been added and some
% bugs were fixed.
%
%
% NEW FUNCTIONS
%
%
%  SERIAL COMPUTING
%   * balancedCollectionQ                 - Checking Kohlberg's criterion.
%   * CddCoreCoverPlot                    - Plots the core cover of a TU game.
%   * CddCoreCoverVertices                - Computes all vertices of the core cover of a TU game.
%   * CddCoreMovie                        - Creates a movie w.r.t. the strong epsilon-cores.
%   * CddStrongCorePlot                   - Plots a strong epsilon core. 
%   * CddWeberSet                         - Computes the vertices of the Weber Set.
%   * CddWeberSetPlot                     - Plots the Weber set.
%   * holler                              - Computes the Holler index.
%   * MLextension                         - Computes the multilinear extension.
%   * PermutationGame                     - Computes from an assignment matrix the permutation game.
%   * Potential                           - Determines the potential of a TU game (recursive).
%   * PrenuclQ                            - Checks if a payoff vector is the pre-nucleolus using Kohlberg's criterion.
%   * ShapleyValueML                      - Computes the Shapley value using multilinear extension. 
%
%
%  PARALLEL COMPUTING
%   * p_holler                            - Computes the Holler index.
%   * p_DecomposeGame                     - Computes the unique decomposition of a TU-game.
%   * p_PermutationGame                   - Computes from an assignment matrix the permutation game.
%   * p_potential                         - Determines the potential of a TU game (basis).
%   * p_ShapleyValue                      - Computes the Shapley value (potential).
%
%
%  MATHEMATICA SYMBOLIC TOOLBOX FUNCTIONS
%
%   * tug_DetRandCoord                    - Returns random unanimity coordinates.
%   * tug_MLExtension                     - Computes the multilinear extension.
%   * tug_ShapleyValueML                  - Determines the Shapley value using multilinear extension.
%
%
% BUG FIXES
%   * TuProp                         - Fix an incorrect line break.
%   * TuRep                          - Fix a typo for a subcase.
%   * p_TuRep                        - Fix a typo for a subcase.     
%
% MODIFIED FUNCTIONS
% 
% SERIAL
%  * coeff_linearbasis, cplex_nucl, cplex_prenucl (compatible with
%    CPLEX Version 12.5.1), CddCorePlot, CorePlot,
%    CddCorePlot, CddCoreVertices, CddImputationVertices,
%    CddPreKernel, game_space, game_space, LeastCore, linear_basis,msk_prekernel,
%    msk_prenucl, PreNucl, nucl, PrekernelQ, ShapleyValue 
%
% PARALLEL
%   * p_coeff_linearbasis, p_linear_basis, p_game_space,
%     p_game_space_red, p_msk_prekernel, p_PrekernelQ, p_replicate_prk
%
% These functions have now a modified behavior:
%  * CddCorePlot                - Five input arguments are allowed. View point can now be adjusted.
%  * CddCoreVertices            - Four output and three input arguments. Requires now Multi-Parametric Toolbox 3.
%  * CddImputationVertices      - Four output and two input arguments. Requires now Multi-Parametric Toolbox 3.
%  * coeff_linearbasis          - Two input arguments are allowed. 
%                                 The second argument invokes the
%                                 computation of a full or sparse basis matrix.
%  * linear_basis               - Same as coeff_linearbasis.
%  * ShapleyValue               - Does not compute anymore the potential. Use the 
%                                 function potential() instead.
%  * p_coeff_linearbasis        - Same as coeff_linearbasis.
%  * p_linear_basis             - Same as coeff_linearbasis.
%  * p_replicate_prk            - Five input arguments are allowed. The computation 
%                                 can be invoked using sparse matrices only.
%  * p_game_space               - Same as p_replicate_prk. 
%  * p_game_space_red           - Same as p_replicate_prk.
%
%
% About Version 0.4
% -----------------
% Release Date: 02-Oct-2013
%
% SUMMARY
% Some code optimization for the Shapley value functions. 
% Performance increase is about 15 percent. A set of new
% functions were added to compute fairness or related values like
% Aumann-Dreze, Owen, Myerson, position, solidarity and coalition
% solidarity value. Moreover, some functions to replicate the 
% Shapley value for related games are incorporated. In addition,
% two new subclass objects have been added. Some functions and
% methods have been modified as well as some bugs have been fixed. 
% The documentation has been revised and updated.
%
%
% NEW FUNCTIONS
%
%  SERIAL COMPUTING
%   * ADvalue                             - Computes the Aumann-Dreze value.
%   * apex_game                           - Creates an apex game.
%   * apu_SolidarityValue                 - Determines the solidarity value w.r.t. a priori unions.
%   * basis_coordinates                   - Determines the basis coordinates of a Tu game.
%   * basis_game                          - Determines bases games.
%   * CoalitionSolidarity                 - Determines the coalition solidarity value.
%   * coeff_linearbasis                   - Determines the coefficients (dividends) of a linear basis from a TU game.
%   * genUnionStable                      - Creates a union stable system.
%   * hypergraphQ                         - Checks whether the system is a hypergraph communication situation.
%   * InteractionSets                     - Determines a system of interaction sets.
%   * linear_basis                        - Determines the linear basis of the n-person TU game space.
%   * LS_Nucl                             - Computes the least square nucleolus of a game.
%   * LS_PreNucl                          - Computes the least square pre-nucleolus of a game.
%   * market_game                         - Determines from two disjoint sets a market game.
%   * mediation_game                      - Determines a mediation game.
%   * monotonic_cover                     - Computes the monotonic cover of a TU game.
%   * MyersonValue                        - Computes the Myerson value of a Tu game.
%   * OwenValue                           - Computes the Owen value.
%   * PartitionSA                         - Computes a partition of S w.r.t. a hypergraph communication situation.
%   * PartitionSL                         - Computes a partition of S w.r.t. a communication situation.
%   * PositionValue                       - Computes the position value.
%   * PowerSet                            - Computes all subsets from a set representation.
%   * production_game                     - Creates an affine production game.
%   * production_game_sq                  - Creates a quadratic production game.
%   * pure_overhead                       - Creates the matrix of pure overhead games.
%   * quotas                              - Determines the quotas of a game.
%   * replicate_Shapley                   - Replicates the Shapley value for a game space. 
%   * ShapleyValueLB                      - Computes the Shapley value from the linear basis.
%   * SolidarityShapleyValue              - Determines the solidarity Shapley value. 
%   * SolidarityValue                     - Determines the solidarity value. 
%   * streps_value                        - Determines the strong epsilon-game.
%   * union_stableQ                       - Checks whether a system is union stable.
%   * weightedOwen                        - Computes the weighted Owen value.
%   * weightedSolidarity                  - Computes the weighted solidarity value.
%
% Class Objects
% -------------
%   * TuVal                               - subclass object of TuGame (fairness and related values).
%   * TuShRep                             - subclass object of TuSol (Shapley value replication).
%
%
%  PARALLEL COMPUTING
%   * p_ADvalue                           - Computes the Aumann-Dreze value.
%   * p_apu_SolidarityValue               - Determines the solidarity value w.r.t. a priori unions.
%   * p_coeff_linearbasis                 - Determines the coefficients (dividends) of a linear basis from a TU game.
%   * p_CoalitionSolidarity               - Determines the coalition solidarity value.
%   * p_linear_basis                      - Determines the linear basis of the n-person TU game.
%   * p_LS_Nucl                           - Computes the least square nucleolus of a game.
%   * p_LS_PreNucl                        - Computes the least square pre-nucleolus of a game.
%   * p_genUnionStable                    - Creates a union stable system.
%   * p_Kernel                            - Computes a kernel point using the optimization toolbox.
%   * p_MyersonValue                      - Computes the Myerson value of a Tu game.
%   * p_OwenValue                         - Computes the Owen value.
%   * p_PositionValue                     - Computes the position value.
%   * p_pure_overhead                     - Creates the matrix of pure overhead games.
%   * p_replicate_Shapley                 - Replicates the Shapley value for a game space. 
%   * p_ShapleyValueLB                    - Computes the Shapley value from the linear basis.
%   * p_SolidarityShapleyValue            - Determines the solidarity Shapley value. 
%   * p_SolidarityValue                   - Determines the solidarity value. 
%   * p_union_stableQ                     - Checks whether a system is union stable.
%   * p_weightedOwen                      - Computes the weighted Owen value.
%   * p_weightedShapley                   - Computes the weighted Shapley value.
%   * p_weightedSolidarity                - Computes the weighted solidarity value.
% 
% Class Objects
% -------------
%   * p_TuVal                             - subclass object of TuGame (fairness and related values).
%   * p_TuShRep                           - subclass object of p_TuSol (Shapley value replication).
%
% BUG FIXES
%   * weightedShapley                     - Fix incorrect computation (completely revised). 
%
% MODIFIED FUNCTIONS
% 
% SERIAL
% cls_kernel,nucl,PreNucl,ShapleyValue,unanimity_games
%
% PARALLEL
% p_unanimity_games
%
%
% About Version 0.3
% -----------------
% Release Date: 04-Jun-2013
%
% SUMMARY
% Code revision of the pre-kernel functions. The evaluation is in average
% 30 times faster than under the previous version. A set of new
% functions were added to compute the pre-nucleolus, nucleolus and
% the kernel. To increase the set of available solvers various interfaces to
% third party solvers have been added. Moreover, to assist the user
% in relative complex and comprehensive game investigations the
% class object TuGame with several subclasses was designed.
% Finally, the code of many functions have been modified and some bugs 
% have been fixed, see below for the details.
%
%
% NEW FUNCTIONS
%
%  SERIAL COMPUTING
%   * Anti_Kernel                         - Computes an anti-kernel point.
%   * Anti_kernelQ                        - Checks if an imputation is an anti kernel point.
%   * Anti_PrekernelQ                     - Checks if an imputation is an anti prekernel point.
%   * additive_game                       - Creates an additive game.
%   * belongToAntiCoreQ                   - Checks if a payoff vector belongs to the anti-core.
%   * CddAntiCorePlot                     - Plots the anti-core of a game using cddmex. 
%   * CddAntiCoreQ                        - Check if the anti-core exists (cddmex).
%   * CddAntiCoreVertices                 - Computes the vertices of the anti-core (cddmex).
%   * CddLeastCore                        - Computes the leastcore (cddmex).
%   * clp_kernel                          - Computes a kernel point using the CLP solver.
%   * cls_kernel                          - Computes a kernel point using the CLS solver.
%   * cplex_kernel                        - Computes a kernel point using the CPLEX solver.
%   * cplex_nucl                          - Computes the nucleolus using the CPLEX solver.
%   * cplex_prekernel                     - Computes a prekernel point using the CPLEX solver.
%   * cplex_prenucl                       - Computes the prenucleolus using the CPLEX solver.
%   * cvx_kernel                          - Computes a kernel point using the CVX solver. 
%   * cvx_prekernel                       - Computes a prekernel point using the CVX solver. 
%   * glpk_kernel                         - Computes a kernel point using the GLPK solver. 
%   * glpk_nucl                           - Computes the nucleolus using the GLPK solver.
%   * glpk_prekernel                      - Computes a prekernel point using the GLPK solver. 
%   * glpk_prenucl                        - Computes the prenucleolus using the GLPK solver.
%   * gurobi_kernel                       - Computes a kernel point using the GUROBI solver. 
%   * gurobi_nucl                         - Computes the nucleolus using the GUROBI solver.
%   * gurobi_prekernel                    - Computes a prekernel point using the GUROBI solver. 
%   * gurobi_prenucl                      - Computes the prenucleolus using the GUROBI solver.
%   * hsl_prekernel                       - Computes a prekernel point using HSL solvers. 
%   * ipopt_kernel                        - Computes a kernel point using the IPOPT solver. 
%   * ipopt_prekernel                     - Computes a prekernel point using the IPOPT solver.
%   * Kernel                              - Computes a kernel point using optimization toolbox.
%   * kernelQ                             - Checks if an imputation is a kernel point.
%   * k_Reduced_game_propertyQ            - Checks the k-RGP.
%   * LeastCore                           - Computes the leastcore using optimization toolbox.
%   * lin_prekernel                       - Computes a prekernel point using optimization toolbox.
%   * msk_kernel                          - Computes a kernel point using the MOSEK solver.
%   * msk_nucl                            - Computes the nucleolus using the MOSEK solver.
%   * msk_prekernel                       - Computes a prekernel point using the MOSEK solver.
%   * msk_prenucl                         - Computes the prenucleolus using the MOSEK solver.
%   * nucl                                - Computes the nucleolus using optimization toolbox.
%   * oases_kernel                        - Computes a kernel point using the OASES solver.
%   * oases_prekernel                     - Computes a prekernel point using the OASES solver.
%   * ols_prekernel                       - Computes a prekernel point using optimization toolbox.
%   * PreNucl2                            - Computes the prenucleolus using optimization toolbox.
%   * PreNucl                             - Computes the prenucleolus using optimization toolbox.
%   * qpBB_kernel                         - Computes a kernel point using the QPBB solver.
%   * qpc_kernel                          - Computes a kernel point using the QPC solver.
%   * qpc_prekernel                       - Computes a prekernel point using the QPC solver.
%
% Class Objects
% -------------
%   * TuGame                              - to perform several computations for retrieving and modifying game data.
%   * TuProp                              - subclass object of TuGame (game properties).
%   * TuSol                               - subclass object of TuGame (game solutions).
%   * TuCore                              - subclass object of TuSol (core).
%   * TuACore                             - subclass object of TuSol (anti-core).
%   * TuVert                              - subclass object of TuSol (core vertices).
%   * TuRep                               - subclass object of TuSol (prk replication).
%   * TuCons                              - subclass object of TuSol (consistency).
%   * TuKcons                             - subclass object of TuSol (generalized consistency).
%
%
%  PARALLEL COMPUTING
%   * p_Anti_PrekernelQ                   - Checks if an imputation is an anti prekernel point.
%   * p_clp_kernel                        - Computes a kernel point using the CLP solver.
%   * p_cls_kernel                        - Computes a kernel point using the CLS solver.
%   * p_coreQ                             - Checks the non-emptiness of the core.
%   * p_cplex_kernel                      - Computes a kernel point using the CPLEX solver.
%   * p_cplex_prekernel                   - Computes a prekernel point using the CPLEX solver.
%   * p_cvx_kernel                        - Computes a kernel point using the CVX solver. 
%   * p_cvx_prekernel                     - Computes a prekernel point using the CVX solver.
%   * p_glpk_kernel                       - Computes a kernel point using the GLPK solver. 
%   * p_glpk_prekernel                    - Computes a prekernel point using the GLPK solver. 
%   * p_gurobi_kernel                     - Computes a kernel point using the GUROBI solver. 
%   * p_gurobi_prekernel                  - Computes a prekernel point using the GUROBI solver. 
%   * p_hsl_prekernel                     - Computes a prekernel point using HSL solvers. 
%   * p_ipopt_kernel                      - Computes a kernel point using the IPOPT solver.
%   * p_ipopt_prekernel                   - Computes a prekernel point using the IPOPT solver.
%   * p_k_Reduced_game_propertyQ          - Checks the k-RGP.
%   * p_lin_prekernel                     - Computes a prekernel point using optimization toolbox.
%   * p_msk_kernel                        - Computes a kernel point using the MOSEK solver.
%   * p_msk_prekernel                     - Computes a prekernel point using the MOSEK solver.
%   * p_oases_kernel                      - Computes a kernel point using the OASES solver.
%   * p_oases_prekernel                   - Computes a prekernel point using the OASES solver.
%   * p_ols_prekernel                     - Computes a prekernel point using optimization toolbox.
%   * p_qpBB_kernel                       - Computes a kernel point using the QPBB solver.
%   * p_qpc_kernel                        - Computes a kernel point using the QPC solver.
%   * p_qpc_prekernel                     - Computes a prekernel point using the QPC solver.
%
% Class Objects
% -------------
%   * p_TuProp                            - subclass object of TuGame (game properties).
%   * p_TuSol                             - subclass object of TuGame (game solutions).
%   * p_TuRep                             - subclass object of p_TuSol (prk replication).
%   * p_TuCons                            - subclass object of p_TuSol (consistency).
%   * p_TuKcons                           - subclass object of p_TuSol (generalized consistency).
%
%
% BUG FIXES
%   * Converse_RGP_Q                 - Fix closed loop under method 'PRN', and 'HMS_PN'. 
%   * Reconfirmation_propertyQ       - Fix closed loop under method 'PRN', and 'HMS_PN'. 
%   * Reduced_game_propertyQ         - Fix incorrect computation and closed loop under methods 'PRN', and 'HMS_PN'. 
%   * StrConverse_RGP_Q              - Fix closed loop under method 'PRN', and 'HMS_PN'. 
%   * k_Converse_RGP_Q               - Fix incorrect computation and closed loop under methods 'PRN', and 'HMS_PN'. 
%   * k_Reconfirmation_propertyQ     - Fix closed loop under method 'PRN', and 'HMS_PN'. . 
%   * k_StrConverse_RGP_Q            - Fix incorrect computation and closed loop under methods 'PRN', and 'HMS_PN'. 
%   * p_Converse_RGP_Q               - Fix closed loop under method 'PRN', and 'HMS_PN'. 
%   * p_Reconfirmation_propertyQ     - Fix closed loop under method 'PRN', and 'HMS_PN'. 
%   * p_Reduced_game_propertyQ       - Fix incorrect computation and closed loop under methods 'PRN', and 'HMS_PN'. 
%   * p_StrConverse_RGP_Q            - Fix closed loop under method 'PRN', and 'HMS_PN'. 
%   * p_k_Converse_RGP_Q             - Fix incorrect computation and closed loop under methods 'PRN', and 'HMS_PN'. 
%   * p_k_Reconfirmation_propertyQ   - Fix closed loop under method 'PRN', and 'HMS_PN'. . 
%   * p_k_StrConverse_RGP_Q          - Fix incorrect computation and closed loop under methods 'PRN', and 'HMS_PN'. 
%   
%
% MODIFIED FUNCTIONS
% 
% SERIAL
% AllMarginalContributions,Anti_PreKernel,average_convexQ,bankruptcy_game,banzhaf,BestCoalitions,
% cardinality_game,CddCoreQ,CddCoreVertices,CddImputationVertices,CddPreKernel,Converse_RGP_Q,
% convex_gameQ,coreQ,critical_value1,critical_value2,critical_value_star,DM_Reduced_game,excess,
% game_basis,game_space,gameToMama,gameToMatlab,Gap,getgame,HMS_RedGame,HMS_Reduced_game,
% homogeneous_representationQ,ImputationVertices,k_Converse_RGP_Q,k_convexQ,k_cover,
% k_Reconfirmation_propertyQ,k_StrConverse_RGP_Q,minimal_winning,monotone_gameQ,RedGame,PreKernel,
% PrekernelQ,Reconfirmation_propertyQ,reasonable_outcome,savings_game,select_starting_pt,
% semi_convexQ,separable_cost_allocation,ShapleyValue,ShapleyValueM,SortMg,sortsets,
% StrConverse_RGP_Q,SubGame,SubSets,Sum_Marg_Contributions,super_additiveQ,TauValue,
% unanimity_games,veto_players,weightedShapley,zero_normalization;
%
% PARALLEL 
% p_AllMarginalContributions,p_Anti_PreKernel,p_assignment_game,p_average_convexQ,p_banzhaf,
% p_BestCoalitions,p_Converse_RGP_Q,p_convex_gameQ,p_coreQ,p_DM_Reduced_game,p_excess,
% p_game_basis,p_Gap,p_getgame,p_HMS_Reduced_game,p_homogeneous_representationQ,
% p_k_Converse_RGP_Q,p_k_convexQ,p_k_cover,p_k_Reconfirmation_propertyQ,p_k_StrConverse_RGP_Q,
% p_minimal_winning,p_monotone_gameQ,p_PreKernel,p_PrekernelQ,p_reasonable_outcome,
% p_Reconfirmation_propertyQ,p_RedGame,p_Reduced_game_propertyQ,p_replicate_prk,
% p_select_starting_pt,p_semi_convexQ,p_ShapleyValueM,p_StrConverse_RGP_Q,p_SubSets,
% p_TauValue,p_unanimity_games,p_weighted_majority
%
%
% About Version 0.2
% -----------------
% Release Date: 08-Oct-2012
%
% SUMMARY
% Code optimization for serial as well as parallel computing. Performance
% gain is about 40 percent in serial, whereas in parallel mode we observe
% a performance gain up to 60 percent and even more. The manual has been 
% revised and extended. Some new functions have been added and some bugs 
% have been fixed, see below for the details.
%
% NEW FUNCTIONS
%
%  SERIAL COMPUTING
%   * Anti_PreKernel               - Computes an anti-pre-kernel point.
%   * clToMatlab                   - Computes the unique integer representation of coalitions.
%   * critical_value1              - Computes the biggest gain of any group of players.
%   * critical_value2              - Computes a critical value w.r.t. the strong epsilon-core.
%   * critical_value_star          - Computes a critical value which contains the intersection of the imputation and reasonable set.
%   * homogeneous_representationQ  - Checks if the weighted majority game possesses a homogeneous representation.
%   * minimal_winning              - Computes the minimal winning coalitions.
%   * reasonable_outcome           - Determines the reasonable outcome.
%   * select_starting_pt           - Selects a starting point for the pre-kernel computation.
%   * vclToMatlab                  - Computes a Tu-game and the corresponding unique integer representation of coalitions
%   * winning_coalitions           - Determines the whole set of winning coalitions.
%
%
%  PARALLEL COMPUTING
%   * p_Anti_PreKernel                  - Computes an anti-pre-kernel element. 
%   * p_excess                          - Computes the excesses.
%   * p_homogeneous_representationQ     - Checks if the weighted majority game possesses a homogeneous representation.
%   * p_minimal_winning                 - Computes the minimal winning coalitions.
%   * p_reasonable_outcome              - Determines the reasonable outcome.
%   * p_select_starting_pt              - Selects a starting point for the pre-kernel computation.
%
%
% BUG FIXES
%   * BestCoalitions               -  Fix incorrect computation with null-vector as starting point. 
%   * belongToCoreQ                -  Fix incorrect computation with non-efficient payoff vectors.
%   * coreQ                        -  Fix incorrect computation with non-efficient payoff vectors.
%   * PreKernel                    -  Fix incorrect computation with null-vector as starting point. 
%   * PreKernelQ                   -  Fix incorrect computation with null-vector as starting point. 
%   * p_BestCoalitions             -  Fix incorrect computation with null-vector as starting point. 
%   * p_k_Converse_RGP_Q           -  Fix a problem with undefined variables.
%   * p_k_Reconfirmation_propertyQ -  Fix a problem with undefined variables.
%   * p_k_StrConverse_RGP_Q        -  Fix a problem with undefined variables.
%   * p_PreKernel                  -  Fix incorrect computation with null-vector as starting point. 
%   * p_PrekernelQ                 -  Fix incorrect computation with null-vector as starting point. 
%
%
% MODIFIED FUNCTIONS
% These functions have now a modified behavior
%   * coreQ                     - A tolerance value can be provided as input. 
%                               - Two input arguments are allowed. 
%   * CorePlot                  - Five input arguments are allowed. 
%                               - A tolerance value can be provided as input. 
%   * CoreVertices              - A tolerance value can be provided as input.
%                               - Two input arguments are allowed. 
%   * CddCorePlot               - Only four input arguments are allowed. 
%                               - Superfluous methods have been removed. 
%                               - A tolerance value can be provided.
%   * CddCoreVertices           - A tolerance value can be provided as input.
%                               - Two input arguments are allowed. 
%   * CddImputationVertices     - Superfluous methods have been removed.
%                               - Only one input argument is allowed.
%
%  Almost all MatTuGames Toolbox functions were modified to clean up or simplify
%  the code. Exceptions are:
%  AllMarginalContributions,AllSubGames,cardinality_game,dual_game,game_space,getgame,harsanyi_dividends,
%  profit_matrix,replicate_prk,ShapleyValueM,SortMg,sortsets,StandardSolution,SubDual,SubGame,
%  Talmudic_Rule,weakly_super_additiveQ,zero_monotonicQ,
%  p_AllMarginalContributions,p_game_space,p_game_space_red,p_getgame,p_harsanyi_dividends,
%  p_ShapleyValueM,p_super_additiveQ,p_zero_monotonicQ, 
%  and all functions in the mama folder.
%
%
% DELETED FUNCTIONS
%   * qhull2mama                 - Superfluous now.
%   * p_bestCoalitions           - Removed, no longer used by the toolbox.
%   * tug_AvConvexCoord          - Replaced by tug_DetUCoord.
%   * tug_CddGmpVerticesCore     - Removed, no longer used by the toolbox.
%
% 
% NEW SCRIPT
%   * getting_started            - Checks the installation.
%   * ReleaseNote                - Invokes MatTuGames_Version_0.2.
% 
% About Version 0.1
% -------------------
% Release date: 03-Apr-2012
%
%
% Updated installation instruction
% Revised documentation
% A boundary value corrected in CddPreKernel.m
%
%
%
