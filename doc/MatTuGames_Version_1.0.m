% ChangeLog file for the MatTuGames Toolbox
%
%
%  Author:        Holger I. Meinhardt 
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%
%
% This release of MatTuGames was developed and tested using the Matlab
% R2018a release. A set of functions uses the Optimization Toolbox
% and the cdd-library by Komei Fukuda, which can be found under the URL:
%
% http://www.inf.ethz.ch/personal/fukudak/cdd_home/
%
% as well as the Matlab interface to the cdd solver (cddmex): 
%
% http://control.ee.ethz.ch/~hybrid/cdd.php.
%
% Alternatively on can also install the MPT3 toolbox that can be
% downloaded from 
% 
% http://control.ee.ethz.ch/~mpt/3/Main/Installation
%
% which includes the Cdd library, and to get full scope of
% operation of the graphical features.
%
% To run the toolbox even in parallel mode Matlab's Parallel Computing 
% Toolbox is needed. For connecting the Mathematica Package TuGames
% where Version 2.4 is available at the URL:
%
% http://library.wolfram.com/infocenter/MathSource/5709/
%
% the Mathematica Symbolic Toolbox is required, which can be found under
% the URL:
%
% http://www.mathworks.com/matlabcentral/fileexchange/6044-mathematica-symbolic-toolbox-for-matlab-version-2-0
%
% The most recent version can be made available by the author of this
% toolbox upon request.
%
% The MatTuGames Toolbox should work with all platforms. Moreover, the 
% Toolbox works also with the game theory toolbox written by Jean Derks. 
% This toolbox can be downloaded from:
%
% http://www.personeel.unimaas.nl/Jean-Derks/downlds/CGinMatlab20050731.zip
%
% In order to get full operation of the consistency functions we
% provide a set of adjusted and optimized files of the Derks
% toolbox. This set of files have been adjusted by us to conduct a 
% proper and more rapid consistency investigation. It also fixes a 
% problem with closed loops under certain game classes. This set of
% files can be made available upon request. Alternatively, one can 
% rely on the PreNucl() function, but this requires a license of 
% Matlab's optimization toolbox.
%
%
% Installation
% ------------
% To install the MatTuGames Toolbox, unzip the zip-file mat_tugV0d8b.zip, 
% and place the folder containing the functions on a local hard drive or 
% a network drive accessible to your computer. In the next step rename 
% the folder mat_tugV0d8b to mat_tug before including the folder location 
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
% Copyright (c) 2012-18, Holger I. Meinhardt
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
% The MatTuGames Toolbox was not created by the MathWorks Inc., nor is it 
% supported by MathWorks Inc.
%
%
% About Version 1.0
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
% As a consequence, these functions are not anymore backward compatible.
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
%   * MLextension                         - Computes the multi-linear extension.
%   * PermutationGame                     - Computes from an assignment matrix the permutation game.
%   * Potential                           - Determines the potential of a TU game (recursive).
%   * PrenuclQ                            - Checks if a payoff vector is the pre-nucleolus using Kohlberg's criterion.
%   * ShapleyValueML                      - Computes the Shapley value using multi-linear extension. 
%
%
%  PARALLEL COMPUTING
%   * p_holler                            - Computes the Holler index.
%   * p_PermutationGame                   - Computes from an assignment matrix the permutation game.
%   * p_potential                         - Determines the potential of a TU game (basis).
%   * p_ShapleyValue                      - Computes the Shapley value (potential).
%
%
%  MATHEMATICA SYMBOLIC TOOLBOX FUNCTIONS
%
%   * tug_DetRandCoord                    - Returns random unanimity coordinates.
%   * tug_MLExtension                     - Computes the multi-linear extension.
%   * tug_ShapleyValueML                  - Determines the Shapley value using multi-linear extension.
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


