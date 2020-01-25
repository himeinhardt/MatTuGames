disp('Checking Basic Installation of MatTuGames.'); 
disp('The claims vector is specified by:');
d=[40 32 11 73.3 54.95 81.1]
disp('and the estate by:');
E=ceil((3/5)*sum(d))
disp('The Talmudic distribution is given by:');
disp('tlm_rl=Talmudic_Rule(E,d)');
tlm_rl=Talmudic_Rule(E,d)
disp('The corresponding bankruptcy game is given by:');
disp('bv=bankruptcy_game(E,d)');
bv=bankruptcy_game(E,d)
disp('Is the game convex?');
disp('cvQ=convex_gameQ(bv)');
cvQ=convex_gameQ(bv)
disp('Is the core non-empty?');
disp('crQ=coreQ(bv)');
crQ=coreQ(bv)
disp('Is the game monotone?');
disp('mQ=monotone_gameQ(bv)');
mQ=monotone_gameQ(bv)
disp('Is the game zero-monotone?');
disp('zmQ=zero_monotonicQ(bv)');
zmQ=zero_monotonicQ(bv)
disp('Is the game average convex?');
disp('acvQ=average_convexQ(bv)');
acvQ=average_convexQ(bv)
disp('Is the game super additive?');
disp('sadQ=super_additiveQ(bv)');
sadQ=super_additiveQ(bv)
disp('Is the game semi convex?');
disp('scvQ=semi_convexQ(bv)');
scvQ=semi_convexQ(bv)
disp('The Harsanyi dividends are given by:');
disp('hd=harsanyi_dividends(bv)');
hd=harsanyi_dividends(bv)
disp('We get the following game from the Harsanyi dividends:');
disp('v=getgame(hd)');
v=getgame(hd)
tol=10^8*eps;
disp('Coincides this game with the original game bv?');
geqQ1=all(abs(v-bv)<tol)
disp('The Shapley value of the game is:');
disp('sh_v=ShapleyValue(bv)');
sh_v=ShapleyValue(bv)
disp('The Tau value of the game is:');
disp('tau_v=TauValue(bv)');
tau_v=TauValue(bv)
disp('The solidarity value of the game is:');
disp('sl_vl=SolidarityValue(bv)');
sl_vl=SolidarityValue(bv)
disp('The pre-kernel element of the game is:');
disp('prk_v=PreKernel(bv)');
prk_v=PreKernel(bv)
disp('The kernel element of the game is:');
disp('kr_v=Kernel(bv)');
kr_v=Kernel(bv)
disp('The pre-nucleolus of the game is:');
disp('prn_v=PreNucl(bv)');
prn_v=PreNucl(bv)
disp('Is the solution found a pre-nucleolus?:');
disp('prnQ_v=PrenuclQ(bv,prn_v)');
prnQ_v=PrenuclQ(bv,prn_v)
disp('Checking Property II (Kohlberg). Is the induced collection balanced?:');
disp('bcQ_v=balancedCollectionQ(bv,prn_v)');
bcQ=balancedCollectionQ(bv,prn_v)
disp('The nucleolus of the game is:');
disp('nuc_v=nucl(bv)');
nuc_v=nucl(bv)
disp('The vector of excesses w.r.t. the pre-kernel is:');
disp('ex_prk=excess(bv,prk_v)');
ex_prk=excess(bv,prk_v)
disp('Is the vector prk_v a pre-kernel element?');
disp('prkQ=PrekernelQ(bv,prk_v)');
prkQ=PrekernelQ(bv,prk_v)
disp('The anti-pre-kernel element of the game is given by:');
disp('aprk=Anti_PreKernel(bv)');
aprk=Anti_PreKernel(bv)
disp('The dual game is:');
disp('dv=dual_game(bv)');
dv=dual_game(bv)
disp('The greedy bankruptcy game is specified by:');
disp('gv=greedy_bankruptcy(E,d)');
gv=greedy_bankruptcy(E,d)
disp('Is the greedy game equal to the dual?');
geqQ2=all(abs(gv-dv)<tol)
disp('A pre-kernel element of the dual is given by:');
disp('prk_dv=PreKernel(dv)');
prk_dv=PreKernel(dv)
disp('An anti-pre-kernel element of the dual is given by:');
disp('aprk_dv=Anti_PreKernel(dv)');
aprk_dv=Anti_PreKernel(dv)
disp('Satisfies the pre-kernel consistency?');
% disp('[RGP RGPC]=Reduced_game_propertyQ(bv,prk_v,PRK)');
echo on
[RGP RGPC]=Reduced_game_propertyQ(bv,prk_v,'PRK');
echo off;
RGP
disp('Satisfies the pre-kernel converse consistency?');
echo on; 
[CRGP CRGPC]=Converse_RGP_Q(bv,prk_v,'PRK'); 
echo off;
CRGP
disp('Satisfies the pre-kernel element the reconfirmation property?');
echo on; 
[RCP RCPC]=Reconfirmation_propertyQ(bv,prk_v,'PRK'); 
echo off;
RCP
disp('Is the pre-kernel element replicable for some related games?');
disp('RepSol=replicate_prk(bv,prk_v,2,1)');
RepSol=replicate_prk(bv,prk_v,2,1);
RepSol
disp('Select a partition of the player set by:');
disp('P={[1 3],[2 4 5],[6]}');
P={[1 3],[2 4 5],[6]}
disp('Transcribe it to its unique integer representation through:');
disp('pm=clToMatlab(P)');
pm=clToMatlab(P)
disp('The Owen value w.r.t. a priori unions pm is:');
disp('ow_vl=OwenValue(bv,pm)');
ow_vl=OwenValue(bv,pm)
disp('The weighted Owen value w.r.t. a priori unions pm is:');
disp('wow_vl=weightedOwen(bv,pm)');
wow_vl=weightedOwen(v,pm)
disp('The coalition solidarity value w.r.t. a priori unions pm is:');
disp('csl_vl=CoalitionSolidarity(bv,pm)');
csl_vl=CoalitionSolidarity(bv,pm)
disp('Define a communication structure by:');
disp('CS={[1 2],[1 3],[1 4],[2 3],[2 4],[3 4],[4 5],[4 6],[5 6]}');
CS={[1 2],[1 3],[1 4],[2 3],[2 4],[3 4],[4 5],[4 6],[5 6]}
disp('Transform it to its unique integer representation by:');
disp('csm=clToMatlab(CS)');
csm=clToMatlab(CS)
disp(['The Myerson value w.r.t. communication structure csm is specified by:']);
disp('my_vl=MyersonValue(bv,csm)');
my_vl=MyersonValue(bv,csm,'cs')



