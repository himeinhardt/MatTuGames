# *Matlab* Toolbox *MatTuGames* Version 1.7.0

## 1. Introduction

The game theoretical *Matlab* toolbox *MatTuGames* provides about 230 functions
for modeling, and calculating some solutions as well as properties of cooperative
games with transferable utilities. In contrast to existing Matlab toolboxes to
investigate TU-games, which are written in a C/C++ programming style with the consequence
that these functions are executed relatively slowly, we heavily relied on vectorized
constructs in order to write more efficient Matlab functions. In particular, the toolbox
provides functions to compute the (pre-)kernel, (pre-)nucleolus, and anti (pre-)kernel
as well as game values like the Banzhaf, Myerson, Owen, position, Shapley, solidarity,
and coalition solidarity value and much more. In addition, we will discuss how one can
use *Matlab's Parallel Computing Toolbox* in connection with this toolbox to benefit
from a gain in performance by launching supplementary Matlab workers. Some information
are provided how to call our *Mathematica* package *TuGames* within a running Matlab session.

## 2. Getting started:
In order to get some insight how to analyze a cooperative game,
a so-called transferable utility game with the Game Theory Toolbox
*MatTuGames*, we discuss a small example to demonstrate how one can
compute some game properties or solution concepts, like convexity,
the Shapley value, the (pre-)nucleolus or a pre-kernel element.

For this purpose, consider a situation where an estate is insufficient
to meet simultaneously all of the debts/claims of a set of claimants,
such a situation is known in game theory as a bankruptcy problem.
The problem is now to find a fair/stable distribution in the sense that
no claimant/creditor can find an argument to obstruct the proposed division
to satisfy at least partly the mutual inconsistent claims of the creditors.
In a first step, we define a bankruptcy situation while specifying
the debts vector and the estate that can be distributed to the
creditors. We restrict our example to a six-person bankruptcy problem
with a debts vector given by

```
>> d = [40.0000 32.0000 11.0000 73.3000 54.9500 81.1000];
```

and an estate value which is equal to

```
>> E = 176;
```

We observe immediately that the estate `E` is insufficient to meet all
of the claims simultaneously. It should be obvious that with these values
we do not have defined a cooperative game, however, these information
are enough to compute a proposal how to divide the estate between the
creditors. A fair division rule which is proposed by the Babylonian Talmud,
is given by

```
>> tlm_rl=Talmudic_Rule(E,d)
>>
tlm_rl =

20.0000 16.0000 5.5000 48.3500 30.0000 56.1500
```

However, this distribution rule does not incorporate the coalition formation
process. Thus, we might get a different outcome when we consider the
possibility that agents can form coalitions to better enforce their claims.
This means, we have to study the corresponding cooperative game. This can
be constructed while calling the following function

```
>> bv=bankruptcy_game(E,d);
```

Having generated a game, we can check some game properties like convexity

```
>> cvQ=convex_gameQ(bv)
>>
cvQ =

1
```

The returned logical value indicates that this game is indeed convex. This must
be the case for bankruptcy games. In addition, we can also verify if the
core of the game is non-empty or empty. To see this one needs just to invoke

```
>> crQ=coreQ(bv)
>> Optimization terminated.

crQ =

1
```

which is answered by affirmation. This result confirms our expectation, since each
convex game has a non-empty core.

After this short introduction of game properties, we turn our attention now
to some well known solution concepts from game theory. We start with the
Shapley value, which can be computed by

```
>> sh_v=ShapleyValue(bv)
>>
sh_v =

23.5175 18.7483 6.4950 44.3008 33.3317 49.6067
```

A pre-kernel element can be computed with the function

```
>> prk_v=PreKernel(bv)
>>
prk_v =

20.0000 16.0000 5.5000 48.3500 30.0000 56.1500
```

which must be identical to the distributional law of justice proposed by the Talmudic
rule. Moreover, it must also coincides with the nucleolus due to the convexity
of the game. To see this, let us compute first the nucleolus and in the next
step the pre-nucleolus

```
>> nc_bv=nucl(bv)

nc_bv =

20.0000 16.0000 5.5000 48.3500 30.0000 56.1500

>> pn_bv=PreNucl(bv)

pn_bv =

20.0000 16.0000 5.5000 48.3500 30.0000 56.1500
```

We observe that both solutions coincide, which must be the case for zero-monotonic games.
To check that these solutions are indeed the pre-nucleolus can be verified by Kohlberg's
criterion

```
>> balancedCollectionQ(bv,pn_bv)

ans =

1

>> balancedCollectionQ(bv,nc_bv)

ans =

1
```

Finally, to verify that the solution found is really a pre-kernel element can be done while typing

```
>> prkQ=PrekernelQ(bv,prk_v)
>>
prkQ =

1
```
For a deeper discussion of the function set provided by the toolbox consult the Manual
or type help mat_tug to get a short overview.

## 4. Custom Installation

To install the toolbox we recommend a custom installation while following the instructions that are given at the URL: 

[mltbx](https://mathworks.com/matlabcentral/answers/242430-how-do-i-install-a-mltbx-file-from-the-filesharing-site-into-r2015a)

A mltbx file is provided in the Release section. 

## 5. Acknowledgment


The author acknowledges support by the state of Baden-Württemberg through bwHPC.

Of course, the usual disclaimer applies.   

## 6. License

This project is licensed under the FreeBSD License - see the [LICENSE](LICENSE.md) file.

## Author

** Holger I. Meinhardt **
Institute of Operations Research
University of Karlsruhe (KIT) 
E-mail: Holger.Meinhardt@wiwi.uni-karlsruhe.de
        holger.meinhardt@partner.kit.edu
