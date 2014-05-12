NB-QFT-learning
===============

Quantum Field Theory based model-independent learning of continuous
distributions; following Nemenman and Bialek, PRE, 2002

Contributors: Ilya Nemenman

(c) Ilya Nemenman, 2000-2014

The algorithm behind this code is described in

Nemenman, I. & Bialek, W. Occam factors and model independent
      Bayesian learning of continuous distributions. Physical Review E
      65, 026137 (2002).

The code itself was written in 2000-2001, with some edits in 2005, to
support the analysis in the PRE paper. It was only released in 2014
due to an oversight.

The software is written in Matlab.


Functions/fiiles in package:
	0. bcsn.m -- main routine; given the data and auxiliary
        parameters, performs optimization over the smoothness scale
        and find the a posteriori optimal probability distribution.
    1. qclass.m  --  Solving for classical solution  of the learning
        problem (as in BCS'96), given histogramed data and the chosen
        smoothness scale
	2. QNcorrelator.m --  Negative action of Q(x_i), i=1..N.  refer to
    	Bialek et.al. 1996; needs the classical solution, the data
    	samples, and the chosen smoothness scale.
    3. Rfunctdet.m -- calculates the log of R, functional determinant,
        for the given classical solution and the smoothness scale.
	4. action.m -- calculates the action for minimization over l, the
    	smoothness scale, given the data and the smoothness scale
	5. cdiffl.m and cdiffr.m -- cyclic differences functions
	6. histQ.m -- histogramming data as needed for the algorithm to
        work. 
