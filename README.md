Code for the paper "Query Complexity of Boolean Functions under Indefinite Causal Order" by Alastair A. Abbott, Mehdi Mhalla, and Pierre Pocreau

The Matlab project requires:
- [CVX](https://github.com/cvxr/CVX): a Matlab package for convex optimizaton. 
- [QETLab](https://qetlab.com/Main_Page): a Matlab toolbox for exploring quantum entanglement theory.
- [qsuperop](https://github.com/alastair-abbott/qsuperops): a Matlab package to work with quantum superoperators and their causal structure.

Scripts of the Matlab project:
- **try_4bit_functions.m**: performs an exhaustive search for a gap between 2-query FO-supermaps and general supermaps, on all 4-bit boolean functions.
- **compute_analytic_primal.m**: computes the minimum bounded error for 2-query FO-supermap and general supermaps on the Boolean function exhibiting the highest separation, that of truth table: 0001011011101001.
- **compute_analytic_dual.m**: dual formulation of the SDP computed in _compute_analytic_primal.m_.
- **extraction.m**: extracts rigorous lower and upper bounds by rationalising the solutions obtained from the SDPs of primal.m and dual.m.

The folders data, extraction and utility should be added to the Matlab path for the code to function.