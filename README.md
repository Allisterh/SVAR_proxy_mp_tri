# proxy_mp_tri
 Code for proxy SVAR paper: [**Herwartz, Rohloff and Wang (2020): Proxy SVAR identification of monetary policy shocks - Monte Carlo evidence and insights for the US**](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3714542)

## Notes: 
### Code:
1. For Monte Carlo Simulation: 
A complete simulation with full sample size and iteration steps could take a while. 
Please consider using ``Code/Simulation/Showcase.R`` to run a small-scale Monte Carlo experiment. 
For replication of the results, please use ``Code/Simulation/Replica.R`` to extract the files from the folder ``Code/Simulation/Server``. 
These results are generated on the StatOek3 Server. Server information: 128 CPUs each Intel(R) Xeon(R) CPU E5-4660 v4 @ 2.20GHz.

2. For Application:
The main file is ``Code/Main.R``, which contains most of our calculations such as, inter alia, specifications and diagonastics of reduced-form VAR, identifications of structural VAR, computation of IRFs, boostrap inference.
Supplement files:
- ``Code/Fa.R`` performs factor-augumented VAR analysis.
- ``Code/Hist_decomp.R`` performs historical decompositions.
- ``Code/Volcker.R`` computes the cumulative comtributions of Vocker's monetary policy to disinflation.
- ``Code/barplot.R`` plots the cumulative comtributions of Vocker's disinflation as a beautiful barplot.
**Please always run ``Code/Main.R`` first, in order to run other supplement files.**

3. All involved functions are collected in the local package ``Code/Functions``. 

### Data:
1. ``data/USA_Tri.csv`` contains variables in the VAR system.
2. Folder ``data/Instruments` contains all employed monetary policy proxies.
3. ``data/Factors/fred-database_code/current.csv`` contains informational variables, which latent factors are extracted from.
