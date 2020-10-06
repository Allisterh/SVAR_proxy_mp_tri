# proxy_mp_tri
 Code for proxy SVAR paper: **Herwartz, Rohloff and Wang (2020): Proxy SVAR identification of monetary policy shocks - Monte Carlo evidence and insights for the US**

## notes: 
1. For Monte Carlo simulations: A complete simulation with full sample size and iteration steps could take a while. Please use ``Showcase.R`` to run a small-scale Monte Carlo experiment. For replication of the results, please use ``Replica.R`` to extract the files from the folder Replica. These results are generated on the StatOek3 Server.
2. For application: In order to increase the readability of the code, we did not put ``set.seed`` all over the place. For replication of the results, please use ``set.seed(1234)``. The main file is ``main.R``.
3. All involved functions are collected in the local package ``Functions``. 
