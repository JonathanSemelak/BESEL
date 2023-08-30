![plot](./img1.png)

# FENEB
A tool for performing nudged elastic bands simulations on the free energy surface using Amber.

# How it works
Run molecular dynamics with AMBER, feed this code with .rst7 and .nc files, and you will get new .rst7 files to continue the optimization. Additionally, every time this code is executed, an output file will be generated with thermodynamic information.

![plot](./img2.png)

Using a bash script is recommended for automatization. For this reason, a feneb_wizard is provided, which helps you running in parallel the required MD simulations, executing feneb, and iterating until the desired convergence criteria is reached.

The feneb code uses an input file "feneb.in". Below and example feneb.in file:


```
 prefix my_system    !Prefix of .prmtop .nc and .rst7 files
 nrep 25             !Number of replicas (including the extremes)
 tangoption 1        !Tangent option definition
 ksring 200          !Spring force constant (kcal/mol A**2)
 kref 500            !Bias force constant (kcal/mol A**2) for V(q)=k/2(q-qref)**2
 nrestr 5            !Number of atoms forming the Reaction Coordinate Space
 mask 1 2 5-7        !Indexes of the atoms forming the Reaction Coordinate Space
 skip 1000           !Number of initial frames ommited for analysis
```

Below all the variables available for feneb.in with their default options:

```
 prefix                  !Prefix of .prmtop .nc and .rst7 files. No default for this variable. Must be defined.
 nrep                    !Number of replicas (including the extremes). No default for this variable. Must be defined.
 per                   T !System is periodic
 velin                 F !Input .rst7 files have velocity info
 velout                F !Include dummy velocity info for output files
 nscycle           50000 !Max. optimization steps in string optimization
 rextrema              F !Use info from previous extremes optimization (requires feneb.reactants and feneb.products)
 steep_spring    0.001d0 !Step length (in A)  for string optimization
 steep_size       0.01d0 !Step length (in A)  for steepest descent optimization
 smartstep             T !Ignore steep_size and use a gradient based criteria for deciding the step lenght
 skip                  0 !Number of initial frames ommited for analysis
 typicalneb            F !Instead of an uncoupled feneb optimization, run a feneb optimization with the NEB force
 usensteps             F !Use only a certain number of steps from the trajectory files
 nstepsexternal     5000 !Specifiy how many steps will be used
 wtemp                 F !Write temperature for the atoms in the Reaction Coordinate Space
 wtempfrec             1 !Frequency of temperature writing
 dt                0.001 !Time step in the trajectory file (depends on how often the coordinates were written to the .nc file!)
 dostat                F !Make some stats to check if the free energy gradient reaches a stationary value (Mann-Kendall test)
 minsegmentlenght    100 !Has to do with the Mann-Kendall test. Only modify if you know what you're doing.
 nevalfluc          1000 !Number of minimum frames to evaluate fluctuations (Has to do with the Mann-Kendall test)
 tangoption            1 !Tangent option definition: 0, classic; 1, improved; 2, free energy based tangent (uses tg+ and tg-)
 tangrecalc            T !Recalculate the tangent every uncoupled string optimization step
 maxdist         0.001d0 !Criteria for the replicas to be "equispaced" (A)
 stopifconverged       F !Stop if the convergence criteria is achived (if so, no new band is generated)
 wgrad                 F !Write the free energy gradient of each coordinate of the Reaction Coordinate Space
 ftol               2.25 !Convergence criteria in terms of maximum root squeare for the perpendicular free energy gradient (kcal/mol AA)

```

If you want to use this code for a full optimization on the free energy surface, set nrep=1. All keywords related with NEB will be ignored in this case.

It is strongly recommended to use the feneb_wizard (and feopt_wizard). This will automatice the whole procedure. To learn how to use it, just type:

```
feneb_wizard -h

```

Or :

```
feopt_wizard -h

```

Check the tutorial folder for examples and a short tutorial!

# Requirements
This code requires NETCDF Libraries.

# Compilation
Execute

```
  make
```

NETCDF_DIR is assumed to be /usr. If not, override with:

```
  make NETCDF_DIR=/your/path/to/netcdf`
```

After installation, export PATH=$PATH:/your/path/to/feneb/bin

# References

- Semelak, et. al. (2023). Minimum Free Energy Pathways of Reactive Processes with Nudged Elastic Bands. Journal of Chemical Theory and Computation.
- Bohner, et. al. (2014). Nudged-elastic band used to find reaction coordinates based on the free energy. The Journal of Chemical Physics, 140(7), 074109.
