# Execution order of the programs

For each folder, there is a README.md providing in which order the program should be executed. Here I provide only the order of the folders execution

1. **clim60sec**, to create the climate database
2. **createData**, to create the tree database (merge then with climate)
3. **growth**, to estimate the individual radial growth
4. **mortality**, to estimate the individual mortality
5. **little1971**, *nothing to execute there, the data are already formatted*
6. **createMatlabData**, to create the data for the calculus of $R_0$
7. **R0**, to calculate $R_0$
8. **randomForest**, to execute the sdm
9. **validation**, to calculate the correlation $R_0$ and presence
10. **time**, to calculate the time it takes to grow up to a certain height without competition. *This folder is optional, I did it as a posteriori to check that it takes an 'infinite' time to grow up to 45 m height*
11. **toolFunction.R** is a script containing regularly called functions and is therefore at the 'root folder' for ease

