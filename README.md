# MAMEBUS
Meridionally-Averaged Model of Eastern Boundary Upwelling Systems

For details pertaining to MAMEBUS, please refer to:
Moscoso, J.E, Stewart, A.L., Bianchi D., and McWilliams, J.C. A Meridionally Averaged Model of Eatern Boundary Upwelling Systems (MAMEBUSv1.0). Submitted to Geoscientific Model Development. 

Note:
The setup and analysis code uses cmocean for visualization, and assumes that 
the file is stored in the utils folder:
Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle, H. M., & DiMarco, S. F. (2016). True colors of oceanography. Oceanography, 29(3), 10.

To run the model:
1. Open the setparams.m code and update the model_code_dir to point to where the MAMEBUS code is stored on the local machine.

2. Run setparams.m code with a test function. This may looks like:
setparams('~/Desktop','test') 
to run a test file on the users local desktop

3. Open the terminal and change directories to the "test" folder
cd ~/Desktop/test

4. Build the model, and run. Enter the following into the terminal
sh Build_MAMEBUS.sh
sh run.sh

5. Use the plotSolution.m funciton to visualize the code at any time during the model run, for the test run, the command would be:
plotSolution('~/Desktop','test',true,2,5)
This specific command plots the temperature tracer (2) and averages the solution over the final five years of the model run.



Please contact Jordyn Moscoso: jmoscoso@atmos.ucla.edu with any questions
