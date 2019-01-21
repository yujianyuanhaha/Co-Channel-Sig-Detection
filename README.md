Joint MAP detection of co-channel signals
author: Daniel M. Jakubisim, Jet Yu  
Contact: Jet Yu (jianyuan@vt.edu)  


# code dependency
1. matlab communication toolbox
    * student licence of matlab does not include extra communication toolbox, please contact __Brandon__ (rbrand7@vt.edu) for access to install department licenced matlab
    * for Mac OS user, follow steps below, and install version 2015b
        * Open Finder, open the Go menu and select Connect to Server.
        * In the Server Address input, type: smb://computing.ece.vt.edu/matlab
        * Log in with your ECE account creditials.
        * Run the install script: osx_VERSION.command 

2. [dimple](https://github.com/analog-garage/dimple)  
    * the toolbox for  probabilistic modeling, inference, and learning. 
    * install on matlab:
        * download source file [version 0.07](https://github.com/analog-garage/dimple/tree/release_0.07)
        * in matlab prompt
            ```
                cd <Path -to -Dimple >
                startup
            ```    
    * copy the missing `.jar` file at another hidden [repository](https://s3.amazonaws.com/files.dimple.probprog.org/dimple_v0_07.zip), to be exact, copy folder `/solver/lib` to the published one.
    * full manual: https://s3.amazonaws.com/files.dimple.probprog.org/DimpleUserManual_v0.07_MATLAB_API.pdf 




# run codes
1. setup matlab with toolbox
2. setup dimple
3. export external java path e.g. `export MATLAB_JAVA=/Library/Java/JavaVirtualMachines/jdk1.7.0_40.jdk/Contents/Home/jre/`
4. launch matlab from terminal e.g. `cd /Applications/MATLAB_R2015a.app/bin/`, `./matlab`
5. change path to working folder, in matlab prompt, `start.m`
6. run main code `sim_cochannel_damp04_06_10.m`


# Appendix
1. matlab call java function  
Matlab support call external java function. e.g. [hello word tutorial](https://www.mathworks.com/matlabcentral/answers/37185-cannot-call-java-class-from-matlab) 



2. version differen in matlab and operation system  
However, matlab got its own java version, that may not be same as operation system. e.g. Matlab 2015b apply internal java 1.7 (by `version -java` to check in prompt), while my Mac OS is 1.8 (by `java -version` in terminal). Thus, it is necessary to call external java to fix it. According to [solution](https://www.mathworks.com/matlabcentral/answers/103056-how-do-i-change-the-java-virtual-machine-jvm-that-matlab-is-using-on-macos).  
    1. export external java path e.g. `export MATLAB_JAVA=/Library/Java/JavaVirtualMachines/jdk1.7.0_40.jdk/Contents/Home/jre/`  
    2. launch matlab from terminal e.g. `cd /Applications/MATLAB_R2015a.app/bin/`, `./matlab`
    3. test the hellow word example as above.
