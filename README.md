# HornerK
Horner's method as accurate is if computed in k-fold precision and then rounded into the working precision

## Authors
*[Thomas R. Cameron](https://thomasrcameron.com)
	* Mathematics Department, Penn State Behrend
	* [Email: trc5475@psu.edu](mailto:trc5475@psu.edu)
	
* [Stef Graillat](stef.graillat@sorbonne-universite.fr)
	* LIP6, Sorbonne University
	* [Email: stef.graillat@sorbonne-universite.fr](mailto:stef.graillat@sorbonne-universite.fr)
	
## Instructions
Below are instructions for the installation of HornerK, the testing of Horner, and compiling of TeX figures. These instructions have been tested on macOS High Big Sur 11.7.8.

### Installation
First, open the make.inc file to specify the C compiler and flags. The default settings use the *gcc* compiler with flags *-O3*. Once these parameters are set, the tests in the C folder can be installed by running *make install test_name* in the terminal. The command *make uninstall* can be used to remove all executable files. 

### Run Tests
When running tests, all output will be placed in a csv file inside of the csv_files folder. 

### Compile TeX
The data stored in the csv files can be viewed by compiling the corresponding tex file in the TeX folder. 