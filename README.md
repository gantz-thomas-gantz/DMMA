# Final Project - Direct Meet-in-the-Middle Attack (DMMA)

## Overview

This project explores the **Direct Meet-in-the-Middle Attack** (DMMA) algorithm, which is a cryptographic attack method designed to break symmetric key ciphers. The project involves implementing and optimizing the algorithm using different parallelization strategies to improve computation time and scalability.

The project contains several versions of the algorithm, each building on the previous one by incorporating additional parallelization techniques, such as MPI, OpenMP, and AVX2.

### Versions:

1. **v1**: Minimal MPI usage to distribute computations across multiple nodes.
2. **v2**: Full MPI usage to fully parallelize the algorithm across multiple nodes.
3. **v3**: Full MPI and OpenMP usage for both distributed and shared-memory parallelism.
4. **v4**: Full MPI, OpenMP, and AVX2 usage.

### Directory Structure:

- `src/`: Contains the source code files.
  - `communication.h`: Header for MPI communication-related functions.
  - `mitm_original.c`: Original version of the Direct Meet-in-the-Middle Attack.
  - `mitm_version1.c`: Version 1 with minimal MPI.
  - `mitm_version2.c`: Version 2 with full MPI.
  - `mitm_version3.c`: Version 3 with MPI and OpenMP.
  - `mitm_version4.c`: Version 4 with MPI, OpenMP, and AVX2.
  - `utilities.h`: Header for utility functions, containing mostly functions for dynamic arrays.
  
- `build/`: Contains the compiled build files.
  - `original`, `version1`, `version2`, `version3`, `version4`: Directories for compiled versions.

- `scripts/`:
  - `weak_scaling.sh`: Script for running weak scaling experiments on one specific version.
  - `version_comparison.sh`: Script for running weak scaling experiments on multiple versions to compare execution times.
  - `run_stress_test.sh`: Script for running stress tests.

- `stress_test_output.txt`: Output of stress tests.
- `weak_scaling.txt`: Results from weak scaling experiments.
- `versions_weak_scaling.txt`: Results from weak scaling experiments on different versions.

### Running the Code:

#### Weak Scaling

To run a weak scaling experiment, use the following command when connected to a Grid5000 site:

```bash
oarsub -p <cluster> -l host=<maximum number of nodes>/core=<number of cores per node>,walltime=<maximum expected runtime> -O weak_scaling_console_output.txt -E weak_scaling_error.txt "./weak_scaling.sh"
```

The weak\_scaling.sh script is responsible for executing the different versions of the algorithm with a range of arguments. It runs the selected version doubling the amount of computing nodes as the problem size doubles to measure the performance.
To run it, the maximum amount of nodes determined by the script (MPI processes) should be used. Adjust the file accordinly before running. 

#### Stress Test
To run a stress test on a specific version of the algorithm, use:
```bash
oarsub -p <cluster> -l host=<number of nodes>/core=<number of cores per node>,walltime=<maximum expected runtime> -O stress_test_output.txt -E stress_test_error.txt "./run_stress_test.sh"
```
The run\_stress\_test.sh script runs the algorithm with selected arguments to evaluate the performance of a selected version under heavy computational load.
Adjust the file accordingly and run it with the amount of nodes equal to the amount of MPI processes chosen in the script.

### Build Instructions
From the root of the project:
**Create the build directory and make file, compile and execute**:

   ```bash
   cmake -B  build
   cmake --build build
   mpiexec <mpiexec arguments> build/<executable> <executable arguments>
	 ```

### Requirements


- MPI (Message Passing Interface for all but the original version)
- OpenMP (for parallel processing in the 3th and 4th versions)
- AVX2 (for SIMD instructions in the 4th version)
- CMake for building the project

