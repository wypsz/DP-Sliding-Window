# DPSW-Sketch: A Differentially Private Sketch Framework for Frequency Estimation over Sliding Windows
This repository hosts the implementation code for the "DPSW-Sketch: A Differentially Private Sketch Framework for Frequency Estimation over Sliding Windows" as introduced in our research paper. The codebase includes implementations for various mechanisms, including U-Sketch, PCC-MG, and BLMZ-Sketch, alongside the DPSW-Sketch.

## Dataset

The dataset required for running the experiments is expected to be in a text file format, with each line representing an item. An example dataset file is provided in the `dataset` directory. Ensure your dataset follows this format for compatibility.

## Mechanisms Implementation

The `Mechanism` directory contains the implementation code for the U-Sketch, PCC-MG, BLMZ-Sketch, and the DPSW-Sketch algorithms.

## Requirements

- CMake 3.16 or higher
- C++14 compatible compiler (GCC 4.9+, Clang 3.4+, MSVC 19.0+)

## Getting Started

Follow these steps to compile and run the DPSW-Sketch implementation:

1. **Uncomment the DPSW-Sketch Code Block**: In `main.cpp`, find and uncomment the section of code corresponding to the DPSW-Sketch execution.

2. **Create Build Directory**: Run `mkdir DPSW` to create a new directory for building the project.

3. **Navigate to Build Directory**: Change directory to the newly created build directory with `cd ./DPSW`.

4. **CMake Configuration**: From within the `DPSW` directory, run `cmake ..` to configure the project build with CMake.

5. **Compile the Code**: Execute `make -j 16` to compile the code using multiple cores for faster compilation.

6. **Return to Project Root**: Once the compilation is complete, navigate back to the project root directory with `cd ../`.

7. **Execute the Script**: Run `./DPSW.sh` to execute the DPSW-Sketch.

The execution steps for other algorithms are similar to those outlined above.

