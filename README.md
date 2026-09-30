# Student progression model

Fortran implementation of a stage-based student progression model. The
program calculates performance and graduation-related quantities, generates
data for different persistence scenarios, and estimates numerical bounds.

## Source files

- `progression_model.f90`: Core model, stage probabilities, persistence,
  graduation and dropout quantities.
- `funciones_4.f90`: Analytical functions for the limiting cases of full
  persistence and zero persistence.
- `optimizacion_1.f90`: Numerical optimization of the bounds for a specified
  performance rate.
- `operaciones_4.f90`: Generates simulation data, analytical curves and
  numerical bound files.
- `main_4.f90`: Main program and input parameters. Choose the calculation
  by changing `operation_case` near the beginning of this file.
- `experiments.py`: generates figures 3 and 4.

## Compile

A Fortran compiler such as gfortran is required. From the directory
containing the five source files, run:

gfortran -O2 -o progression_model.exe progression_model.f90 \
  funciones_4.f90 optimizacion_1.f90 operaciones_4.f90 main_4.f90

## Run

On Windows:

progression_model.exe

On macOS or Linux:

./progression_model.exe

## Available calculations

Set `operation_case` in `main_4.f90` to:

1. Sample performance rate and graduation-rate data.
2. Calculate analytical bounds for full persistence.
3. Calculate analytical bounds for zero persistence.
4. Generate a curve with homogeneous pass probabilities and persistence.
5. Estimate numerical bounds for a selected persistence profile.

The remaining parameters and output filenames are defined in the
INPUT PARAMETERS block of `main_4.f90`. Output is written to `.dat` files
in the directory from which the program is run.
