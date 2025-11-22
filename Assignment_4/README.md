# Assignment 4: J2 Plasticity with Radial Return

A MATLAB project demonstrating the implementation of a J2 plasticity model using a radial return algorithm for a single quadrilateral element.

## Background

This final assignment for the ME621 course was the culmination of the semester's work, focusing on a fundamental and widely used concept in computational solid mechanics: J2 plasticity. The goal was to implement a complete material model from scratch, including the crucial radial return algorithm for integrating the elastoplastic constitutive equations. This project provided a hands-on understanding of how non-linear material behavior is modeled in finite element simulations.

## Overview

This project analyzes the behavior of a single 4-node quadrilateral element under a prescribed displacement path. The material is modeled using J2 plasticity with isotropic hardening.

The core of the project is the implementation of the radial return algorithm within the main script, `Main.m`. This algorithm is an iterative process used to enforce the yield condition after an elastic trial stress is calculated. The script handles:

1.  **Initialization**: Defining material properties (Young's modulus, Poisson's ratio, yield stress, hardening modulus) and the prescribed displacement path.
2.  **Loading Loop**: Incrementally applying displacements.
3.  **Stress Update**: For each increment:
    *   Calculate an elastic trial stress.
    *   Check for yielding against the von Mises yield criterion.
    *   If yielding occurs, perform the **radial return** to bring the stress state back to the yield surface.
    *   Update the plastic strain and stress.
4.  **Results**: Plotting the resulting stress-strain curve and the stress path.

## Quick Start

To run the analysis, simply open and execute the `Main.m` script in MATLAB.

```bash
# Open Main.m in MATLAB and run it.
```

The script will run the simulation and generate the plots showing the material response.

## Results

The output of the simulation clearly demonstrates the expected elastoplastic behavior. The stress-strain curve shows the initial elastic region, the yield point, and the subsequent plastic hardening. The loading-unloading cycle plot illustrates the permanent deformation characteristic of plastic behavior.

![Stress vs. Strain Curve](latex/img/stress_vs_strain.png)
![Loading and Unloading Cycle](latex/img/loading_unloading_cycle.png)

## Future Improvements

-   **Full FEA Integration**: The current implementation is for a single element. A major extension would be to integrate this material model into a larger finite element code to analyze complex structures.
-   **Hardening Models**: The project uses simple isotropic hardening. It could be extended to include more advanced models like kinematic hardening or a combination of both.
-   **Performance**: For larger simulations, the efficiency of the radial return algorithm could be analyzed and potentially improved.

---

Have a nice coding day,
Tommaso :panda_face:
