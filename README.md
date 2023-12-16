# Online Dynamic Learning and Robust Control for Nonlinear Systems

## Overview

This repository contains the implementation of an innovative approach for Online Dynamic Learning and Robust Control. The proposed method is implemented on two systems, namely CDPR1 and CDPR2, each corresponding to a specific case study.

## Simulation Structure

The simulation folder is organized into two main subfolders:

- **CDPR_with_8_Cables:**
  - This subfolder includes the code for the implementation of the proposed method on CDPR1, representing the first case study.
- **CDPR_with_12_Cables:**
  - This subfolder contains the code for the implementation of the proposed method on CDPR2, representing the second case study.

### CDPR_with_8_Cables Subfolder

For the CDPR_with_8_Cables, each subfolder (Scenario_1, Scenario_2, and Scenario_3) contains variations of the simulations with and without disturbances. Here is a breakdown of the subfolders:

- **Scenario_1 (with exponential trajectory desired):**
  - **Scenario_1_without_dis:**
    - `scenario_1_without_dis.m`: Main code for the project without disturbance.
    - `Ploting_file.m`: Common file for plotting results.
  - **Scenario_1_with_dis:**
    - `scenario_1_with_dis.m`: Main code for the project with disturbance.
    - `Ploting_file.m`: Common file for plotting results.

- **Scenario_2 (with sinusoidal trajectory desired):**
  - **Scenario_2_without_dis:**
    - `scenario_2_without_dis.m`: Main code for the project without disturbance.
    - `Ploting_file.m`: Common file for plotting results.
  - **Scenario_2_with_dis:**
    - `scenario_2_with_dis.m`: Main code for the project with disturbance.
    - `Ploting_file.m`: Common file for plotting results.

- **Scenario_3 (with cylindrical trajectory desired):**
  - **Scenario_3_without_dis:**
    - `scenario_3_without_dis.m`: Main code for the project without disturbance.
    - `Ploting_file.m`: Common file for plotting results.
  - **Scenario_3_with_dis:**
    - `scenario_3_with_dis.m`: Main code for the project with disturbance.
    - `Ploting_file.m`: Common file for plotting results.

### CDPR_with_12_Cables Subfolder

The structure of the CDPR_with_12_Cables subfolder mirrors that of CDPR_with_8_Cables with similar subfolders and files.

## Usage

All codes have been implemented in MATLAB. To run a specific scenario, navigate to the corresponding subfolder and execute the main code (e.g., `scenario_1_without_dis.m`). Ensure MATLAB is properly configured, and all dependencies are installed.

## Important Note

- Please refer to the associated documentation for additional details on the scenarios, methodologies, and the proposed control approach. The Scenario_1 and Scenario_3 (Scenario 2 in the report) have been reported in the document. The Scenario_2 was not reported in the document due to page limitations.

- All parts of the codes have been explained using comments in the `.m` files.

- The codes of the literature methods have been incorporated in the Literature folder.

