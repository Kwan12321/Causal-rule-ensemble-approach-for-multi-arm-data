# Simulation and Real-Data Analysis Code

This folder contains R scripts used to generate data, implement comparison methods, run simulation studies, and conduct the real-data analysis corresponding to Sections 4–5 and the Supporting Information of the manuscript.

## Overview of Workflow ##

1. **Data generation**
   - `data_gen.R`, `data_gen_overlap.R`, and `data_gen_overlap2.R` 
     - Generate synthetic datasets under different overlap and design settings.
     - `data_gen.R` for main simulations; `data_gen_overlap.R` for additional simuation (scenario 1); `data_gen_overlap2.R` for additional simulation (scenario 2).

2. **Method implementation**
   - The scripts `sim_method_*.R` define the estimation procedures for each method (BART, random forest, XGBoost, and the proposed method).
   - `mcre_agenet.R` contains the **core function for the proposed method** (estimation procedure and some helpful functions).  
   - `imp.R` contains functions for **variable importance calculation** used in the proposed method.
     
3. **Main simulations (Section 4)**
   - The scripts `Sec4_simulation_*.R` run the simulation studies for each method and save the results.
   - `Sec4_Figure.R` creates the main figures in Section 4 from the saved simulation results.

4. **Real-data analysis (Section 5)**
   - `Sec5_real_data_analysis.R` performs the real-data analysis and produces the numerical results and figures reported in Section 5.

5. **Supporting Information**
   - `Support_info_sec2_simulation_setting1.R`, `Support_info_sec2_simulation_setting2.R`  
     - Run additional simulation settings (weak and moderate overlap scenarios) reported in Section 2 of the Supporting Information.
   - `Support_Info_sec1_Figure.R`, `Support_Info_sec1_Table.R`  
     - Generate figures and tables for Section 1 of the Supporting Information.
   - `Support_Info_sec2_Figure.R`, `Support_Info_sec2_Table.R`  
     - Generate figures and tables for Section 2 of the Supporting Information.
   - `Support_Info_sec3_Figure.R`, `Support_Info_sec3_Table.R`  
     - Run the parameter tuning for proposed method in real-data analysis and generate figures and tables for Section 3 of the Supporting Information.
   - `Support_info_sec4_Figure.R`, `Support_info_sec4_simulation.R`  
     - Run the additional simulations (conformal prediction) and create the corresponding figures for Section 4 of the Supporting Information.
