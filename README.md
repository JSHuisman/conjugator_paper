# Conjugator Paper

This repository contains the code and data behind the the manuscript:<br/>
[Jana S. Huisman, Fabienne Benz, Sarah J.N. Duxbury, J. Arjan G.M. de Visser, Alex R. Hall, Egil A.J. Fischer, Sebastian Bonhoeffer. Estimating the rate of plasmid transfer in liquid mating cultures. *BioRxiv* 2020.](https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1)

For related content, please see the 
- [R package](https://github.com/JSHuisman/conjugator)
- [Shiny app](https://github.com/JSHuisman/conjugator_shiny)


## Experimental Data

Time course assay in LB medium:
- Raw data of colony counts, dilution factors and cell densities are provided in the file ‘Time_course_DRT_raw_colonycounts_dilutionfactors.csv’. 
- Raw data of optical density measurements for monoculture growth rates in LB medium combined for all strains are provided in the file ‘Growth_OD_data_D_R_T_LBmedium.csv’. 
- *Processed data* used for conjugation rate estimation in the time course assay are provided in the file ‘Time_course_DRT.csv’. 

Full protocol in VL medium:
*Note: some of this data has been re-used as a small subset from [Duxbury et al., 2021](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2021.2027), and is available in the corresponding [Dryad repository](https://datadryad.org/stash/dataset/doi:10.5061/dryad.2jm63xsp7). The strains from Duxbury et al., 2021, relate to the strains used in the current study as follows: D = D’, R = R, T= T, ESBL375 = D.*
- Raw data of colony counts, dilution factors and cell densities are provided in the file ‘Full_protocol_DRT_TRT_raw_colonycounts_dilutionfactors.csv’. This data for the TRT assay has been re-used from the medium condition ‘1.0x_VL’ for t = 0, 4 in dataset ‘Dec19_raw_colony_counts_conjug_assay_20210326.csv’ ([Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.2jm63xsp7)). 
- Raw datasets of optical density data used as input for maximum growth rate estimation are provided on [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.2jm63xsp7) (the ‘1.0x_VL’ medium condition).
- *Processed data* used for conjugation rate estimation in the full protocol assay are provided in the file ‘Full_protocol_DRT_TRT.csv’. The growth rate values (psi values) for strains D, R and T presented in this dataset have been reused from the control medium condition ‘1.0x_VL’ in datasets ‘20201027_Growth_rates_D_R_T.csv’ and ‘20201027_Growth_rates_ESBL375.csv’ provided on [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.2jm63xsp7). The cell density values presented in ‘Full_protocol_DRT_TRT.csv’ have been reused as a small subset from the control medium condition ‘1.0x_VL’ for t = 0, 4 in dataset ‘Dec19_conjug_assay_cell_densities_t0_4_allreps.csv’.

The R script for logistic growth model fitting from OD data for the time course assay is provided as the file ‘Growth_Model_Logistic.R’. The R script for maximum growth rate fitting from OD data for the full protocol assay is provided as ‘calculateslope_new_growthrate.py’ in the Dryad repository. For the cell density values for t = 24 presented in ‘Full_protocol_DRT_TRT.csv’, two different dilution factors were used per strain per agar type for the two agar plate replicates. Rather than taking a mean estimate, we assumed a Poisson error structure of this colony count data (from file ‘Full_protocol_DRT_TRT_raw_colonycounts_dilutionfactors.csv’) and used negative log likelihood values to select the optimal cell density estimate within a small parameter space. The R script is provided as file ‘Poisson estimates.R’, from Mark Zwart (Department of Microbial Ecology, The Netherlands Institute of Ecology (NIOO-KNAW), Wageningen, The Netherlands). 

