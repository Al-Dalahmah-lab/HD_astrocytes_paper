# HD Astrocyte Paper code

Github repository for the code used in the manuscript Multi-OMIC analysis of Huntington disease reveals a neuroprotective astrocyte state. 

## Usage

The directory contains two folders, analysis and function. The analysis section contains portion of codes used for the various analysis performed in our study. The function folder contains helper functions for our analysis, but does not require any further usage other than calling it in your R environment to use. 

## Additional data
We have three main dataset which contain the original cells used in the snRNAseq analysis. They are available at this link  which takes you to a google drive[https://drive.google.com/drive/folders/1FB8IrtYNyJnaWh7X7q68m8uZBlWsZK4p?usp=share_link](). There are three objects, the "master_hd_object_paper.rds" which contains the original nuclei of all cell types before QC (as described in paper). "comb_hd_cd44_neg_paper.rds" contains only CD44 negative (aka protoplasmic) astrocytes and "comb_hd_cd44_pos_paper.rds" contain the CD44 positive (aka fibrous-like) astrocytes both after QC. 

## Contact
Fahad Paryani, fp2409@cumc.columbia.edu

## License

[MIT](https://choosealicense.com/licenses/mit/)

