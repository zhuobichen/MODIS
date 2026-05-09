'''
Author: longshicheng longshicheng@hycx-gd.cn
Date: 2025-11-14 12:34:22
LastEditors: longshicheng longshicheng@hycx-gd.cn
LastEditTime: 2025-11-17 01:39:50
FilePath: /shixiansheng/MODIS/MEGEFP32/MEGAN_EFP.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
""" MEGAN EF Processor

Description:
    The MEGAN Emission Factor Processor (MEGAN EFP) can be used to integrate growth form, ecotype, 
    plant species composition and species-specific emission factor data to generate landscape average 
    emission factors and light dependent factors for your model domain. 
    
    This program is built on Python (version 3.6) and SQLite.

Required Python Packages:
    pandas
    numpy

Local Modules (located in ./src/):
    run_M3EFP
    M3GEFP

Usage:
    $ python MEGAN_EFP.py

Development:
    This program was developed by Alex Guenther (alex.guenther@uci.edu) in ACCESS VBA  
    and then converted to python + SQLite by Andy Wentland (Ramboll) in 2017
    Revised by Alex Guenther in July 2019 (This is Version 1.9)    
    Updated by Yuge Shi in May 2021 (Version 3.2)

"""
# Import libraries and modules
import sys
sys.path.append("./src")
import run_M3EFP as efp



# ~~~~~~~~~~~~~ BEGIN USER OPTIONS ~~~~~~~~~~~~~ #

# ~~~~ SCENARIO NAME ~~~~ #

scen_name = "GD_cn27"

# ~~~~ DATABASE ~~~~ #
# Specify path and name of M3VTEF, M3GEFP, M3LDF databases to be created
M3GEFP_database = "./database/M3GEFP_database."+scen_name+".db"


# ~~~~ INPUT TABLES ~~~~ #
# Directories with premade CSV files
# to be converted into database tables
GridDB_tables = "./inputs/EFP/"

# input file names
# Under GridDB_tables
EF = "EFv210806.csv"
Ecotype_Crop = "SpeciationCrop210806.csv"
Ecotype_Herb = "SpeciationHerb210806.csv"
Ecotype_Shrub = "SpeciationShrub210806.csv"
Ecotype_Tree = "SpeciationTree210725.csv"
grid_ecotype = "grid_ecotype." + scen_name +".csv"
grid_growth_form = "grid_growth_form." + scen_name + ".csv"

# ~~~~ SETTINGS ~~~~ #
# Total number of classes to loop through
# Default: EFClasses = 18, LDFClasses = 3-6
# If non-default number used, updates must be made in all submodules' concat_*_tables function
#TotalClasses = 20
EFClasses = 19
LDFClasses0 = 3
LDFClasses1 = 6


# ~~~~ OUTPUTS ~~~~ #
# Output directory
outputs_path = "./outputs/"

# ~~~~~~~~~~~~~ END USER OPTIONS ~~~~~~~~~~~~~ #

# Call m3efp_driver function in the run_M3EFP module to run program
if __name__ == "__main__":
    efp.m3efp_driver(scen_name, GridDB_tables,EF, Ecotype_Crop, Ecotype_Herb, Ecotype_Tree, Ecotype_Shrub, grid_ecotype, grid_growth_form, M3GEFP_database, 1,EFClasses,LDFClasses0,LDFClasses1, outputs_path)
