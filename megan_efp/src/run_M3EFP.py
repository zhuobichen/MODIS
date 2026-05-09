# Import libraries, packages, and modules
from __future__ import division
import sqlite3
import os
import errno
import warnings
import pandas as pd
import numpy as np
import M3GEFP as mgefp
warnings.simplefilter(action="ignore", category=FutureWarning)


def make_dir(path):
    """

    Function to check if output path exists
    and create directories if needed
    :param path: directory path
    :return: new directory if one does not previously exists

    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def validate_EF(indir,EF,Ecotype_Crop, Ecotype_Shrub,Ecotype_Herb, Ecotype_Tree):
    df_ef = pd.read_csv(indir+EF)
    df_gf = pd.concat(map(pd.read_csv, [indir+Ecotype_Crop, 
                                        indir+Ecotype_Shrub, 
                                        indir+Ecotype_Herb, 
                                        indir+Ecotype_Tree]))
    plants = df_gf['VegID'].unique()
    all_plants = df_ef['VegID'].unique()
    missing_plants = [x for x in plants if x not in all_plants]
    if len(missing_plants)>0:
        print("There are missing emission factor data for plant species:")
        print('\n'.join(missing_plants))
        skip = input("Do you want to continue without setting them? The missing emission factors will be set to zero by default. (Y/N)")
        
        if skip in ['Y','y']: 
            for mp in missing_plants:
                df_mp = pd.Series(0, index=df_ef.columns)
                df_mp['VegID'] = mp
                df_ef = df_ef.append(df_mp, ignore_index=True)
                df_ef.to_csv(indir+EF+'tmp',index=False)
        else: 
            exit(1)
    else:
        df_ef.to_csv(indir+EF+'tmp',index=False)
    return EF+'tmp'

def m3efp_driver(scen_name, GridDB_tables, EF,
                 Ecotype_Crop, Ecotype_Herb, Ecotype_Tree,
                 Ecotype_Shrub, grid_ecotype, grid_growth_form,
                 M3GEFP_database_path, EFa, EFz, LDFa, LDFz, output_path):
    """
    Main function to aggregate user options and call functions to create databases and output
    :param scen_name: Scenario Name
    :param GridDB_tables: Path to the CSVs to import into the M3GEFP database
    :Ecotype_Crop: input csv file name
    :Ecotype_Herb: input csv file name
    :Ecotype_Tree: input csv file name
    :Ecotype_Shrub: input csv file name
    :grid_ecotype: input csv file name
    :grid_growth_form: input csv file name
    :param M3GEFP_database_path: Path to SQLite database to be created for M3GEFP DB
    :param TotalClasses: Total number of classes to loop through when generating output
    :param output_path: Path to output CSV

    :return: 1 Grid EF CSV with all classes
    """

    print("\n MEGAN EF Processor")
    # Print user specification
    print("\n\nUser Settings")
    print("Scenario Name: %s" % scen_name)
    print("M3GEFP Input Directory: %s" % GridDB_tables)
    print("M3GEFP Database: %s" % M3GEFP_database_path)
    print("CSV Output Directory: %s" % output_path)

    print("Creating grid EF database")

    try:
        os.remove(M3GEFP_database_path)
    except OSError:
        pass


    # Check EF data completeness
    EF_tmp = validate_EF(GridDB_tables,EF,Ecotype_Crop, Ecotype_Shrub,Ecotype_Herb, Ecotype_Tree)
    
    # Create base database and initialize functions
    M3GEFP_connection = sqlite3.connect(M3GEFP_database_path)
    M3GEFP_connection.text_factory = str

    # Load CSV files as tables into the base database
    mgefp.make_M3GEFP_tables(M3GEFP_connection, GridDB_tables, Ecotype_Crop, Ecotype_Shrub,
                    Ecotype_Herb, Ecotype_Tree, grid_ecotype, grid_growth_form,EF_tmp)
    print("M3GEFP Database created: %s" % M3GEFP_database_path)


        
    # Run M3GEF database driver
    mgefp.run_M3GEFP_DB(M3GEFP_connection, EFa, EFz, LDFa, LDFz, output_path+"OutputGridEF.%s.csv" % scen_name)
    
    # Close database connection
    M3GEFP_connection.close()
    print("\n M3GEFP Database connection closed\n\n")

    print("grid EF Output: "+output_path+"OutputGridEF.%s.csv" % scen_name)

    os.remove(GridDB_tables+EF_tmp)
