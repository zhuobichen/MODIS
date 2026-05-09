# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 13:43:58 2021

@author: YSHI
"""

## validation

import pandas as pd

def validate_EF(EF,Ecotype_Crop, Ecotype_Shrub,Ecotype_Herb, Ecotype_Tree):
    df_ef = pd.read_csv(EF)
    df_gf = pd.concat(map(pd.read_csv, [Ecotype_Crop, Ecotype_Shrub, Ecotype_Herb, Ecotype_Tree]))
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
            df_ef.to_csv(EF+'tmp',index=False)
    else: 
        exit(1)
    
    return EF+'tmp'



if __name__ == '__main__':
    EF = r"\\rosie\disk46\aqrp_20007\MEGEFP32\inputs\EFP\EF_new_mtest.csv"
    Ecotype_crop = r"\\rosie\disk46\aqrp_20007\MEGEFP32\inputs\EFP\SpeciationCrop.csv"
    df_ef = pd.read_csv(EF)
    df_gf = pd.read_csv(Ecotype_crop)
    newEF = validate_EF(EF,Ecotype_crop)
    
    print('no exit')