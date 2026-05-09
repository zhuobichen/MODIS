#!/bin/csh

foreach dom ( GD_27km )

/bigdata/MEGAN/MEGANv3.2/MEGAN_Prep/Code/prepmegan4cmaq_grwform.x < prepmegan4cmaq.growthform.inp | tee -a grwform.$dom.log
mv ./Guangdong2023/grid_growth_form.csv ./Guangdong2023/grid_growth_form.$dom.csv

/bigdata/MEGAN/MEGANv3.2/MEGAN_Prep/Code/prepmegan4cmaq_cantype.x < prepmegan4cmaq.cantype.inp | tee -a cantype.$dom.log
mv ./Guangdong2023/CT3.csv ./Guangdong2023/CT3.$dom.csv

/bigdata/MEGAN/MEGANv3.2/MEGAN_Prep/Code/prepmegan4cmaq_ecotype.x < prepmegan4cmaq.ecotype.inp | tee -a ecotype.$dom.log
mv ./Guangdong2023/grid_ecotype.csv ./Guangdong2023/grid_ecotype.$dom.csv

/bigdata/MEGAN/MEGANv3.2/MEGAN_Prep/Code/prepmegan4cmaq_arid.x < prepmegan4cmaq.soil.$dom.inp | tee -a arid.$dom.log
mv ./Guangdong2023/grid_arid.csv ./Guangdong2023/grid_arid.$dom.csv

/bigdata/MEGAN/MEGANv3.2/MEGAN_Prep/Code/prepmegan4cmaq_non_arid.x < prepmegan4cmaq.soil.$dom.inp | tee -a non_arid.$dom.log
mv ./Guangdong2023/grid_non_arid.csv ./Guangdong2023/grid_non_arid.$dom.csv

/bigdata/MEGAN/MEGANv3.2/MEGAN_Prep/Code/prepmegan4cmaq_landtype.x < prepmegan4cmaq.soil.$dom.inp | tee -a landtype.$dom.log
mv ./Guangdong2023/grid_LANDTYPE.csv ./Guangdong2023/grid_LANDTYPE.$dom.csv

/bigdata/MEGAN/MEGANv3.2/MEGAN_Prep/Code/prepmegan4cmaq_fert.x < prepmegan4cmaq.soil.$dom.inp | tee -a fert.$dom.log
mv ./Guangdong2023/grid_FERT.csv ./Guangdong2023/grid_FERT.$dom.csv

/bigdata/MEGAN/MEGANv3.2/MEGAN_Prep/Code/prepmegan4cmaq_nitrogen.x < prepmegan4cmaq.soil.$dom.inp | tee -a nitrogen.$dom.log
mv ./Guangdong2023/grid_NITROGEN.csv ./Guangdong2023/grid_NITROGEN.$dom.csv



end
