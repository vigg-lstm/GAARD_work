Rscript focal_gwas_PBS.r Avrankou_coluzzii_Delta > Avrankou_coluzzii_Delta_focal_gwas.log 2>&1

Rscript focal_gwas_PBS.r Baguida_gambiae_Delta > Baguida_gambiae_Delta_focal_gwas.log 2>&1 &
Rscript focal_gwas_PBS.r Baguida_gambiae_PM > Baguida_gambiae_PM_focal_gwas.log 2>&1

Rscript focal_gwas_PBS.r Korle-Bu_coluzzii_Delta > Korle-Bu_coluzzii_Delta_focal_gwas.log 2>&1 &
Rscript focal_gwas_PBS.r Korle-Bu_coluzzii_PM 2La Ace1 > Korle-Bu_coluzzii_PM_focal_gwas.log 2>&1

Rscript focal_gwas_PBS.r Madina_gambiae_Delta > Madina_gambiae_Delta_focal_gwas.log 2>&1 &
Rscript focal_gwas_PBS.r Madina_gambiae_PM Ace1 > Madina_gambiae_PM_focal_gwas.log 2>&1

Rscript focal_gwas_PBS.r Obuasi_gambiae_Delta > Obuasi_gambiae_Delta_focal_gwas.log 2>&1 &
Rscript focal_gwas_PBS.r Obuasi_gambiae_PM Ace1 > Obuasi_gambiae_PM_focal_gwas.log 2>&1
