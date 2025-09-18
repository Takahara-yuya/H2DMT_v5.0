gmt begin result_plot png
  
#gmt makecpt -T-0.5/3.5/0.02 -Cwysiwyg -I -H > colorobem2.cpt

gmt subplot begin 2x1 -Fs15c/5c -M1c

  gmt subplot set 0

  gmt basemap -R-20/1860/-2/100 -JX15c/-5c -Bxaf+l"Distance(km)" -Byaf+l"Depth(km)" -BWSen
  gmt plot ../Result/Inv_model.gmt -L -Ccolorobem.cpt

 #gmt plot ../Result/tmp/mesh/inv_nlcg_itr80.gmt -L -Ccolorobem.cpt
  
 gmt subplot set 1

  gmt basemap -R-20/1860/-2/800 -JX15c/-5c -Bxaf+l"Distance(km)" -Byaf+l"Depth(km)" -BWSen
  gmt plot ../Result/Inv_model.gmt -L -Ccolorobem.cpt

  gmt colorbar -Ccolorobem2.cpt  -B+l"Ohm-m"

gmt subplot end

gmt end show

