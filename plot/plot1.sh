gmt begin result_plot png
  
#gmt makecpt -T1/3/0.02 -Cwysiwyg -H > color.cpt
  gmt basemap -R-0/10/0/10 -JX15c/-15c -Bxaf+l"Distance(km)" -Byaf+l"Depth(km)" -BWSen
  gmt plot ../Result/Inv_model.gmt -L -Ccolor.cpt

  gmt colorbar -Ccolor.cpt  -B+l"Ohm-m"


gmt end show

