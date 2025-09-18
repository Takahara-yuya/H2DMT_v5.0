gmt begin result_plot png
  
#gmt makecpt -T0.5/2.5/0.02 -Cwysiwyg -H -I > color1.cpt
  gmt basemap -R0/237/0/250 -JX15c/-5c -Bxaf+l"Distance(km)" -Byaf+l"Depth(km)" -BWSen
  gmt plot ../Result/Inv_model.gmt -L -Ccolor1.cpt
#gmt plot ../Result/tmp/mesh/inv_nlcg_itr17.gmt -L -Ccolor1.cpt


  gmt colorbar -Ccolor1.cpt  -B+l"Ohm-m"


gmt end show

