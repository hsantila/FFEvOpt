
slice=$(head -n2 electronDENSITY.xvg | awk '{dz=$1-prev;prev=$1}END{print dz}')

minbox=$(head -n 1 electronDENSITYsol.xvg | awk '{print $1}')
maxbox=$(tail -n 1 electronDENSITYsol.xvg | awk '{print $1}')
bulkDENS=$(awk -v minb=$minbox -v maxb=$maxbox '{if ($1<0.33+minb || $1>maxb-0.33)  {n=n+1; s=s+$2}} END{print s/n}' electronDENSITYsol.xvg)

cat electronDENSITY.xvg | awk -v slice=$slice -v bulkDENS=$bulkDENS 'BEGIN{scale=0.01;}{for(q=0;q<2000;q=q+1){Fa[q]=Fa[q]+($2-bulkDENS)*cos(scale*q*$1)*slice;Fb[q]=Fb[q]+($2-bulkDENS)*sin(scale*q*$1)*slice}}END{for(q=0;q<1000;q=q+1){print 0.1*q*scale" "0.01*sqrt(Fa[q]*Fa[q]+Fb[q]*Fb[q])}}'
mv electronDENSITY.xvg Electron_Density_From_Simulation.dat
