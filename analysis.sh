#/bin/bash

plot_energy_data="plot_energy.gnu"
plot_cluster_stat="plot_cluster_stat.gnu"
plot_rings_dist="plot_rings_dist.gnu"
plot_q_dist="plot_q_dist.gnu"
plot_nq_dist="plot_nq_dist.gnu"

plotter="/usr/bin/gnuplot"

$plotter $plot_energy_data
$plotter $plot_cluster_stat
#$plotter $plot_rings_dist
$plotter $plot_q_dist
$plotter $plot_nq_dist
