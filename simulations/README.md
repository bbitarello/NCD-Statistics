#!/bin/bash
#14.05.2015
#I cloned this from Cesare, in order to make sure I don't accidentally corrupt his simulation files.
#Last mofified: May 2019

#Original dir

#/mnt/sequencedb/PoPGen/cesare/bs_genomescan_ncv_test/
#and
#/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/

#add tajd to the files

echo "tajd_AFR" > tmpN
awk '$2==0{print $10}' neutral_n100.msms_3000bp.msstats >> tmpN
paste  neutral_n100.msms_3000bp.ncv+hka.out tmpN > neutral_n100.msms_3000bp.ncv+hka+tajd.out
echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs5_f0.5_n100.msms_3000bp.msstats >>  tmp5
paste Tbs5_f0.5_n100.msms_3000bp.ncv+hka.out tmp5 > Tbs5_f0.5_n100.msms_3000bp.ncv+hka+tajd.out 
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs5_f0.4_n100.msms_3000bp.msstats >>  tmp4
paste Tbs5_f0.4_n100.msms_3000bp.ncv+hka.out tmp4 > Tbs5_f0.4_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs5_f0.3_n100.msms_3000bp.msstats >>  tmp3
paste Tbs5_f0.3_n100.msms_3000bp.ncv+hka.out tmp3 > Tbs5_f0.3_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs5_f0.2_n100.msms_3000bp.msstats >>  tmp2
paste Tbs5_f0.2_n100.msms_3000bp.ncv+hka.out tmp2 > Tbs5_f0.2_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs5_f0.1_n100.msms_3000bp.msstats >>  tmp1
paste Tbs5_f0.1_n100.msms_3000bp.ncv+hka.out tmp1 > Tbs5_f0.1_n100.msms_3000bp.ncv+hka+tajd.out
###################################################
#repat commands for 6000bp

echo "tajd_AFR" > tmpN
awk '$2==0{print $10}' neutral_n100.msms_6000bp.msstats >> tmpN
paste  neutral_n100.msms_6000bp.ncv+hka.out tmpN > neutral_n100.msms_6000bp.ncv+hka+tajd.out
echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs5_f0.5_n100.msms_6000bp.msstats >>  tmp5
paste Tbs5_f0.5_n100.msms_6000bp.ncv+hka.out tmp5 > Tbs5_f0.5_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs5_f0.4_n100.msms_6000bp.msstats >>  tmp4
paste Tbs5_f0.4_n100.msms_6000bp.ncv+hka.out tmp4 > Tbs5_f0.4_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs5_f0.3_n100.msms_6000bp.msstats >>  tmp3
paste Tbs5_f0.3_n100.msms_6000bp.ncv+hka.out tmp3 > Tbs5_f0.3_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs5_f0.2_n100.msms_6000bp.msstats >>  tmp2
paste Tbs5_f0.2_n100.msms_6000bp.ncv+hka.out tmp2 > Tbs5_f0.2_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs5_f0.1_n100.msms_6000bp.msstats >>  tmp1
paste Tbs5_f0.1_n100.msms_6000bp.ncv+hka.out tmp1 > Tbs5_f0.1_n100.msms_6000bp.ncv+hka+tajd.out
###########################################

#repeat commands for 12000bp

echo "tajd_AFR" > tmpN
awk '$2==0{print $10}' neutral_n100.msms_12000bp.msstats >> tmpN
paste  neutral_n100.msms_12000bp.ncv+hka.out tmpN > neutral_n100.msms_12000bp.ncv+hka+tajd.out

echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs5_f0.5_n100.msms_12000bp.msstats >>  tmp5
paste Tbs5_f0.5_n100.msms_12000bp.ncv+hka.out tmp5 > Tbs5_f0.5_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs5_f0.4_n100.msms_12000bp.msstats >>  tmp4
paste Tbs5_f0.4_n100.msms_12000bp.ncv+hka.out tmp4 > Tbs5_f0.4_n100.msms_12000bp.ncv+hka+tajd.out
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs5_f0.3_n100.msms_12000bp.msstats >>  tmp3
paste Tbs5_f0.3_n100.msms_12000bp.ncv+hka.out tmp3 > Tbs5_f0.3_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs5_f0.2_n100.msms_12000bp.msstats >>  tmp2
paste Tbs5_f0.2_n100.msms_12000bp.ncv+hka.out tmp2 > Tbs5_f0.2_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs5_f0.1_n100.msms_12000bp.msstats >>  tmp1
paste Tbs5_f0.1_n100.msms_12000bp.ncv+hka.out tmp1 > Tbs5_f0.1_n100.msms_12000bp.ncv+hka+tajd.out
####################################################################
###############################################################################
#tbs=3
echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs3_f0.5_n100.msms_3000bp.msstats >>  tmp5
paste Tbs3_f0.5_n100.msms_3000bp.ncv+hka.out tmp5 > Tbs3_f0.5_n100.msms_3000bp.ncv+hka+tajd.out 
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs3_f0.4_n100.msms_3000bp.msstats >>  tmp4
paste Tbs3_f0.4_n100.msms_3000bp.ncv+hka.out tmp4 > Tbs3_f0.4_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs3_f0.3_n100.msms_3000bp.msstats >>  tmp3
paste Tbs3_f0.3_n100.msms_3000bp.ncv+hka.out tmp3 > Tbs3_f0.3_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs3_f0.2_n100.msms_3000bp.msstats >>  tmp2
paste Tbs3_f0.2_n100.msms_3000bp.ncv+hka.out tmp2 > Tbs3_f0.2_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs3_f0.1_n100.msms_3000bp.msstats >>  tmp1
paste Tbs3_f0.1_n100.msms_3000bp.ncv+hka.out tmp1 > Tbs3_f0.1_n100.msms_3000bp.ncv+hka+tajd.out
###################################################
#repat commands for 6000bp

echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs3_f0.5_n100.msms_6000bp.msstats >>  tmp5
paste Tbs3_f0.5_n100.msms_6000bp.ncv+hka.out tmp5 > Tbs3_f0.5_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs3_f0.4_n100.msms_6000bp.msstats >>  tmp4
paste Tbs3_f0.4_n100.msms_6000bp.ncv+hka.out tmp4 > Tbs3_f0.4_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs3_f0.3_n100.msms_6000bp.msstats >>  tmp3
paste Tbs3_f0.3_n100.msms_6000bp.ncv+hka.out tmp3 > Tbs3_f0.3_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs3_f0.2_n100.msms_6000bp.msstats >>  tmp2
paste Tbs3_f0.2_n100.msms_6000bp.ncv+hka.out tmp2 > Tbs3_f0.2_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs3_f0.1_n100.msms_6000bp.msstats >>  tmp1
paste Tbs3_f0.1_n100.msms_6000bp.ncv+hka.out tmp1 > Tbs3_f0.1_n100.msms_6000bp.ncv+hka+tajd.out
###########################################
#repat commands for 12000bp

echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs3_f0.5_n100.msms_12000bp.msstats >>  tmp5
paste Tbs3_f0.5_n100.msms_12000bp.ncv+hka.out tmp5 > Tbs3_f0.5_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs3_f0.4_n100.msms_12000bp.msstats >>  tmp4
paste Tbs3_f0.4_n100.msms_12000bp.ncv+hka.out tmp4 > Tbs3_f0.4_n100.msms_12000bp.ncv+hka+tajd.out
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs3_f0.3_n100.msms_12000bp.msstats >>  tmp3
paste Tbs3_f0.3_n100.msms_12000bp.ncv+hka.out tmp3 > Tbs3_f0.3_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs3_f0.2_n100.msms_12000bp.msstats >>  tmp2
paste Tbs3_f0.2_n100.msms_12000bp.ncv+hka.out tmp2 > Tbs3_f0.2_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs3_f0.1_n100.msms_12000bp.msstats >>  tmp1
paste Tbs3_f0.1_n100.msms_12000bp.ncv+hka.out tmp1 > Tbs3_f0.1_n100.msms_12000bp.ncv+hka+tajd.out
###############################
#TBs1


echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs1_f0.5_n100.msms_3000bp.msstats >>  tmp5
paste Tbs1_f0.5_n100.msms_3000bp.ncv+hka.out tmp5 > Tbs1_f0.5_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs1_f0.4_n100.msms_3000bp.msstats >>  tmp4
paste Tbs1_f0.4_n100.msms_3000bp.ncv+hka.out tmp4 > Tbs1_f0.4_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs1_f0.3_n100.msms_3000bp.msstats >>  tmp3
paste Tbs1_f0.3_n100.msms_3000bp.ncv+hka.out tmp3 > Tbs1_f0.3_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs1_f0.2_n100.msms_3000bp.msstats >>  tmp2
paste Tbs1_f0.2_n100.msms_3000bp.ncv+hka.out tmp2 > Tbs1_f0.2_n100.msms_3000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs1_f0.1_n100.msms_3000bp.msstats >>  tmp1
paste Tbs1_f0.1_n100.msms_3000bp.ncv+hka.out tmp1 > Tbs1_f0.1_n100.msms_3000bp.ncv+hka+tajd.out
###################################################
###############################
#repat commands for 6000bp

echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs1_f0.5_n100.msms_6000bp.msstats >>  tmp5
paste Tbs1_f0.5_n100.msms_6000bp.ncv+hka.out tmp5 > Tbs1_f0.5_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs1_f0.4_n100.msms_6000bp.msstats >>  tmp4
paste Tbs1_f0.4_n100.msms_6000bp.ncv+hka.out tmp4 > Tbs1_f0.4_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs1_f0.3_n100.msms_6000bp.msstats >>  tmp3
paste Tbs1_f0.3_n100.msms_6000bp.ncv+hka.out tmp3 > Tbs1_f0.3_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs1_f0.2_n100.msms_6000bp.msstats >>  tmp2
paste Tbs3_f0.2_n100.msms_6000bp.ncv+hka.out tmp2 > Tbs1_f0.2_n100.msms_6000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs1_f0.1_n100.msms_6000bp.msstats >>  tmp1
paste Tbs1_f0.1_n100.msms_6000bp.ncv+hka.out tmp1 > Tbs1_f0.1_n100.msms_6000bp.ncv+hka+tajd.out
###########################################
#repat commands for 12000bp
echo "tajd_AFR" > tmp5
awk '$2==0{print $10}' Tbs1_f0.5_n100.msms_12000bp.msstats >>  tmp5
paste Tbs1_f0.5_n100.msms_12000bp.ncv+hka.out tmp5 > Tbs1_f0.5_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp4
awk '$2==0{print $10}' Tbs1_f0.4_n100.msms_12000bp.msstats >>  tmp4
paste Tbs1_f0.4_n100.msms_12000bp.ncv+hka.out tmp4 > Tbs1_f0.4_n100.msms_12000bp.ncv+hka+tajd.out
echo "tajd_AFR" > tmp3
awk '$2==0{print $10}' Tbs1_f0.3_n100.msms_12000bp.msstats >>  tmp3
paste Tbs1_f0.3_n100.msms_12000bp.ncv+hka.out tmp3 > Tbs1_f0.3_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp2
awk '$2==0{print $10}' Tbs1_f0.2_n100.msms_12000bp.msstats >>  tmp2
paste Tbs1_f0.2_n100.msms_12000bp.ncv+hka.out tmp2 > Tbs1_f0.2_n100.msms_12000bp.ncv+hka+tajd.out
#
echo "tajd_AFR" > tmp1
awk '$2==0{print $10}' Tbs1_f0.1_n100.msms_12000bp.msstats >>  tmp1
paste Tbs1_f0.1_n100.msms_12000bp.ncv+hka.out tmp1 > Tbs1_f0.1_n100.msms_12000bp.ncv+hka+tajd.out
###############################



