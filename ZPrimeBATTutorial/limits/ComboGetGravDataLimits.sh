./templateLimitBATcombSys 0 0 0 0 1 1 1 0 1 1 0 1> combo_data_output_0.txt
./templateLimitBATcombSys 0 0 0 1 1 1 1 0 1 1 0 1> combo_data_output_1.txt
./templateLimitBATcombSys 0 0 0 2 1 1 1 0 1 1 0 1> combo_data_output_2.txt
./templateLimitBATcombSys 0 0 0 3 1 1 1 0 1 1 0 1> combo_data_output_3.txt
./templateLimitBATcombSys 0 0 0 4 1 1 1 0 1 1 0 1> combo_data_output_4.txt
./templateLimitBATcombSys 0 0 0 5 1 1 1 0 1 1 0 1> combo_data_output_5.txt
./templateLimitBATcombSys 0 0 0 6 1 1 1 0 1 1 0 1> combo_data_output_6.txt
./templateLimitBATcombSys 0 0 0 7 1 1 1 0 1 1 0 1> combo_data_output_7.txt
./templateLimitBATcombSys 0 0 0 8 1 1 1 0 1 1 0 1> combo_data_output_8.txt
./templateLimitBATcombSys 0 0 0 9 1 1 1 0 1 1 0 1> combo_data_output_9.txt

grep "limit" combo_data_output*.txt | awk '{print $2,$3}' | sort > GravitonCombo_ObservedLimit.txt
