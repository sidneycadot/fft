
.PHONY : clean run

round_fixed_point_testbench e~round_fixed_point_testbench.o : work-obj93.cf round_fixed_point.o
	ghdl -e round_fixed_point_testbench

work-obj93.cf round_fixed_point.o : round_fixed_point.vhdl
	ghdl -a round_fixed_point.vhdl

run : round_fixed_point_testbench
	./round_fixed_point_testbench

clean :
	$(RM) round_fixed_point_testbench e~round_fixed_point_testbench.o work-obj93.cf round_fixed_point.o
