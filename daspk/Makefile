# Makefile for DASPK3

daspk: examples/DM/Makefile

examples/DM/Makefile: ddaspk30.tar.gz
	gunzip -c ddaspk30.tar.gz | tar -xf -
	patch examples/DM/Makefile omu_ddaspk30.patch
	cd examples/DM; make heat heatilu; cd ../..

clean:
	rm -rf solver preconds examples *~
