# Makefile for DASPK3

daspk:
	tar -xzf ddaspk30.tar.gz
	patch examples/DM/Makefile omu_ddaspk30.patch
	cd examples/DM; make heat heatilu; cd ../..

clean:
	rm -rf solver preconds examples *~
