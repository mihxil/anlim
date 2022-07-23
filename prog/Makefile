FC = g77 -O2 -g


ani:    ani.f correlation.f parameters.f
	$(FC) ani.f -o ani

testcorrelation: correlation.f testcorrelation.f
	$(FC) testcorrelation.f -o testcorrelation

clean:
	rm *~ ani testcorrelation