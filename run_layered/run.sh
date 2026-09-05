gfortran create_acquisition_oneline.f90
./a.out

mpirun -np 1 ../main mode=1 \
	fvm=0 \
	airbc=1 \
	preco=1 \
	algopt=1 \
	tol=1e-8 \
	niter=1000 \
	freqs=1.25 0.25,0.75,1.25 \
	chsrc=Ex \
	chrec=Ex,Ey,Ez \
	nb=10 \
	fx1=fx1 \
	fx2=fx2 \
	fx3=fx3 \
	frho11=frho11 \
	frho22=frho22 \
	frho33=frho33 \
	fsrc=sources.txt \
	frec=receivers.txt \
	fsrcrec=src_rec_table.txt \
	param_extended=1

