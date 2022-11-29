PROG = ggmcalc
FC = gfortran
FCFLAGS = -O3 -ffree-form -ffree-line-length-none -ffast-math -Ofast -fopenmp
LD = gfortran
FCMOD = nrtype.o coordinates_mod.o legendre_mod.o c_2n_mod.o gamma_mod.o w_mod.o u_mod.o undulation_mod.o height_anomaly_mod.o gravity_disturbance_mod.o gravity_anomaly_mod.o progressbar_mod.o date_sub.o duration.o

all: $(PROG)

$(PROG): $(FCMOD) $(PROG).o
	$(LD) $(FCFLAGS) $^ -o $@

.f.o:
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod $(PROG)
