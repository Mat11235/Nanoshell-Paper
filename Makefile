CC = g++
CFLAGS = -Wall -I/usr/include/ -L/usr/local/lib
LIBS = -lgsl -lgslcblas -lm -larmadillo

all: bin/anl bin/nss bin/oGp bin/eV2nm bin/sgl bin/eV2ex bin/nsISS bin/sGp bin/sap bin/oap bin/rho2ome_sp bin/nsh bin/vrb bin/Esat bin/Gap bin/lyns bin/nom bin/sfr bin/fro

bin/anl: src/QS_OB_anl_time.cxx
	$(CC) $(CFLAGS) src/QS_OB_anl_time.cxx -o bin/anl $(LIBS)

bin/nss: src/nanoshell_ss_spe.cxx
	$(CC) $(CFLAGS) src/nanoshell_ss_spe.cxx -o bin/nss $(LIBS)

bin/oGp: src/nanoshell_omeG_p3.cxx
	$(CC) $(CFLAGS) src/nanoshell_omeG_p3.cxx -o bin/oGp $(LIBS)

bin/eV2nm: src/eV2nm.cxx
	$(CC) $(CFLAGS) src/eV2nm.cxx -o bin/eV2nm $(LIBS)

bin/sgl: src/single.cxx
	$(CC) $(CFLAGS) src/single.cxx -o bin/sgl $(LIBS)

bin/eV2ex: src/eV2ex.cxx
	$(CC) $(CFLAGS) src/eV2ex.cxx -o bin/eV2ex $(LIBS)

bin/nsISS: src/nanoshell_ISS.cxx
	$(CC) $(CFLAGS) src/nanoshell_ISS.cxx -o bin/nsISS $(LIBS)

bin/sGp: src/single_omeG_p3.cxx
	$(CC) $(CFLAGS) src/single_omeG_p3.cxx -o bin/sGp $(LIBS)

bin/sap: src/single_ome_al_p3.cxx
	$(CC) $(CFLAGS) src/single_ome_al_p3.cxx -o bin/sap $(LIBS)

bin/oap: src/nanoshell_ome_al_p3.cxx
	$(CC) $(CFLAGS) src/nanoshell_ome_al_p3.cxx -o bin/oap $(LIBS)

bin/rho2ome_sp: src/rho2ome_sp.cxx
	$(CC) $(CFLAGS) src/rho2ome_sp.cxx -o bin/rho2ome_sp $(LIBS)

bin/nsh: src/nanoshell_num.cxx
	$(CC) $(CFLAGS) src/nanoshell_num.cxx -o bin/nsh $(LIBS)

bin/vrb: src/variables.cxx
	$(CC) $(CFLAGS) src/variables.cxx -o bin/vrb $(LIBS)

bin/Esat: src/Esat.cxx
	$(CC) $(CFLAGS) src/Esat.cxx -o bin/Esat $(LIBS)

bin/Gap: src/nanoshell_G_p3.cxx
	$(CC) $(CFLAGS) src/nanoshell_G_p3.cxx -o bin/Gap $(LIBS)

bin/lyns: src/lycurguseV_ns.cxx
	$(CC) $(CFLAGS) src/lycurguseV_ns.cxx -o bin/lyns $(LIBS)

bin/nom: src/numOme_gsl.cxx
	$(CC) $(CFLAGS) src/numOme_gsl.cxx -o bin/nom $(LIBS)

bin/sfr: src/sfrohlich.cxx
	$(CC) $(CFLAGS) src/sfrohlich.cxx -o bin/sfr $(LIBS)

bin/fro: src/frohlich.cxx
	$(CC) $(CFLAGS) src/frohlich.cxx -o bin/fro $(LIBS)

clean:
	rm -f bin/anl bin/nss bin/oGp bin/eV2nm bin/sgl bin/eV2ex bin/nsISS bin/sGp bin/sap bin/oap bin/rho2ome_sp bin/nsh bin/vrb bin/Esat bin/Gap bin/lyns bin/nom bin/sfr bin/fro
