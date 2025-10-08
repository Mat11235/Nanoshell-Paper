CC = g++
CFLAGS = -Wall -I../src/headers -L../src/lib
LIBS = -lgsl -lgslcblas -lm -larmadillo

# List of test programs
TARGETS = time_behavior threshold steady_state tau_characterization oscillating_electricf electric_field_test

# Define how to build each target
all: $(TARGETS)

time_behavior: time_behavior.cxx
	$(CC) $(CFLAGS) time_behavior.cxx -o time_behavior $(LIBS)

threshold: threshold.cxx
	$(CC) $(CFLAGS) threshold.cxx -o threshold $(LIBS)

steady_state: steady_state.cxx
	$(CC) $(CFLAGS) steady_state.cxx -o steady_state $(LIBS)

tau_characterization: tau_characterization.cxx
	$(CC) $(CFLAGS) tau_characterization.cxx -o tau_characterization $(LIBS)

oscillating_electricf: oscillating_electricf.cxx
	$(CC) $(CFLAGS) oscillating_electricf.cxx -o oscillating_electricf $(LIBS)

electric_field_test: electric_field_test.cxx
	$(CC) $(CFLAGS) electric_field_test.cxx -o electric_field_test $(LIBS)

clean:
	rm -f $(TARGETS)
