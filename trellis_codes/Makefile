
# -------------------------------------------------------------------------------

.PHONY : clean

# -------------------------------------------------------------------------------

# -fno-stack-protector is required when compiling on platform with older 
# versions of glibc and gcc not supporting the stack protector

CFLAGS = -Wall #-fno-stack-protector

INCLUDES = -I/usr/local/include/

LFLAGS = -L/usr/local/lib/

LIBS = -lgsl -lgslcblas #-lqd

# -------------------------------------------------------------------------------
all : definitions.o basic_functions.o sync_functions.o decoders.o gen_rect_graph.o gen_rect_graph gen_rect_code.o gen_rect_code \
	perf_rect.o perf_rect opt_rect.o opt_rect concat.o concat sync_conv.o sync_conv

definitions.o : definitions.cc definitions.h
	g++ -c definitions.cc  $(CFLAGS) -o definitions.o

basic_functions.o : basic_functions.cc basic_functions.h definitions.h
	g++ -c basic_functions.cc $(CFLAGS) -o basic_functions.o

sync_functions.o : sync_functions.cc sync_functions.h basic_functions.h \
			definitions.h
	g++ -c sync_functions.cc $(CFLAGS) -o sync_functions.o

decoders.o : decoders.cc decoders.h sync_functions.h basic_functions.h \
		definitions.h
	g++ -c decoders.cc $(CFLAGS) -o decoders.o

gen_rect_graph.o : gen_rect_graph.cc definitions.o basic_functions.o
	g++ -c gen_rect_graph.cc $(CFLAGS) -o gen_rect_graph.o

gen_rect_graph : gen_rect_graph.o
	g++ basic_functions.o gen_rect_graph.o $(CFLAGS) $(LIBS) -o gen_rect_graph

gen_rect_code.o : gen_rect_code.cc definitions.o basic_functions.o
	g++ -c gen_rect_code.cc $(CFLAGS) -o gen_rect_code.o

gen_rect_code : gen_rect_code.o
	g++ basic_functions.o gen_rect_code.o $(CFLAGS) $(LIBS) -o gen_rect_code


all_functions = definitions.o basic_functions.o sync_functions.o decoders.o

# -------------------------------------------------------------------------------



# -------------------------------------------------------------------------------

perf_rect.o : perf_rect.cc $(all_functions)
	g++ -c perf_rect.cc $(CFLAGS) -o perf_rect.o

perf_rect : perf_rect.o
	g++ $(all_functions) perf_rect.o $(CFLAGS) $(LIBS) -o perf_rect

opt_rect.o : opt_rect.cc $(all_functions)
	g++ -c opt_rect.cc $(CFLAGS) -o opt_rect.o

opt_rect : opt_rect.o
	g++ $(all_functions) opt_rect.o $(CFLAGS) $(LIBS) -o opt_rect

# -------------------------------------------------------------------------------

concat.o : concat.cc $(all_functions)
	g++ -c concat.cc $(CFLAGS) -o concat.o

concat : $(all_functions) concat.o
	g++ $(all_functions) concat.o $(CFLAGS) $(LIBS) -o concat

clean:
# Le - permet d'éviter des erreurs si la commande échoue.
	rm -f *.o aaa concat perf_rect opt_rect sync_conv gen_rect_graph gen_rect_code

# -------------------------------------------------------------------------------
