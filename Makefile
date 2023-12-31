CFLAGS := -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format


all:
	g++ -pthread $(CFLAGS) main.cpp basic_funcs.cpp matrinit.cpp msr_funcs.cpp my_thread.cpp operations.cpp residuals.cpp


clear:
	rm ./a.out
