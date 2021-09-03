#ERR = -fsanitize=address -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wno-suggest-attribute=format
#DOP_ERR = -Wcast-qual
OPT = -O3 -ffast-math -march=native  -mfpmath=sse  
a.out: main.o jordan.o time.o
	mpicxx $(OPT) $(ERR) $(DOP_ERR) jordan.o main.o time.o
main.o: main.cpp jordan.h
	mpicxx $(OPT) $(ERR) $(DOP_ERR) -c main.cpp
jordan.o: jordan.cpp jordan.h
	mpicxx $(OPT) $(ERR) $(DOP_ERR) -c jordan.cpp
time.o: time.cpp jordan.h
	mpicxx $(OPT) $(ERR) $(DOP_ERR) -c time.cpp
