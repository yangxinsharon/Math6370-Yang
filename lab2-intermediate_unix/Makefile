
driver.exe:driver.o one_norm.o vector_difference.o vector_product.o vector_sum.o
	clang++ -o driver.exe *.o -lm
driver.o:driver.cpp
	clang++ -c driver.cpp
one_norm.o:one_norm.cpp
	clang++ -c one_norm.cpp
vector_difference.o:vector_difference.cpp
	clang++ -c vector_difference.cpp
vector_product.o:vector_product.cpp 
	clang++ -c vector_product.cpp
vector_sum.o:vector_sum.cpp
	clang++ -c vector_sum.cpp
clean:
	rm -f *.o