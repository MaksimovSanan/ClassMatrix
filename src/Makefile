all: clean test

.PHONY: test
test:
	g++ -std=c++17 -Wall -Werror -Wextra tests_class_matrix.cc class_matrix.cc -o test -lgtest -lgtest_main
	./test

style:
	cp ./../materials/linters/.clang-format ./
	clang-format -i *.cc *.h
	rm .clang-format

.PHONY: class_matrix.a
class_matrix.a : class_matrix.o
	ar rc libclass_matrix.a *.o
	ranlib libclass_matrix.a
	cp libclass_matrix.a class_matrix.a

.PHONY: class_matrix.o
class_matrix.o:
	g++  -std=c++17 -Wall -Werror -Wextra -c class_matrix.cc

clean:
	rm -rf *.o *.out *.gch *.dSYM *.gcov *.gcda *.gcno *.a tests_class_matrix *.css *.html vgcore* report *.info *.gz *.log test
