OS := $(shell uname -s)

ifeq ($(OS), Darwin)
	FLAGS_PLATFORM = -lcheck
else
	FLAGS_PLATFORM = -lcheck -lsubunit -lrt -lm -lpthread -g
endif

all: clean build

build: s21_matrix.a

s21_matrix.a: s21_matrix.c
	@gcc -c -Wall -Werror -Wextra s21_matrix.c
	@ar rcs s21_matrix.a s21_matrix.o

test: s21_matrix.a
	@gcc tests.c s21_matrix.a $(FLAGS_PLATFORM) -o tests.o
	@./tests.o

gcov_report:
	@gcc --coverage s21_matrix.c tests.c $(FLAGS_PLATFORM) -o gсov_report.o
	@./gсov_report.o
	@lcov -t s21_matrix_tests -o s21_matrix_tests.info -c -d .
	@genhtml -o report s21_matrix_tests.info

check:
	@cp ../materials/linters/CPPLINT.cfg ./
	@python3 ./../materials/linters/cpplint.py --extensions=c s21_matrix.c tests.c
	@python3 ./../materials/linters/cpplint.py --extensions=h s21_matrix.h
	@rm -f CPPLINT.cfg

clean:
	@rm -rf *.a
	@rm -rf *.o
	@rm -rf report
	@rm -rf *.gcda
	@rm -rf *.gcdo
	@rm -rf *.gcno
	@rm -rf *.info
