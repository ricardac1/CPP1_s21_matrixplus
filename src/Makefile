GCC_FLAGS = -Wall -Wextra -Werror -std=c++17 -fprofile-arcs -ftest-coverage
SANITIZER = -g -fsanitize=address
EXE = test.out

# Папки поиска
SRC_DIRS = .
OBJ_DIR = ./objs
TEST_DIR = ./tests
TEST_OBJ_DIR = ./tests/objs
BUILD_PATH = gcov_report/
REPORT_PATH = $(BUILD_PATH)report/
TEST_CPP_FILES := $(wildcard $(TEST_DIR)/*.cpp)

# Поиск всех файлов .cpp во всех SRC_DIRS
SRC_CPP_FILES := $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))

# Создание списка объектных файлов, сохраняя структуру директорий
OBJS_CPP_FILES := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SRC_CPP_FILES))

# Создание списка объектных файлов для тестов
OBJS_TEST_FILES := $(patsubst %.cpp, $(TEST_OBJ_DIR)/%.o, $(notdir $(TEST_CPP_FILES)))

OS := $(shell uname)

ifeq ($(OS),Linux)
	OPEN = xdg-open
	TEST_FLAGS =  -L/usr/local/lib -lgtest -lgtest_main -lsubunit -lrt -lm -pthread
endif
ifeq ($(OS),Darwin)
	OPEN = open
	TEST_FLAGS =  -lgtest -lgtest_main -L/usr/local/lib
endif


all: s21_matrix_oop.a test

rebuild: clean build

build: clean
	g++ $(GCC_FLAGS) -o ${EXE} $(SRC_CPP_FILES) $(TEST_CPP_FILES) $(TEST_FLAGS) 
s21_matrix_oop.a: $(OBJS_CPP_FILES)
	ar rcs s21_matrix_oop.a $(OBJS_CPP_FILES)
	ranlib s21_matrix_oop.a

test: $(TEST_OBJ_DIR)/test_main.o $(OBJS_TEST_FILES) s21_matrix_oop.a
	g++ ./tests/objs/test_main.o s21_matrix_oop.a  $(TEST_FLAGS) -lgcov -o test.out
	./$(EXE)

# Компиляция исходных файлов в объектные
$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	mkdir -p $(dir $@)
	g++ $(GCC_FLAGS) -c $< -o $@

# Компиляция тестовых файлов в объектные
$(TEST_OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp | $(TEST_OBJ_DIR)
	mkdir -p $(dir $@)
	g++ $(GCC_FLAGS)  -c $< -o $@

# Создание директории для объектных файлов, если она не существует
$(OBJ_DIR) $(TEST_OBJ_DIR):
	mkdir -p $@

gcov_report: build_dir
	mkdir -p $(REPORT_PATH)
	g++ $(GCC_FLAGS) --coverage $(SRC_CPP_FILES) $(TEST_CPP_FILES) -o $(BUILD_PATH)$(EXE) $(TEST_FLAGS)  
	chmod +x $(BUILD_PATH)$(EXE)
	$(BUILD_PATH)$(EXE)
	lcov -t "Report" --directory . --capture --output-file $(BUILD_PATH)coverage.info --ignore-errors empty,inconsistent --no-external
	genhtml $(BUILD_PATH)coverage.info --output-directory $(REPORT_PATH)
	$(OPEN) $(REPORT_PATH)index.html

build_dir:
	mkdir -p $(BUILD_PATH)

clean:
	rm -rf $(OBJ_DIR) $(TEST_OBJ_DIR) s21_matrix_oop.a $(EXE) *.out *.gc* *.info $(BUILD_PATH)

docker_build:
	cp ../materials/linters/.clang-format ./
	docker build -t vault:v1 .
	docker create -it --name matrix_vault vault:v1

docker_start:
	docker start -i matrix_vault

docker_delete:
	docker rm matrix_vault
	docker image rm vault:v1

valgrind: rebuild
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=leaks_log.txt ./test.out
	echo --- Valgrind summary --- && cat leaks_log.txt | grep 'total heap usage' && cat leaks_log.txt | grep 'ERROR SUMMARY'

style_test:
	cp ../materials/linters/.clang-format .
	clang-format -n */*.cpp *.h *.cpp 
	rm .clang-format

format_style:
	cp ../materials/linters/.clang-format .
	clang-format -i */*.cpp *.h *.cpp 