# Compiler
CXX = clang++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra

# Include path for Catch2 (if installed globally, this might not be necessary)
# Uncomment and adjust the path if needed
# INCLUDES = -I/path/to/catch2

# Directories
SRC_DIR = src/cythonbiogeme/cpp
TEST_DIR = tests

# Name of the test executable
TEST_EXEC = test_runner

# Source files for the test
TEST_SRC = $(TEST_DIR)/test_bioExprBelongs_to.cc $(SRC_DIR)/bioExprBelongsTo.cc

all: $(TEST_EXEC)

$(TEST_EXEC): $(TEST_SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^

test: $(TEST_EXEC)
	./$(TEST_EXEC)

clean:
	rm -f $(TEST_EXEC)
