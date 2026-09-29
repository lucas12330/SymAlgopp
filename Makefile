# ==========================================
# Makefile pour le projet SymAlgo++
# ==========================================

# Compilateur et options
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -Iinclude -Ivendor/eigen

# Répertoires
SRC_DIR = src
INCLUDE_DIR = include
TEST_DIR = tests
BENCH_DIR = benchmarks
BUILD_DIR = build
BIN_DIR = bin

# Fichiers sources et objets de la bibliothèque
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Cibles de tests
TEST_AST_SRC = $(TEST_DIR)/test_ast.cpp
TEST_DIFF_SRC = $(TEST_DIR)/test_differentielle.cpp

TEST_AST_BIN = $(BIN_DIR)/test_ast
TEST_DIFF_BIN = $(BIN_DIR)/test_differentielle

# Cibles de benchmark
BENCH_SRC = $(wildcard $(BENCH_DIR)/*.cpp)
BENCH_BIN = $(BIN_DIR)/benchmark_suite

.PHONY: all clean run_tests bench

# Cible par défaut
all: $(TEST_AST_BIN) $(TEST_DIFF_BIN)

# Création des dossiers
$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

# Règle pour compiler les fichiers objets de la bibliothèque
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Règles pour compiler les tests en liant les objets de la bibliothèque
$(TEST_AST_BIN): $(TEST_AST_SRC) $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TEST_DIFF_BIN): $(TEST_DIFF_SRC) $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Lancer tous les tests
run_tests: all
	@echo "\n--- EXECUTION DU TEST AST ---"
	@./$(TEST_AST_BIN)
	@echo "\n--- EXECUTION DU TEST EQUATION DIFFERENTIELLE ---"
	@./$(TEST_DIFF_BIN)

# Règle pour compiler les benchmarks (avec optimisation maximale)
$(BENCH_BIN): $(BENCH_SRC) $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -O3 -DNDEBUG $^ -o $@ -lbenchmark -lpthread -lginac -lcln

bench: $(BENCH_BIN)
	@echo "\n--- EXECUTION DES BENCHMARKS ---"
	@./$(BENCH_BIN)

# Nettoyage
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
