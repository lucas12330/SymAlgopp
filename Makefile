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

.PHONY: all clean directories run_tests

# Cible par défaut
all: directories $(TEST_AST_BIN) $(TEST_DIFF_BIN)

# Création des dossiers temporaires et de sortie
directories:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

# Règle pour compiler les fichiers objets de la bibliothèque
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Règles pour compiler les tests en liant les objets de la bibliothèque
$(TEST_AST_BIN): $(TEST_AST_SRC) $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TEST_DIFF_BIN): $(TEST_DIFF_SRC) $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Lancer tous les tests
run_tests: all
	@echo "\n--- EXECUTION DU TEST AST ---"
	@./$(TEST_AST_BIN)
	@echo "\n--- EXECUTION DU TEST EQUATION DIFFERENTIELLE ---"
	@./$(TEST_DIFF_BIN)

# Nettoyage
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
