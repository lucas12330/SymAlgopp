# ==========================================
# Makefile pour le projet SymAlgo++
# ==========================================

# Compilateur et options
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O2 -Iinclude -isystem vendor/eigen
# Génère un fichier .d par objet pour recompiler quand un header change
DEPFLAGS = -MMD -MP

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
TEST_SRCS = $(wildcard $(TEST_DIR)/test_*.cpp)
TEST_BINS = $(patsubst $(TEST_DIR)/%.cpp,$(BIN_DIR)/%,$(TEST_SRCS))

# Démonstration
DEMO_BIN = $(BIN_DIR)/demo

# Cibles de benchmark
BENCH_SRC = $(wildcard $(BENCH_DIR)/*.cpp)
BENCH_BIN = $(BIN_DIR)/bench_suite

.PHONY: all clean run_tests bench

# Cible par défaut
all: $(TEST_BINS) $(DEMO_BIN)

# Création des dossiers
$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

# Règle pour compiler les fichiers objets de la bibliothèque
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

# Règles pour compiler les tests en liant les objets de la bibliothèque
$(BIN_DIR)/test_%: $(TEST_DIR)/test_%.cpp $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $^ -o $@

$(DEMO_BIN): demo.cpp $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $^ -o $@

# Lancer tous les tests : s'arrête (code de retour non nul) au premier échec
run_tests: $(TEST_BINS)
	@for t in $(TEST_BINS); do \
		printf '\n--- EXECUTION DE %s ---\n' "$$t"; \
		./$$t || exit 1; \
	done

# Règle pour compiler les benchmarks (avec optimisation maximale)
$(BENCH_BIN): $(BENCH_SRC) $(SRCS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -O3 -DNDEBUG $^ -o $@ -lbenchmark -lpthread -lginac -lcln

bench: $(BENCH_BIN)
	@printf '\n--- EXECUTION DES BENCHMARKS ---\n'
	@./$(BENCH_BIN)

# Nettoyage
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

-include $(OBJS:.o=.d) $(TEST_BINS:$(BIN_DIR)/%=$(BIN_DIR)/%.d) $(DEMO_BIN).d
