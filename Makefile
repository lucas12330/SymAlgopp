# ==========================================
# Makefile pour le projet SymAlgo++
# ==========================================

# Compilateur et options
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O2 -Iinclude -isystem vendor/eigen
# Génère un fichier .d par objet pour recompiler quand un header change
DEPFLAGS = -MMD -MP
# Arithmétique exacte de taille arbitraire (GMP)
LDLIBS = -lgmpxx -lgmp

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

# Compilation instrumentée (AddressSanitizer + UndefinedBehaviorSanitizer)
SAN_FLAGS = -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all -g -O1
SAN_BUILD_DIR = $(BUILD_DIR)/san
SAN_BIN_DIR = $(BIN_DIR)/san
SAN_OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(SAN_BUILD_DIR)/%.o,$(SRCS))
SAN_TEST_BINS = $(patsubst $(TEST_DIR)/%.cpp,$(SAN_BIN_DIR)/%,$(TEST_SRCS))

.PHONY: all clean run_tests check bench

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
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $^ -o $@ $(LDLIBS)

$(DEMO_BIN): demo.cpp $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $^ -o $@ $(LDLIBS)

# Lancer tous les tests : s'arrête (code de retour non nul) au premier échec
run_tests: $(TEST_BINS)
	@for t in $(TEST_BINS); do \
		printf '\n--- EXECUTION DE %s ---\n' "$$t"; \
		./$$t || exit 1; \
	done

# Tests sous sanitizers : toute fuite mémoire, lecture hors limites ou
# comportement indéfini fait échouer la cible
$(SAN_BUILD_DIR) $(SAN_BIN_DIR):
	mkdir -p $@

$(SAN_BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(SAN_BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(SAN_FLAGS) $(DEPFLAGS) -c $< -o $@

$(SAN_BIN_DIR)/test_%: $(TEST_DIR)/test_%.cpp $(SAN_OBJS) | $(SAN_BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SAN_FLAGS) $(DEPFLAGS) $^ -o $@ $(LDLIBS)

check: $(SAN_TEST_BINS)
	@for t in $(SAN_TEST_BINS); do \
		printf '\n--- EXECUTION SOUS SANITIZERS DE %s ---\n' "$$t"; \
		ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1 ./$$t || exit 1; \
	done

# Règle pour compiler les benchmarks (avec optimisation maximale)
$(BENCH_BIN): $(BENCH_SRC) $(SRCS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -O3 -DNDEBUG $^ -o $@ $(LDLIBS) -lbenchmark -lpthread -lginac -lcln

bench: $(BENCH_BIN)
	@printf '\n--- EXECUTION DES BENCHMARKS ---\n'
	@./$(BENCH_BIN)

# Nettoyage
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

-include $(OBJS:.o=.d) $(SAN_OBJS:.o=.d) $(TEST_BINS:%=%.d) $(SAN_TEST_BINS:%=%.d) $(DEMO_BIN).d
