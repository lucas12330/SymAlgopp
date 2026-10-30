/**
 * @file incertitudes.cpp
 * @brief Cas d'usage : propagation d'incertitudes d'une mesure de g avec un pendule.
 *
 * On mesure la longueur L du fil, la période T et l'amplitude th d'un pendule ; la valeur de g
 * s'en déduit par une formule non linéaire (avec la correction d'amplitude). Quelle incertitude
 * sur g, et quelle mesure faut-il améliorer en priorité ?
 *
 * La formule est lue depuis du texte, ses dérivées partielles sont calculées exactement, puis
 * évaluées par un seul programme compilé qui partage les sous-expressions communes. Une
 * simulation de Monte-Carlo vérifie le résultat. Le programme renvoie un code non nul si
 * l'écart dépasse la tolérance, ce qui en fait aussi un test.
 */

#include <symalgopp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace symalgo;

namespace {

// Une mesure : valeur moyenne et incertitude-type
struct Mesure {
    double valeur;
    double incertitude;
};

const std::string& nomDe(const ExprPtr& variable) { return comme<Variable>(variable)->getNom(); }

} // namespace

int main() {
    // g = 4 pi^2 L (1 + th^2/16)^2 / T^2, T = 2 pi sqrt(L/g) (1 + th^2/16)
    const ExprPtr g = lire("4*pi^2*L*(1 + th^2/16)^2/T^2", {"L", {"T", "th"}});

    // Chronométrage à la main : la période est la mesure la moins précise
    const std::map<std::string, Mesure> mesures = {
        {"L", {0.9950, 0.0020}},  // m
        {"T", {2.0070, 0.0500}},  // s
        {"th", {0.150, 0.010}},   // rad
    };

    const std::vector<ExprPtr> vars = variables(g); // L, T, th
    const std::size_t n = vars.size();
    std::vector<double> moyennes, sigmas;
    for (const ExprPtr& v : vars) {
        moyennes.push_back(mesures.at(nomDe(v)).valeur);
        sigmas.push_back(mesures.at(nomDe(v)).incertitude);
    }

    // Dérivées exactes : gradient (sensibilités) et diagonale de la hessienne (biais)
    const std::vector<ExprPtr> gradientG = gradient(g, vars);
    const auto H = hessienne(g, vars);

    // Un seul programme pour g, son gradient et les dérivées secondes pures
    std::vector<ExprPtr> sorties = {g};
    sorties.insert(sorties.end(), gradientG.begin(), gradientG.end());
    for (std::size_t i = 0; i < n; ++i) sorties.push_back(H[i][i]);
    const ProgrammeEvaluation programme(sorties, vars);
    const std::vector<double> v = programme.evaluerEn(moyennes);

    const double g0 = v[0];
    double variance = 0.0, biais = 0.0;
    std::vector<double> contribution(n);
    for (std::size_t i = 0; i < n; ++i) {
        contribution[i] = v[1 + i] * v[1 + i] * sigmas[i] * sigmas[i];
        variance += contribution[i];
        biais += 0.5 * v[1 + n + i] * sigmas[i] * sigmas[i];
    }
    const double sigma = std::sqrt(variance);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "g = " << g << "\n\n";
    std::cout << "Sensibilités exactes (dérivées partielles) :\n";
    for (std::size_t i = 0; i < n; ++i) {
        std::cout << "  dg/d" << nomDe(vars[i]) << " = " << gradientG[i] << "\n";
    }
    std::cout << "\nEn (L, T, th) = (" << moyennes[0] << ", " << moyennes[1] << ", " << moyennes[2] << ") :\n";
    std::cout << "  g = " << g0 << " m/s²,  incertitude-type (au premier ordre) = " << sigma << " m/s²\n";
    std::cout << "  part de chaque mesure dans la variance :\n";
    for (std::size_t i = 0; i < n; ++i) {
        std::cout << "    " << std::setw(2) << nomDe(vars[i]) << " : " << std::setw(6) << std::setprecision(1)
                  << 100.0 * contribution[i] / variance << " %\n";
    }
    std::cout << std::setprecision(4) << "  biais de la formule (hessienne) = " << biais << " m/s²\n";
    std::cout << "  programme compilé : " << programme.nombreInstructions() << " instructions pour "
              << programme.nombreSorties() << " sorties\n";

    // Vérification par Monte-Carlo : tirages gaussiens indépendants des trois mesures
    const ProgrammeEvaluation seulG({g}, vars);
    std::mt19937_64 alea(2026);
    std::normal_distribution<double> loi(0.0, 1.0);
    const int tirages = 1000000;
    double somme = 0.0, sommeCarres = 0.0;
    std::vector<double> point(n);
    for (int k = 0; k < tirages; ++k) {
        for (std::size_t i = 0; i < n; ++i) point[i] = moyennes[i] + sigmas[i] * loi(alea);
        double valeur;
        seulG.evaluerEn(point.data(), &valeur);
        somme += valeur;
        sommeCarres += valeur * valeur;
    }
    const double moyenneMC = somme / tirages;
    const double sigmaMC = std::sqrt(sommeCarres / tirages - moyenneMC * moyenneMC);
    const double erreurMoyenne = sigmaMC / std::sqrt(static_cast<double>(tirages));

    std::cout << "\nMonte-Carlo (" << tirages << " tirages) :\n";
    std::cout << "  moyenne - g = " << moyenneMC - g0 << " ± " << erreurMoyenne << "   (biais prévu " << biais << ")\n";
    std::cout << "  écart-type  = " << sigmaMC << "   (prévu " << sigma << ")\n";

    const bool ecartTypeOk = std::abs(sigmaMC / sigma - 1.0) < 0.02;
    const bool biaisOk = std::abs(moyenneMC - g0 - biais) < 6 * erreurMoyenne;
    std::cout << (ecartTypeOk && biaisOk ? "\nLa propagation analytique est confirmée.\n"
                                         : "\nÉCART : la propagation analytique n'est pas confirmée.\n");
    return ecartTypeOk && biaisOk ? 0 : 1;
}
