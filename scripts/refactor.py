import re

def process_hpp():
    with open('include/ASTNode.hpp', 'r') as f:
        content = f.read()

    # 1. Add virtual methods to ASTNode
    astnode_methods = """
    /*
     * Nom : integrer
     * Description : Calcule la primitive symbolique de l'expression.
     * Utilisation : ExprPtr p = noeud->integrer();
     */
    virtual ExprPtr integrer() const = 0;
    
    /*
     * Nom : limite
     * Description : Calcule la limite symbolique quand x tend vers a.
     * Utilisation : ExprPtr l = noeud->limite(a);
     */
    virtual ExprPtr limite(double a) const = 0;
"""
    content = re.sub(r'(virtual ExprPtr clone\(\) const = 0;)', r'\1\n' + astnode_methods, content)

    # 2. Add override methods to subclasses
    classes_to_update = ['Constante', 'Variable', 'Addition', 'Soustraction', 'Multiplication', 'Division', 'Puissance', 'Sinus', 'Cosinus']
    override_methods = """
    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;
"""
    for cls in classes_to_update:
        pattern = r'(class ' + cls + r' : public .*?\{.*?)(bool estEgal\(const ASTNode& autre\) const override;)'
        content = re.sub(pattern, r'\1\2\n' + override_methods, content, flags=re.DOTALL)

    # 3. Add Logarithme class
    logarithme_class = """
class Logarithme : public FonctionUnaire {
public:
    explicit Logarithme(ExprPtr arg);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;
};
"""
    content = content.replace('// --- Fonctions Helpers', logarithme_class + '\n// --- Fonctions Helpers')

    # 4. Add ast_ln helper
    ln_helper = """
/*
 * Nom : ast_ln
 * Description : Fonction pour créer le noeud logarithme d'une expression.
 * Utilisation : ExprPtr l = ast_ln(expr);
 */
ExprPtr ast_ln(ExprPtr arg);
"""
    content = content + '\n' + ln_helper

    with open('include/ASTNode.hpp', 'w') as f:
        f.write(content)

def process_cpp():
    with open('src/ASTNode.cpp', 'r') as f:
        content = f.read()

    # 1. Update Puissance::derivee
    puissance_derivee_new = """ExprPtr Puissance::derivee() const {
    if (m_droite->estConstante()) { // (u^n)' = n*u^{n-1}*u'
        double n = m_droite->getValeurConstante();
        if (n == 0) return cst(0.0);
        return cst(n) * ast_pow(m_gauche, cst(n - 1)) * m_gauche->derivee();
    }
    // Cas general : (u^v)' = u^v * (v' * ln(u) + v * u' / u)
    auto ln_u = ast_ln(m_gauche);
    auto terme1 = m_droite->derivee() * ln_u;
    auto terme2 = m_droite * (m_gauche->derivee() / m_gauche);
    return ast_pow(m_gauche, m_droite) * (terme1 + terme2);
}"""
    # Replace the old Puissance::derivee
    content = re.sub(r'ExprPtr Puissance::derivee\(\) const \{.*?\n\}', puissance_derivee_new, content, flags=re.DOTALL)

    # We will append the implementations of integrer, limite and Logarithme at the end of the file
    with open('src/ASTNode.cpp', 'w') as f:
        f.write(content)

process_hpp()
process_cpp()
print("Refactoring ASTNode done.")
