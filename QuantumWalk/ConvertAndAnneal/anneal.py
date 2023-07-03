import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

from dimod import BinaryQuadraticModel

import neal
import dimod

from misc import HamiltonianToPauli, PauliToString

def quadratizeExpr(sym_expr, expr_type):

    print('Quadratizing ',sym_expr)
    if expr_type == "SPIN":
        #convert to binary
        #spin = 2*binary - 1
        # spin -1 (down) : binary  0
        #spin 1 (up) : binary 1
        spins = sym_expr.free_symbols
        #print(spins)
        for spin in spins:
            #if str(spin)[0] == 'z' or str(spin)[0] == 'Z':
            bin_sym = sp.Symbol('b'+str(spin)[1:])
            sym_expr = sym_expr.subs(spin,2*bin_sym-1)
            sym_expr =sp.expand(sym_expr)
        print('Spin to Binary (z = 2b -1): ',sym_expr)

   # sym_expr = parse_expr(string_input)
    coeff_from_terms = sym_expr.as_coefficients_dict()
    #print(coeff_from_terms)

    done = True
    M = 1
    for term in coeff_from_terms:
        M += abs(coeff_from_terms[term])
        if len(term.atoms(sp.Symbol)) > 2:
            done = False

    var_num = len(sym_expr.free_symbols)+1
    while not done:
        done = True
        terms_to_delete = []
        terms_to_add = {}
        for term in coeff_from_terms:
            if len(term.atoms(sp.Symbol)) > 2:
                done = False
                #We will break this
    #            print( "Breaking ", term)
                sym_iter = iter(sorted(term.free_symbols,key=str))
                sym_var1 = next(sym_iter)
                sym_var2 =  next(sym_iter)
   #             print(" Choosing ", sym_var1,sym_var2)

                #new_sym_var = sp.symbols(str(sym_var1)+'a')
                #while new_sym_var in sym_expr.free_symbols:
                new_sym_var = sp.symbols(str(sym_var1)[0]+str(var_num))
                var_num += 1
                print(" New symbol ", new_sym_var, " = ",sym_var1, "*",sym_var2 )

                #wherever we find old symbols, replace them with the new symbol
                for term2 in coeff_from_terms:
                    #print("Searching ",term2)
                    if term2.has(sym_var1) and term2.has(sym_var2):
                        others = sp.Mul(*[t for t in term2.args if not (t.has(sym_var1) or t.has(sym_var2))])
                        if str(others) != "1":
                            #print(new_sym_var * others)
                            terms_to_add[new_sym_var * others]  =  coeff_from_terms[term2]
                            terms_to_delete.append(term2)

                if sym_var1 * sym_var2 not in coeff_from_terms:
                    terms_to_add[sym_var1 * sym_var2] = M
                else:
                    coeff_from_terms[sym_var1 * sym_var2] += M
                terms_to_add[sym_var1 * new_sym_var] = -2*M
                terms_to_add[sym_var2 * new_sym_var] = -2*M
                terms_to_add[new_sym_var] = 3*M
   #             print(" Deleting ", terms_to_delete)
    #            print(" Adding ", terms_to_add)
                break
        for term in terms_to_delete:
            coeff_from_terms.pop(term)

        for term in terms_to_add:
            coeff_from_terms[term] = terms_to_add[term]
    #    print("After first round ",coeff_from_terms)

    sym_expr = None
    for term in coeff_from_terms:
        if not sym_expr:
            sym_expr = coeff_from_terms[term]*term
        else:
            sym_expr += coeff_from_terms[term]*term


    return sym_expr

def constructAnnealBQM(quadratic_expr):
    #print(quadratic_expr)
    linear = {}
    quadratic = {}
    offset = 0
    coeff_from_terms = quadratic_expr.as_coefficients_dict()
    for term in coeff_from_terms:
        term_len = len(term.atoms(sp.Symbol))
        if term_len == 0:
            offset += coeff_from_terms[term]
        elif term_len == 1:
            linear[str(term)] =coeff_from_terms[term]
        else:
            sym_iter = iter(sorted(term.free_symbols, key=str))
            sym_var1 = next(sym_iter)
            sym_var2 = next(sym_iter)
            quadratic[(str(sym_var1),str(sym_var2))] = coeff_from_terms[term]
   # print(linear, quadratic, offset)
    bqm = BinaryQuadraticModel(linear, quadratic, offset, 'BINARY')
    return bqm

def zs(i,j,n):
    return i+ n*(j-1)

def getIexp(i,j,k,n):
    if j == k:
        return sp.Symbol('1')
    return parse_expr("(1/2)*(1+" + "z" + str(zs(i, j, n)) + "*" + "z" + str(zs(i, k, n)) + ")")

def getXexp(i,j,k,n):
    if j == k:
        return sp.Symbol('0')
    return parse_expr("(1/2)*(1-" + "z" + str(zs(i, j, n)) + "*" + "z" + str(zs(i, k, n)) + ")")

def getYexp(i,j,k,n):
    if j == k:
        return sp.Symbol('0')
    return sp.I*parse_expr("(1/2)*(" + "z" + str(zs(i, k, n)) + "-" + "z" + str(zs(i, j, n)) + ")")

def getZexp(i,j,k,n):
    if j == k:
        return parse_expr("z"+str(zs(i, j, n)))
    return parse_expr("(1/2)*(" + "z" + str(zs(i, j, n)) + "+" + "z" + str(zs(i, k, n)) + ")")

def generateNewHamiltonian(H_expr, num_orig_qubits, num_new_groups):
    Pauli_set = {'x', 'y', 'z', 'i'}
    new_Hamiltonian_exprs = {}
    for group1 in range(1,num_new_groups+1):
        for group2 in range(group1, num_new_groups+1):
            H_new_expr = sp.Symbol('0')

            print("Generating H (",group1,", ",group2, ")")
            H_by_term_expr = sp.Add.make_args(H_expr)
            #print(H_by_term_expr, len(H_by_term_expr))

           # if len(H_by_term_expr)  == 1:
            #    terms = H_by_term_expr
           # else:
            #    terms = H_by_term_expr[1:]
            terms =H_by_term_expr
            print(terms)
            for term in terms:
                print(" Processing term: ",term)
                #Break this into terms to get Pauli terms out
                term_by_sym_exp = sp.Mul.make_args(term)
                #account for initial I's
                expected_pos = 1
                new_ops = sp.Symbol('1')
                new_coeffs = sp.Symbol("1")
                print("  Processing symbols: ",term_by_sym_exp)
                for sym in term_by_sym_exp:
                    print("   Processing symbol: ", sym)
                    #Get coeffs out (coeffs = not a pauli term)
                    if str(sym)[0].lower() not in Pauli_set:
                        print("    Coeff ")
                        new_coeffs *= sym
                    else:
                        print("    Transforming Pauli symbol ")
                        #get_pos
                        if len(str(sym)) > 1:
                            pos_in_sym = int(str(sym)[1:])
                        elif len(str(sym)) ==1 and str(sym).lower() == 'i':
                            #Constant term
                            pos_in_sym = 1
                        else:
                            raise ValueError(
                                "Rogue Identity matrix in term"
                            )
                        #print(pos_in_sym)

                        # account for initial I's
                        while pos_in_sym > expected_pos:
                            #print(getIexp(expected_pos,group1,group2,num_orig_qubits), new_ops)
                            new_ops = new_ops * getIexp(expected_pos,group1,group2,num_orig_qubits)
                            expected_pos += 1
                            print("  Appended Is:", new_ops)
                        #transform the Pauli_sym
                        Pauli_sym  =str(sym)[0].lower()
                        if Pauli_sym == 'x':
                            new_ops *= getXexp(expected_pos,group1,group2,num_orig_qubits)
                        elif Pauli_sym == 'y':
                            new_ops *= getYexp(expected_pos,group1,group2,num_orig_qubits)
                        elif Pauli_sym == 'z':
                            print("    Choosing Z ")
                            new_ops *= getZexp(expected_pos, group1, group2, num_orig_qubits)
                        elif Pauli_sym == 'i':
                            new_ops *= getIexp(expected_pos, group1, group2, num_orig_qubits)
                        expected_pos += 1
                while expected_pos <= num_orig_qubits:
                    new_ops = new_ops * getIexp(expected_pos, group1, group2, num_orig_qubits)
                    expected_pos += 1
                    print("  Appended Is:", new_ops)

                print( "New Term: ", new_ops*new_coeffs)
                H_new_expr += new_coeffs*new_ops
            new_Hamiltonian_exprs[(group1,group2)] = parse_expr(str(H_new_expr))
            print("Hamiltonian (",group1,", ",group2," ) ",new_Hamiltonian_exprs[(group1,group2)])
    #H_new_expr = parse_expr(str(H_new_expr))

    return new_Hamiltonian_exprs




if __name__ == "__main__":
    Hamiltonian = np.array([[0, 6, 0, 0], [6, 2, 6., 0], [0, 6, 4, 6], [0, 0, 6, 6.0]])
    labels,coeffs=HamiltonianToPauli(Hamiltonian)
    Pauli_string, _ = PauliToString(labels,coeffs)
    print(Pauli_string)
    num_qubits = 2
    H_orig = "(-J/2)*(1+g)*x1*x2  + (-J/2)*(1-g)*y1*y2 + (-B)*z1 + (-B)*z2"
    H_orig = Pauli_string
    H_orig_exp = parse_expr(H_orig)
    repeat = 1
    num_new_groups = 4#2**num_qubits
    num_new_qubits = num_new_groups*num_qubits

    Hp_new_sub_exprs = generateNewHamiltonian(H_orig_exp, num_qubits, num_new_groups)
    print(Hp_new_sub_exprs)

    Hp_new_expr = sp.Symbol("0")

    for i in range(1,num_new_groups+1):
        for j in range(1,num_new_groups+1):
            if i < j:
                Hp_new_sub_expr = Hp_new_sub_exprs[(i, j)]
            else:
                Hp_new_sub_expr = Hp_new_sub_exprs[(j, i)]

            if i != j:
                Hp_new_expr += Hp_new_sub_expr * sp.Symbol('S' + str(i)) * sp.Symbol('S' + str(j))
            else:
                Hp_new_expr += Hp_new_sub_expr

    Hp_new_expr = parse_expr(str(Hp_new_expr))
    print(Hp_new_expr)

    Cp_exp = sp.Symbol("0")
    for pmcombo_num in range(2**num_qubits):
        pmcombo_bin = f'{pmcombo_num:0{num_qubits}b}'
        combo_term = sp.Symbol("0")
        for i in range(1,num_new_groups+1):
            group_term = sp.Symbol("1")
            for j in range(1,num_qubits+1):
                if pmcombo_bin[j-1] == '0':
                    group_term *= (1 + sp.Symbol("z" + str(zs(j, i, num_qubits)))) / 2
                else:
                    group_term *= (1 - sp.Symbol("z" + str(zs(j, i, num_qubits)))) / 2
            group_term *= sp.Symbol('S' + str(i))
            combo_term += group_term
        Cp_exp += combo_term ** 2

    Cp_exp= parse_expr(str(Cp_exp))
    print(Cp_exp)

    Hp_new_num_exp = Hp_new_expr.subs([(sp.Symbol('B'), 0.001), (sp.Symbol('J'), -0.01), (sp.Symbol('g'), 0)])
    Cp_num_exp = Cp_exp

    # Choose signs
    Sign_vals = [1,-1,1,1] #should be length of num_new_groups
    if len(Sign_vals) != num_new_groups:
        raise ValueError(
            "Sign values not given for all groups"
        )
    for i in range(1,num_new_groups+1):
        #Sign_vals.append(1)
        Cp_num_exp = Cp_num_exp.subs('S' + str(i), Sign_vals[i-1])
        Hp_new_num_exp = Hp_new_num_exp.subs('S' + str(i), Sign_vals[i-1])

    #Choose lambda, L

    L = sp.Symbol('100')
    D_pl = Hp_new_num_exp - L*Cp_num_exp

    D_pl = D_pl.subs(L,float(str(L)))
    D_pl_expanded = sp.expand(D_pl)
    print('Dpl ', D_pl_expanded)

    D_pl_expanded_sym = D_pl_expanded
    for sym in D_pl_expanded.free_symbols:
        D_pl_expanded_sym = D_pl_expanded_sym.subs(sym ** 2, 1)

    orig_vars = D_pl_expanded_sym.free_symbols
    num_orig_vars = len(orig_vars)
    D_pl_quad = quadratizeExpr(D_pl_expanded_sym,'SPIN')

    print('Quadratized ',D_pl_quad)

    bqm = constructAnnealBQM(D_pl_quad)

  #D-Wave
    #Define the sampler
    sampler = EmbeddingComposite(DWaveSampler())
    # sampler = neal.SimulatedAnnealingSampler()
    #Solve the BQM
    sampleset = sampler.sample(bqm, num_reads=100)
    #Get minimum energy record
    solution = sampleset.first
    print("solution.energy, solution.sample: ", solution.energy, solution.sample, end="\n")
    print("#####################################################################", end="\n")

   #Get the data
    bin_vals_from_var = solution.sample
    spin_vals_from_var = {}
    qubit_vals = []
    for bin_var in bin_vals_from_var:
        spin_vals_from_var['z'+bin_var[1:]] = 2*bin_vals_from_var[bin_var] -1
        if bin_vals_from_var[bin_var]  == 1:
            qubit_vals.append('0')
        else:
            qubit_vals.append('1')
    print(spin_vals_from_var)
    print("New Eigenvector: ", "".join([str(q) for q in qubit_vals[:num_orig_vars]]))
    print("New Cost ", solution.energy)

    Cp_num = Cp_num_exp
    for spin_var in orig_vars:
        Cp_num = Cp_num.subs(spin_var, spin_vals_from_var[str(spin_var)])
        #print(spin_var,spin_vals_from_var[str(spin_var)])
    print("Sum of bi^2: ", Cp_num) #This is sum of bi^2

    actual_eigenvalue = float(str(L)) + solution.energy / Cp_num
    print("Estimated eigenvalue", actual_eigenvalue)

    l, v = np.linalg.eig(Hamiltonian)
    print("Eigen Values:", l)
    print("Eigen Vectors", v)

# coeff_denom =np.sqrt(Cp_num)
#  coeff_num = 1
# print("Actual eigenvector (approx): ")


"""

    H_11 = "-B*z1 -B*z2"
    H_12 = "(1/4)*((-J/2)*(1+g)*(1-z1*z3)*(1-z2*z4) + (-J/2)*(1-g)*(z1-z3)*(z4-z2) + (-B)*(z1+z3)*(1+z2*z4) + (-B)*(z2+z4)*(1+z1*z3))"
    H_21 = H_12
    H_22 = "-B*z3  -B*z4"

    S1 = '1' #1 or -1
    S2 = '1' #1 or -1

    L = sp.Symbol('100')
    H_p = H_11 + '+' + '(' + S1 + '*' + S2 + "*" + H_12 + ')' + '+' + '(' + S2 + '*' + S1 + "*" + H_21 + ')' + '+' + H_22
    H_p_exp = parse_expr(H_p)
    print( 'H ', H_p_exp)

    H_p_num_exp = H_p_exp.subs([(sp.Symbol('B'),0.001),(sp.Symbol('J'),-0.01),(sp.Symbol('g'),0)])
    print('Hnum ', H_p_num_exp)

    Cp ="(1/16)*(((1+z1)*(1+z2)*S1 + (1+z3)*(1+z4)*S2)**2 + ((1+z1)*(1-z2)*S1 + (1+z3)*(1-z4)*S2)**2 + ((1-z1)*(1+z2)*S1 + (1-z3)*(1+z4)*S2)**2 + ((1-z1)*(1-z2)*S1 + (1-z3)*(1-z4)*S2)**2)"

    Cp_exp = parse_expr(Cp)
    #print(Cp_exp)

    D_pl = H_p_num_exp - L*Cp_exp
    D_pl = D_pl.subs([('S1',int(S1)),('S2',int(S2)),(L,int(str(L)))])
    #print(D_pl)
    D_pl_expanded = sp.expand(D_pl)
    print('Dpl ', D_pl_expanded)

    #Replace squares of terms
    D_pl_expanded_sym = D_pl_expanded
    for sym in D_pl_expanded.free_symbols:
        D_pl_expanded_sym = D_pl_expanded_sym.subs(sym ** 2, 1)

    orig_vars = D_pl_expanded_sym.free_symbols
    num_orig_vars = len(orig_vars)
    #print('Dpl after squares ', D_pl_expanded_sym)
    D_pl_quad = quadratizeExpr(D_pl_expanded_sym,'SPIN')

    print('Quadratized ',D_pl_quad)
    #print(D_pl_quad)
 #   simple_expr = parse_expr("128*z1*z2*z3 - 56*z1*z2 -48*z1*z3 +16*z2*z3 -52*z1 -52*z2 -96*z3 +196")
 #   simple_expr = parse_expr("z1 + z2 + z1*z2*z3")
  #  simple_quad_expr = quadratizeExpr(simple_expr,"BINARY")
  #  print(simple_expr)
  #  print(simple_quad_expr)
  #  bqm = constructAnnealBQM(simple_quad_expr)

    bqm = constructAnnealBQM(D_pl_quad)

    #D-Wave
    #Define the sampler
    ######sampler = EmbeddingComposite(DWaveSampler())
    sampler = neal.SimulatedAnnealingSampler()
    #Solve the BQM
    sampleset = sampler.sample(bqm, num_reads=100)
    #Get minimum energy record
    solution = sampleset.first
    print(solution.sample)
    #Get the data
    bin_vals_from_var = solution.sample
    spin_vals_from_var = {}
    qubit_vals = []
    for bin_var in bin_vals_from_var:
        spin_vals_from_var['z'+bin_var[1:]] = 2*bin_vals_from_var[bin_var] -1
        if bin_vals_from_var[bin_var]  == 1:
            qubit_vals.append('0')
        else:
            qubit_vals.append('1')
    print(spin_vals_from_var)
    print("New Eigenvector: ", "".join([str(q) for q in qubit_vals[:num_orig_vars]]))
    print("New Cost ", solution.energy)

    Cp_num =Cp_exp.subs([('S1',int(S1)),('S2',int(S2))])

    for spin_var in orig_vars:
        Cp_num = Cp_num.subs(spin_var,spin_vals_from_var[str(spin_var)])
       # print(spin_var,spin_vals_from_var[str(spin_var)])
    #print(Cp_num) #This is sum of bi^2

    actual_eigenvalue = int(str(L)) + solution.energy/Cp_num
    print("Actual eigenvalue", actual_eigenvalue)

   # coeff_denom =np.sqrt(Cp_num)
  #  coeff_num = 1
   # print("Actual eigenvector (approx): ")


    #print(sampleset.info.keys())
    #samples = sampleset.samples()
    #simple_expr_eval = simple_expr.subs([('z1',1),('z2',0),('z3',1),('z4' ,0)])
    #simple_quad_expr_eval = simple_quad_expr.subs([('z1',1),('z2',0),('z3',1),('z4' ,0)])
    #print(simple_quad_expr_eval,simple_expr_eval)

  #  H_eval = H_p_num_exp.subs([('z1',-1),('z2',-1),('z3',-1),('z4' ,-1)])
  #  Cp_eval = Cp_exp.subs([('z1',-1),('z2',-1),('z3',-1),('z4' ,-1)])
    #H_eval = H_p_num_exp.subs([('z1',1),('z2',1),('z3',1),('z4' ,1)])
 #   D_pl_expanded_sym_eval = D_pl_expanded_sym.subs([('z1',1),('z2',1),('z3',1),('z4' ,1)])
 #   print(D_pl_expanded_sym_eval)

 #   D_pl_quad_eval = D_pl_quad.subs([('b1',1),('b2',1),('b3',1),('b4' ,1),('b5',1),('b6',1),('b7',1),('b8' ,1)])
 #   D_pl_quad_eval = D_pl_quad.subs([('z1',1),('z2',1),('z3',1),('z4' ,1),('z5',1),('z6',1),('z7',1),('z8' ,1)])
    #D_pl_quad_eval = D_pl_quad.subs([('z1',-1),('z2',-1),('z3',-1),('z4' ,-1),('z5',-1),('z6',-1),('z7',-1),('z8' ,-1)])
    #D_pl_eval = D_pl.subs([('z1',-1),('z2',-1),('z3',-1),('z4' ,-1)])
  #  print(D_pl_quad_eval)

"""