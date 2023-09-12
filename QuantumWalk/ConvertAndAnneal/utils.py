import numpy as np
import scipy as sp
import itertools
import functools as ft
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import re

from sympy.parsing.sympy_parser import parse_expr
from sympy import Add
import sympy as sp
import copy

#https://obliviateandsurrender.github.io/blogs/vqe.html

def HamiltonianToPauli(Hamiltonian):
    """Takes a Hamiltonian and returns Pauli strings (String) and coefficient"""
    logdim_H = int(np.log2(len(Hamiltonian)))

    # Sanity Checks: Power of 2 and Hermiticity
    if Hamiltonian.shape != (2**logdim_H, 2**logdim_H):
        raise ValueError(
            "The Hamiltonian should have shape (2**n, 2**n), for any qubit number n>=1"
        )

    if not np.allclose(Hamiltonian, Hamiltonian.conj().T):
        raise ValueError("The Hamiltonian is not Hermitian")

    #Pauli matrices
    I = np.eye(2, 2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0,-1]], dtype=complex)
    Y = complex(0,-1)*np.matmul(Z,X)
    Pauli_matrices = [I, X, Y, Z]
    Pauli_labels = ['I', 'X', 'Y', 'Z']

    combo_labels = []
    combo_coeffs = []

    #Iterate over all Pauli combinations and corresponding strings
    #coeff of a combo = 1/N *trace(H*combo)
    for combination_string, combination_matrices in zip(itertools.product(Pauli_labels,repeat=logdim_H),itertools.product(Pauli_matrices,repeat=logdim_H)):
        # Tensor product of matrices
        joint_matrix = ft.reduce(np.kron,combination_matrices)
        #Coeff
        coeff = (1/2**logdim_H)*np.trace(np.matmul(joint_matrix,Hamiltonian))

        #If coefficient is non-zero, store it
        if not np.allclose(coeff,0):
            combo_coeffs.append(coeff)
            combo_labels.append(''.join(combination_string))

    return combo_labels, combo_coeffs

def PauliToHamiltonian (combo_labels,combo_coeffs,dim):
    max_term_len = len(max(combo_labels, key=len))
    if not np.allclose(2 ** max_term_len, dim):
        raise ValueError("The given dimension ", dim, " mismatches with the dimension from the Pauli term length ",
                         2 ** max_term_len)


    #Pauli matrices
    I = np.eye(2, 2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0,-1]], dtype=complex)
    Y = complex(0,-1)*np.matmul(Z,X)
    Pauli_labels = ['I', 'X', 'Y', 'Z']
    Pauli_matrices = {lab: mat for lab, mat in zip(Pauli_labels,[I,X,Y,Z])}

    Hamiltonian = np.zeros((dim,dim),dtype=complex)
    for combo_label,combo_coeff in zip(combo_labels,combo_coeffs):

        #There should be max_term_len elements in each term. If termis shorter fill in identities to the left
        combo_label = 'I'*(max_term_len-len(combo_label))+combo_label

        #Calculate the corresponding tensor product
        matrices = [Pauli_matrices[m_lab] for m_lab in combo_label]
        Hamiltonian += combo_coeff*ft.reduce(np.kron,matrices)

    #print(Hamiltonian)
    return Hamiltonian


def PauliToString(labels,coeffs):
    """

    :param labels:
    :param coeffs:
    :return: String of sum of Pauli terms with position specified for each operator: X1*X2 etc.
    """
    string_sum = ""
    string_rawsum = ""

    for label,coeff in zip(labels,coeffs):
        if not string_rawsum:
            string_rawsum = str(coeff)+"*"+label
        else:
            string_rawsum += str(coeff)+"*"+label

        new_label = ''
        for position,op in enumerate(label):
            if op != 'I':
                if new_label:
                    new_label += '*'

                new_label += op+str(position+1)

        #If new_label empty it was all I : constant
        if not new_label:
            new_label = 'I'

        if string_sum:
            string_sum += '+'
        string_sum += str(coeff)+ '*' + new_label
    return string_sum,string_rawsum

def isPauliStringQuadratic(Pauli_string):
    """
    String of type (<num>)*Op1*Op2*..+( ...
    :param Pauli_string:
    :return:
    """
    Pauli_expr = parse_expr(Pauli_string)
    coeff_from_terms = Pauli_expr.as_coefficients_dict()

    for term in coeff_from_terms:
        if len(term.atoms(sp.Symbol)) > 2:
            return False
    return True

"""
    terms = Pauli_string.split('+(')
    max_term_len = 0
    for term in terms:
        term_len = 0
        parts = term.split('*')
        term_len = len(parts) -1        #leave first coefficient
        #print(term, term_len)
        max_term_len = max(term_len,max_term_len)
    if max_term_len > 2:
        return  False
    return True
"""

def isPauliStringQuadratizable (Pauli_string):
    #Has only Z terms
    Pauli_expr = parse_expr(Pauli_string)
    for sym in Pauli_expr.atoms(sp.Symbol):
        if str(sym)[0] != 'Z' and str(sym)[0] != 'z':
            print(sym, str(sym)[0], " non ising")
            return False
    return True
   # coeff_from_terms = Pauli_expr.as_coefficients_dict()

  #  for term in coeff_from_terms:
  #      if len(term.atoms(sp.Symbol)) > 2:
  #          return False
  #  return True
   # return True
def quadratizeString(string_input):
    return quadratizeExpr(parse_expr(string_input))

def quadratizeExpr(sym_expr):
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
     #           print( "Breaking ", term)
                sym_iter = iter(sorted(term.free_symbols,key=str))
                sym_var1 = next(sym_iter)
                sym_var2 =  next(sym_iter)
       #         print(" Choosing ", sym_var1,sym_var2)

                #new_sym_var = sp.symbols(str(sym_var1)+'a')
                #while new_sym_var in sym_expr.free_symbols:
                new_sym_var = sp.symbols(str(sym_var1)[0]+str(var_num))
                var_num += 1
               # print(" New symbol ", new_sym_var)

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
     #           print(" Deleting ", terms_to_delete)
      #          print(" Adding ", terms_to_add)
                break
        for term in terms_to_delete:
            coeff_from_terms.pop(term)

        for term in terms_to_add:
            coeff_from_terms[term] = terms_to_add[term]
   #     print("After first round ",coeff_from_terms)

    sym_expr = None
    for term in coeff_from_terms:
        if not sym_expr:
            sym_expr = coeff_from_terms[term]*term
        else:
            sym_expr += coeff_from_terms[term]*term


    return str(sym_expr)

def quadratizePauliString(PauliString):
    if isPauliStringQuadratizable(PauliString):
        return quadratizeString(PauliString)
    else:
        raise ValueError("The Pauli String cannot be quadratized as it is not Ising")


if __name__ == "__main__":
    Hamiltonian = np.array([[0, 6, 0, 0], [6, 2, 6., 0], [0, 6, 4, 6], [0, 0, 6, 6.0]])
    #print(Hamiltonian)
    #Hamiltonian = np.array([[1., 0, 0, 0], [0, 0, -1, 0], [0, -1, 0, 0], [0, 0, 0, 1]])
    labels,coeffs=HamiltonianToPauli(Hamiltonian)
    #print(labels, coeffs)
    #hm = PauliToHamiltonian(labels,coeffs,4)

    Pauli_string, _ = PauliToString(labels,coeffs)
    #Pauli_string = "Z1*z2*Z3"
    print(Pauli_string)
    #if isPauliStringQuadratic(Pauli_string):
    #    print("Pauli string is quadratic")
    #else:
    #    print(quadratizePauliString(Pauli_string))

    #string = "5*X1*X2 -7*X1*X2*X3*X4 +2*X1*X2*X3*X5"
    #print(quadratizeString(string))
