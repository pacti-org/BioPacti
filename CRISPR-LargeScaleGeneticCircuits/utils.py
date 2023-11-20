
from pacti.contracts import PolyhedralIoContract
from pacti.terms.polyhedra import PolyhedralTermList
from pacti.utils.lists import list_intersection, list_union, list_diff


# c1 receives inputs from c2. there is nothing to simplify in the guarantees of c1
def customcompose(c1: PolyhedralIoContract, c2:PolyhedralIoContract):
    sharedio = list_intersection(c1.inputvars, c2.outputvars)

    # compute assumptions
    terms_to_process : PolyhedralTermList = c1.a
    helper_terms : PolyhedralTermList = c2.g
    assumptions = []
    for assumption in terms_to_process.terms:
        for helper_term in helper_terms.terms:
            shared_vars = list_intersection(assumption.vars, helper_terms.vars)
            if len(shared_vars) == 1 and assumption.get_polarity(shared_vars[0]) == helper_term.get_polarity(shared_vars[0]):
                newterm = assumption.substitute_variable(var=shared_vars[0],subst_with_term=helper_term.isolate_variable(var_to_isolate=shared_vars[0]))
                assumptions.append(newterm)
    
    # 
    contract = PolyhedralIoContract(
        input_vars= list_diff(list_union(c1.inputvars,c2.inputvars),sharedio),
        output_vars= list_union(c1.outputvars,c2.outputvars),
        assumptions=c2.a,
        guarantees=PolyhedralTermList(assumptions)
        )


    print(contract)