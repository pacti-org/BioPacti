#!/usr/bin/env python
# coding: utf-8

# # A toggle switch in mammalian cells
# 
# We consider a circuit of the form
# 
# ![image-2.png](attachment:image-2.png)
# 
# The dynamics are given by
# 
# $\dot{m}_i = (\alpha_{0,i} + \alpha \frac{a_i A_i}{d_i}) \frac{D_{i\_tot}}{1 + \frac{a_i A_i}{d_i} + \frac{a_i R_i}{d_i}} − \delta m_i$
# 
# $\dot{R}_i = \kappa m_{2-i} − \gamma R_i$
# 
# $\dot{A}_i = \kappa m_i − \gamma A_i$
# 
# 
# ## Specifications
# 
# We want to analyze the region on the parameter space where the system shows two stable states.
# - Parameters of interest: $\alpha_{0,i}$, $\alpha$, $D_{i,tot}$, $\gamma$
# - Constants: $a_i$, $d_i$, $\delta$
# - Requirement: ($m_1 > m_1^H$ AND $m_2 < m_2^L$) OR ($m_1 < m_1^L$ AND $m_2 > m_2^H$)
# 
# We will assume we are given the additional requirement
# 
# - ($A_1 > A_1^H$ AND $A_2 < A_2^L$) OR ($A_1 < A_1^L$ AND $A_2 > A_2^H$)
# 
# ### Analysis
# 
# In order to guarantee this requirement, we find the steady-state behavior of the system:
# 
# $
# \begin{cases}
# (\alpha_{0,i} + \alpha \frac{a_i A_i}{d_i}) \frac{D_{i\_tot}}{1 + \frac{a_i A_1}{d_i} + \frac{a_i A_2}{d_i}} = \delta m_i \\
# \kappa m_i = \gamma A_i
# \end{cases}
# $
# 
# We will introduce new variables $p_i = \alpha_{0,i} D_{i\_tot}$, $q_i = \alpha D_{i\_tot}$, all of which can be considered new parameters. We also introduce the new constants $b_i = a_i/d_i$. Our system becomes
# 
# $
# \frac{p_i + b_i q_i A_i }{1 + b_i (A_i + A_j)} = \delta m_i,
# $
# 
# where $i \ne j$.
# 
# We focus on the requirement $m_i > m_i^H$ and $m_j < m_j^H$. We can achieve this top-level objective by bounding $m_i$ below and $m_j$ above:
# 
# 1. $\delta m_i = F_i(A_i, A_j) = \frac{p_i + b_i q_i A_i }{1 + b_i (A_i + A_j)} \ge \frac{p_i + b_i q_i A_i^H }{1 + b_i (A_i^H + A_j^L)}$
# 2. $\delta m_j = F_j(A_i, A_j) = \frac{p_j + b_j q_j A_j }{1 + b_j (A_i + A_j)} \le \frac{p_j + b_j q_j A_j^L }{1 + b_j (A_i^H + A_j^L)}$
# 
# To satisfy the LHS of requirements 1 and 2, we need $F_i$ to increase in its first argument and decrease in its second, and we need $F_j$ to decrease in its first argument and increase in its second. The decreasing condition is immediate. To analyze the increasing condition, we make nonnegative the derivative wrt the corresponding arguments. We obtain the conditions:
# 
# - $q_i \ge p_i$
# - $q_j ( 1 + b_j A_i^H) \ge p_j$

# In[7]:


import numpy as np
from pacti.contracts import PolyhedralIoContract

###################################
## REQUIREMENTS
# requirement values first scenario

# Run this code for *randomly chosen* floating point 
# values for A_i_H between 5 to 50.
choices = np.random.uniform(0, 1, 10)
for A_j_L in choices:
    print("A_i_H = ", A_i_H)
    # A_i_H = 10
    A_j_L = 0.5
    m_i_H = 15
    m_j_L = 0.1
    # requirement values second scenario
    A_j_H = 10
    A_i_L = 0.5
    m_j_H = 10
    m_i_L = 0.1

    ###################################
    ## CONSTANTS
    # Where is alpha_i?
    delta = 0.1
    a_i = 1
    d_i = 10
    a_j = 1
    d_j = 9
    kappa = 1

    # computed constants
    b_i = a_i / d_i
    b_j = a_j / d_j


    top_level_sc1 = PolyhedralIoContract.from_strings(
        input_vars=[],
        output_vars=["m_i", "m_j", "A_i", "A_j"],
        assumptions=[],
        guarantees=[f"m_i >= {m_i_H}", f"A_i >= {A_i_H}", f"m_j <= {m_j_L}", f"A_j <= {A_j_L}"]
        )


    # We now code the contracts
    # 
    # 1. $\left(q_i \ge p_i, \delta m_i \ge \frac{p_i + b_i q_i A_i^H }{1 + b_i (A_i^H + A_j^L)}\right)$
    # 2. $\left(q_j ( 1 + b_j A_i^H) \ge p_j, \delta m_j \le \frac{p_j + b_j q_j A_j^L }{1 + b_j (A_i^H + A_j^L)}\right)$

    # In[8]:


    monotonicity_component1_sc1 = PolyhedralIoContract.from_strings(
        input_vars=["p_i","p_j","q_i","q_j"],
        output_vars=["m_i"],
        assumptions=["q_i >= p_i"],
        guarantees=[f"{delta}*m_i >= {1/(1 + b_i*(A_i_H + A_j_L))}*p_i + {A_i_H*b_i/(1 + b_i*(A_i_H + A_j_L))}*q_i"]
    )

    monotonicity_component2_sc1 = PolyhedralIoContract.from_strings(
        input_vars=["p_i","p_j","q_i","q_j"],
        output_vars=["m_j"],
        assumptions=[f"{1 + b_j*A_i_H} * q_j >= p_j"],
        guarantees=[f"{delta}*m_j <= {1/(1 + b_j*(A_i_H + A_j_L))}*p_j + {A_j_L*b_j/(1 + b_j*(A_i_H + A_j_L))}*q_j"]
    )

    A_definition_sc1 = PolyhedralIoContract.from_strings(
        input_vars=["gamma_inv"],
        output_vars=["A_i", "A_j"],
        assumptions=[],
        guarantees=[f"{kappa * m_i_H} * gamma_inv <= A_i ",
                    f"A_j <= {kappa * m_j_L}  * gamma_inv"]
    )

    # compute the parameters that yield the first set of top-level requirements
    test_contract = top_level_sc1.quotient(A_definition_sc1)
    parameters_contract_1 = top_level_sc1.quotient(monotonicity_component1_sc1.compose(monotonicity_component2_sc1).compose(A_definition_sc1))


    ##################
    ## Now we analyze the second scenario: when m_i < m_i_L and m_j > m_J_H

    top_level_sc2 = PolyhedralIoContract.from_strings(
        input_vars=[],
        output_vars=["m_i", "m_j", "A_i", "A_j"],
        assumptions=[],
        guarantees=[f"m_j >= {m_j_H}", f"A_j >= {A_j_H}", f"m_i <= {m_i_L}", f"A_i <= {A_i_L}"]
        )

    monotonicity_component1_sc2 = PolyhedralIoContract.from_strings(
        input_vars=["p_i","p_j","q_i","q_j"],
        output_vars=["m_j"],
        assumptions=["q_j >= p_j"],
        guarantees=[f"{delta}*m_j >= {1/(1 + b_j*(A_j_H + A_i_L))}*p_j + {A_j_H*b_j/(1 + b_j*(A_j_H + A_i_L))}*q_j"]
    )

    monotonicity_component2_sc2 = PolyhedralIoContract.from_strings(
        input_vars=["p_i","p_j","q_i","q_j"],
        output_vars=["m_i"],
        assumptions=[f"{1 + b_i*A_j_H} * q_i >= p_i"],
        guarantees=[f"{delta}*m_i <= {1/(1 + b_i*(A_j_H + A_i_L))}*p_i + {A_i_L*b_i/(1 + b_i*(A_j_H + A_i_L))}*q_i"]
    )

    A_definition_sc2 = PolyhedralIoContract.from_strings(
        input_vars=["gamma_inv"],
        output_vars=["A_i", "A_j"],
        assumptions=[],
        guarantees=[f"{kappa * m_j_H} * gamma_inv <= A_j",
                    f"A_i <= {kappa * m_i_L}  * gamma_inv"]
    )

    # compute the parameters that yield the first set of top-level requirements
    parameters_contract_2 = top_level_sc2.quotient(monotonicity_component1_sc2.compose(monotonicity_component2_sc2).compose(A_definition_sc2))

    # report the parameter space that yields bistability
    parameter_contract = parameters_contract_1.merge(parameters_contract_2)
    # print(parameter_contract)


    # These parameters guarantee that the system will have the desired properties. We observe that $\gamma$ is not bounded by other parameters.
    # Similarly, there are inequalities that relate $p_i$ with $q_i$ and $p_j$ with $q_j$. Thus, the following are the valid parameters:

    # In[9]:


    from pacti.iocontract import Var
    from pacti.terms.polyhedra import PolyhedralTermList
    from pacti.utils.plots import plot_guarantees

    # gamma
    gamma_bounds = parameter_contract.get_variable_bounds("gamma_inv")
    # print(f"Gamma has the bounds {[1/e for e in reversed(list(gamma_bounds))]}")

    # q_i and p_i
    contract_i = PolyhedralIoContract(
        assumptions=PolyhedralTermList(),
        guarantees=parameter_contract.g.get_terms_with_vars([Var('p_i'), Var('q_i')]),
        input_vars=[],
        output_vars=[Var('p_i'), Var('q_i')]
        )

    # print(contract_i)

    _ = plot_guarantees(contract=contract_i,x_var="p_i",y_var="q_i",var_values={},x_lims=[-5,0.05],y_lims=[0,8],show=True)

    # q_j and p_j
    contract_j = PolyhedralIoContract(
        assumptions=PolyhedralTermList(),
        guarantees=parameter_contract.g.get_terms_with_vars([Var('p_j'), Var('q_j')]),
        input_vars=[],
        output_vars=[Var('p_j'), Var('q_j')]
        )

    # print(contract_j)

    import matplotlib.pyplot as plt
    _ = plot_guarantees(contract=contract_j,x_var="p_j",y_var="q_j",var_values={},x_lims=[-5,0.05],y_lims=[0,8],show=True)
    plt.show()

    # In[ ]:




