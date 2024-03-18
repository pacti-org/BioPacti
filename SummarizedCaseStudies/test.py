#  This file uncovered a problem in Pacti: the quotient fails when the user installs SciPy > 1.10


from pacti.contracts import PolyhedralIoContract

contract1 = PolyhedralIoContract.from_strings(
    input_vars=[],
    output_vars=["m_i", "m_j"],
    assumptions=[],
    guarantees=[  "-m_i <= -5", "m_j <= 1"])

contract2 = PolyhedralIoContract.from_strings(
    input_vars=["p_i", "p_j", "q_i", "q_j"],
    output_vars=["m_i", "m_j"],
    assumptions=[],
    guarantees=["-m_i + 0.007634 p_i + 0.229 q_i <= 0",
                "m_j - 0.006211 p_j - 0.6211 q_j <= 0"]
)

cq = contract1.quotient(contract2)

print(f"** Contract1 is \n{contract1}")
print(f"** Contract2 is \n{contract2}")
print(f"** C1/C2 is \n{cq}")

