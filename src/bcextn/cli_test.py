def main():
    from .eval_fofeta import test_f_of_eta
    test_f_of_eta()
    from .eval_phi0 import test_phi0
    test_phi0()
    from .eval_phipi import test_phipi
    test_phipi()
    print("All ok")
