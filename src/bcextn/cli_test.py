import sys

def main():

    if '--updateref' in sys.argv[1:]:
        from .test_eq36 import update_ref
        update_ref()
        return

    if '--benchmark' in sys.argv[1:]:
        from .eval_eq36 import benchmark
        benchmark()
        return

    from .eval_fofeta import test_f_of_eta
    test_f_of_eta()
    from .eval_phi0 import test_phi0
    test_phi0()
    from .eval_phipi import test_phipi
    test_phipi()
    from .test_eq36 import test
    test()
    print("All ok")
