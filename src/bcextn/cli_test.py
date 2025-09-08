import sys

def main():
    if '--updateref' in sys.argv[1:]:
        from .test_eq36 import update_ref
        update_ref()
        return

    if '--updateref2' in sys.argv[1:]:
        from .test_eq36 import gen_test2_data
        gen_test2_data()
        return

    if '--benchmark' in sys.argv[1:]:
        from .eval_eq36 import benchmark
        benchmark()
        return


    from .eval_fofeta import test_f_of_eta
    test_f_of_eta()
    from .eval_fofeta_scnd_fresnel import test_f_of_eta_fresnel
    test_f_of_eta_fresnel()
    from .eval_eq36_scnd import test as test_eq36_scnd
    test_eq36_scnd()
    from .eval_phi0 import test_phi0
    test_phi0()
    from .eval_phipi import test_phipi
    test_phipi()
    from .test_eq36 import test, test2
    test()
    test2()
    print("All ok")
