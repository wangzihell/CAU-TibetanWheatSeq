#!/usr/bin/env python3

import dadi
from optparse import OptionParser

usage = "Usage: %prog [-i <input file>] [-o <output file>] [-s <column number of first sample>]\n"
#
parser = OptionParser(usage)
parser.add_option("-i", dest="infile",
              help="Input file, use STDIN if omit; "
                   "gzip file end with \".gz\".", default = "test")
parser.add_option("-n", dest="randn",
              help="Input file, use STDIN if omit; "
                   "gzip file end with \".gz\".", default = "test")
parser.add_option("-f", dest="wfold",
              help="Input file, use STDIN if omit; "
                   "gzip file end with \".gz\".", default = "test")
#
(options, args) = parser.parse_args()

def sec_contact_sym_mig_size(params, ns, pts):
    """
    Split with no gene flow, followed by size change with symmetrical gene flow.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and present.
    m: Migration between pop 2 and pop 1.
    """
    nu2a, nu2b, T1, T2, m = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, 1, nu2a)

    phi = dadi.Integration.two_pops(phi, xx, T2, 1, nu2b, m12=m, m21=m)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

if __name__ == '__main__':
    # Read in data
    if options.wfold is "no":
        data_un = dadi.Spectrum.from_file(options.infile)
        data = data_un.fold()
    else:
        data = dadi.Spectrum.from_file(options.infile)
    ns = data.sample_sizes
    pts_l = [40,50,60]
    func = sec_contact_sym_mig_size
    upper_bound = [10.0,10.0,10.0,10.0,10.0]
    lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4]
    p0 = [0.1,0.1,0.1,0.1,0.1]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    p0 = dadi.Misc.perturb_params(p0,fold=1,upper_bound=upper_bound,lower_bound=lower_bound)

    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0))
    model = func_ex(popt, ns, pts_l)
    ll_model = dadi.Inference.ll_multinom(model, data)
    theta = dadi.Inference.optimal_sfs_scaling(model, data)

    with open("3epoch_m_"+options.randn+".csv", 'a') as f:
        for p in popt:
            print("{:.5f},".format(p), sep='',end='', file=f)
        print(ll_model,theta, sep=",", file=f)
