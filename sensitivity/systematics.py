def detector_systematics(sig_or_bkg, mass,detvars="detvars"):
    # detector variations only make sense compared to CV
    # so need to do everything in relation to that
    detvars_list = [] # return this
    if detvars == "detvars":
        detvars_names = [
            "CV",
            "wiremod_x",
            "wiremod_yz",
            "wiremod_anglexz",
            "wiremod_angleyz",
            "wiremod_dEdx",
            "ly_down",
            "ly_rayleigh",
            "ly_atten",
            "sce",
            "recomb"
        ]
    else:
        detvars_names = ["quadrature"]

    detvars_dict = {}
    
    ratios_root_file = 'systematics/detvar_def-th_1000kev_2hits_%imev_v08_00_00_61/Quadrature_%s.root'%(mass,sig_or_bkg)

    ur = uproot.open(ratios_root_file)
    for dv in detvars_names:
        try:
            detvars_list.append(ur[dv].values().tolist())
            detvars_dict[dv] = ur[dv].values().tolist()
        except:
            print("No histogram found for detvar",dv)
            detvars_dict[dv] = []

    # returns dictionary of samples : histograms (list of bin values)
    return detvars_dict

def detector_systematics_sigbkg(mass):
    # calls previous, returns both
    sig_syst = detector_systematics("signal",mass)
    bkg_syst = detector_systematics("background",mass)
    return sig_syst,bkg_syst

def detector_systematics_quadrature(mass):
    # calls previous, returns both
    sig_syst = detector_systematics("signal",mass,detvars="quadrature")
    bkg_syst = detector_systematics("background",mass,detvars="quadrature")
    return sig_syst,bkg_syst
