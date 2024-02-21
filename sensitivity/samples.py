from math import sqrt

class Sample:
    data_pot = 1.5e21
    data_eve = 3081362.
    _sig_or_bkg        = "sig"
    _id                = "0"
    _bdt_hist          = ""
    _bdt_hist_unscaled = ""
    _description       = "dummy signal"
    _scaling_type      = "pot" # pot or evt
    _scaling_factor    = 1.0
    _has_systematics   = False
    _syst_hist         = ""
    _err_hist          = ""
    
    def __init__(self,uprooted, samplename="sig0"):
        if "sig" in samplename:
            self._sig_or_bkg = "sig"
        elif "bkg" in samplename:
            self._sig_or_bkg = "bkg"
        self._id       = samplename.replace(self._sig_or_bkg,"")
        self._bdt_hist = uprooted[f'{self._sig_or_bkg}{self._id}_hist']
        self._bdt_hist = self._bdt_hist.values().tolist()
        for k in uprooted['total_pot'].keys():
            if samplename in k:
                self._scaling_type = k.replace(f"_{samplename}","")
                value = float(uprooted['total_pot'][k].array(library="np")[0])
                self._scaling_factor = 1.0/value
        if self._scaling_type == "pot":
            self._scaling_factor *= self.data_pot
        elif self._scaling_type == "evt":
            self._scaling_factor *= self.data_eve
        self._bdt_hist_unscaled = self._bdt_hist
        if self._has_systematics == False:
            self._err_hist = [ sqrt(i) for i in self._bdt_hist ]
        # scale everything now
        
        self._bdt_hist = [ self._scaling_factor * i for i in self._bdt_hist ]
        self._err_hist = [ self._scaling_factor * i for i in self._err_hist ]
        #self._bdt_hist = temp_bdt
        #self._err_hist = temp_err
        
            
    def add_systematics(self,sys_hist):
        self._syst_hist = syst_hist
        self._has_systematics = True
        
    def get_systematics(self):
        return self._syst_hist

    def get_scaling(self):
        return self._scaling_factor

    def get_hist_values(self):
        return self._bdt_hist

    def get_err_values(self):
        return self._err_hist
