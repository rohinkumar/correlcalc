class cosmology:
    def __init__(self):
        print "E(z) and DC (comoving distance) methods"
        #print "Cosmology models for Comoving distance calculation class initiated"

    def Ez(self,zv,Om,Ol,Ok):
    	return 1.0/m.sqrt(Om*(1.0+zv)**3+Ok*(1.0+zv)**2+Ol)

    def DC_LCDM(self,z):
	np.vectorize(Ez)
        return integrate.quad(Ez, 0, z)[0]
