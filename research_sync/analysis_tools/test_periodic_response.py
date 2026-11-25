import unittest
import numpy as np
from analyze_periodic_response import harmonic,metrics,period

class PeriodicChecks(unittest.TestCase):
    def test_jittered_exact_multitone_with_drift(self):
        rng=np.random.default_rng(42)
        t=np.arange(10,30,.01)+rng.uniform(-.001,.001,2000)
        f=np.array([.05,.1,.2,.3,.45,.6,.7,.9,1.2,1.5])
        c=np.arange(1,11)/10+1j*np.arange(10,0,-1)/20
        y=3+.04*t+np.real(np.exp(2j*np.pi*t[:,None]*f)@c)
        np.testing.assert_allclose(harmonic(t,y,f),c,atol=1e-12)
        self.assertEqual(period(f),20.)

    def test_complex_aggregation_and_noise_definition(self):
        U=np.ones((4,1),complex)
        Y=(2+np.array([1,-1,1j,-1j]))[:,None]
        H,Z,snr,adj,coh,p,q,sd=metrics(U,Y)
        np.testing.assert_allclose(H,2)
        np.testing.assert_allclose(coh,.8)
        np.testing.assert_allclose(snr,10*np.log10(3))
        np.testing.assert_allclose(p,.2**3)

    def test_phase_wrap_and_one_segment(self):
        U=np.ones((4,1),complex)
        Y=np.exp(1j*np.deg2rad([179,-179,178,-178]))[:,None]
        H,*_=metrics(U,Y)
        self.assertGreater(abs(np.angle(H[0],deg=True)),179)
        with np.errstate(invalid='ignore',divide='ignore'):
            self.assertTrue(np.isnan(metrics(U[:1],Y[:1])[4]).all())

if __name__=='__main__':unittest.main()
