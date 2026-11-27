"""Verify selectable profiles, renderer parity and measurement UI catalogue."""
import json
from pathlib import Path
import subprocess
import unittest
import numpy as np
from .study_stimuli import ROOT,LABELS,profile,available_conditions
from .trial_stimuli import select_stimulus
from .named_recording import STIMULI
from .low_band import evaluate

class StudyStimuliTests(unittest.TestCase):
    def test_band_sweep_profiles_and_independent_copies(self):
        from .band_sweeps import PROFILES, EDGES
        self.assertEqual(len(PROFILES),6)
        for name in PROFILES:
            c=select_stimulus({},name)
            self.assertEqual(c['duration_s'],100)
            self.assertEqual(c['analysis_interval_s'],[10,95])
            self.assertEqual(c['onset_s']+c['sweep_s']+c['offset_s'],100)
            self.assertIn(name,available_conditions().values())
            self.assertIn(name,STIMULI)
            self.assertEqual(c,json.loads((ROOT/f'{name}.json').read_text(encoding='utf-8')))
            c['analysis_interval_s'][0]=0
            self.assertEqual(profile(name)['analysis_interval_s'],[10,95])
        self.assertAlmostEqual(EDGES[-1],.5)
        self.assertAlmostEqual(EDGES[1]/EDGES[0],EDGES[3]/EDGES[2])
        for old in ('study_sweep90','study_sweep300'):
            self.assertNotIn(old,available_conditions().values())
            with self.assertRaises(ValueError):select_stimulus({},old)
    def test_new_profiles_and_enabled_catalogue(self):
        conditions=available_conditions()
        self.assertEqual(list(conditions.values())[:len(LABELS)],list(LABELS.values()))
        self.assertNotIn('pilot_five_x2',conditions.values())
        self.assertNotIn('pilot_five_x3',conditions.values())
        for condition in LABELS.values():
            c=select_stimulus({},condition)
            self.assertTrue(c['enabled']);self.assertIn(condition,STIMULI)
        for condition in conditions.values():
            if condition.startswith('pilot_'):
                c=json.loads((ROOT.parent/'pilot_profiles'/(condition+'.json')).read_text())
                self.assertTrue(c['enabled'])

    def test_equal_velocity_profile_and_blocked_conditions(self):
        from .low_band import BASE_A, profile as pilot_profile
        condition='study_multicosine10_equal_velocity'
        c=profile(condition)
        self.assertTrue(c['enabled'])
        self.assertIn(condition,available_conditions().values())
        self.assertTrue(select_stimulus({},condition)['enabled'])
        self.assertEqual(len(c['frequencies_hz']),10)
        np.testing.assert_allclose(np.array(c['acceleration_components_m_s2'])/(2*np.pi*np.array(c['frequencies_hz'])),c['velocity_component_m_s'])
        self.assertAlmostEqual(c['maxima']['displacement_m']['value'],c['target_max_displacement_m'],places=7)
        self.assertIn(condition,STIMULI)
        for old in ('study_multicosine10_x4','study_five_x4'):
            self.assertNotIn(old,available_conditions().values())
            with self.assertRaises(ValueError):select_stimulus({},old)
        times=np.linspace(0,100,101)
        code="const fs=require('fs'),w=require('./research_sync/waveform.js'),q=JSON.parse(fs.readFileSync(0,'utf8'));console.log(JSON.stringify(q.t.map(t=>w.evaluate(t,q.c))));"
        result=subprocess.run(['node','-e',code],input=json.dumps(dict(c=c,t=times.tolist())),text=True,capture_output=True,check=True,cwd=ROOT.parents[1])
        js=np.array([[r[k] for k in ('x','v','a')] for r in json.loads(result.stdout)])
        np.testing.assert_allclose(js,np.array(evaluate(times,c)).T,atol=1e-10)
        for condition in ('pilot_five_x2','pilot_five_x3'):
            with self.assertRaises(ValueError):pilot_profile(condition)

    def test_browser_sweep_matches_independent_design_samples(self):
        c=json.loads((ROOT/'study_sweep90.json').read_text(encoding='utf-8'))
        reference=np.array(json.loads((ROOT/'verification_samples.json').read_text()))
        code="const fs=require('fs'),w=require('./research_sync/waveform.js'),q=JSON.parse(fs.readFileSync(0,'utf8'));console.log(JSON.stringify(q.t.map(t=>w.evaluate(t,q.c))));"
        result=subprocess.run(['node','-e',code],input=json.dumps(dict(c=c,t=reference[:,0].tolist())),text=True,capture_output=True,check=True,cwd=ROOT.parents[1])
        actual=np.array([[r[k] for k in ('x','v','a')] for r in json.loads(result.stdout)])
        np.testing.assert_allclose(actual,reference[:,1:],atol=2e-8,rtol=0)

    def test_multicosine_lines_and_browser_derivatives(self):
        c=profile('study_multicosine10');t=np.arange(10,90,.01)
        x,v,a=evaluate(t,c);bins=np.fft.rfftfreq(len(t),.01);s=2*abs(np.fft.rfft(a))/len(t)
        amps=s[np.rint(np.array(c['frequencies_hz'])*80).astype(int)]
        np.testing.assert_allclose(amps,c['acceleration_component_m_s2'],rtol=1e-10)
        self.assertEqual(len(c['frequencies_hz']),10)
        np.testing.assert_allclose(evaluate([0,100],c),0,atol=1e-12)
        code="const fs=require('fs'),w=require('./research_sync/waveform.js'),c=JSON.parse(fs.readFileSync(0,'utf8'));let worst=[0,0];for(let t=.013;t<89.99;t+=.137){const h=1e-5,l=w.evaluate(t-h,c),r=w.evaluate(t+h,c),s=w.evaluate(t,c);worst[0]=Math.max(worst[0],Math.abs((r.x-l.x)/(2*h)-s.v));worst[1]=Math.max(worst[1],Math.abs((r.v-l.v)/(2*h)-s.a));}console.log(JSON.stringify(worst));"
        result=subprocess.run(['node','-e',code],input=json.dumps(json.loads((ROOT/'study_sweep90.json').read_text(encoding='utf-8'))),text=True,capture_output=True,check=True,cwd=ROOT.parents[1])
        self.assertLess(max(json.loads(result.stdout)),1e-7)

if __name__=='__main__':unittest.main()
