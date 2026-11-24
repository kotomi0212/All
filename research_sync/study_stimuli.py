"""Runnable 2026-11-24 visual profiles for the acquisition application."""
import copy
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent/'study_profiles'
LABELS={
    '90秒 スイープ 0.70→0.05 Hz（中心固定）':'study_sweep90',
    '100秒 マルチコサイン10周波数（20秒×4周期）':'study_multicosine10',
    '10周波数 100秒・速度成分一定【最大加速度2.09 m/s²】':'study_multicosine10_equal_velocity',
}

def profile(condition):
    if condition=='study_multicosine10_equal_velocity':
        from .low_band import BASE_A, LIMITS, check_limits
        c=json.loads((ROOT/'study_multicosine10.json').read_text(encoding='utf-8'))
        original=copy.deepcopy(c);original['acceleration_component_m_s2']=4*BASE_A
        target=check_limits(original)[0]['displacement_m']['value']/2
        c.pop('acceleration_component_m_s2',None)
        omega=2*np.pi*np.asarray(c['frequencies_hz'])
        c['acceleration_components_m_s2']=omega.tolist()
        unit=check_limits(c)[0]['displacement_m']['value']
        speed=target/unit
        c['acceleration_components_m_s2']=(speed*omega).tolist()
        c.update(condition=condition,enabled=True,velocity_component_m_s=speed,
            target_max_displacement_m=target,engineering_limits=LIMITS,
            normalization='Equal per-component velocity; maximum displacement equals half of previous 10-frequency 4x stimulus',
            status='runnable_user_authorized',
            limit_override='User authorized acquisition with calculated acceleration peak 2.09 m/s2',
            caution='Provisional engineering limits may be exceeded; human safety not established')
        c['maxima'],c['limit_violations']=check_limits(c)
        return c
    if condition not in LABELS.values():raise ValueError('Unknown study stimulus')
    c=json.loads((ROOT/(condition+'.json')).read_text(encoding='utf-8'))
    if c.get('condition')!=condition or not c.get('enabled'):raise ValueError('Stimulus is not enabled')
    if c['waveform']=='sampled_centered_sweep':
        samples=np.asarray(c['samples_xva'],float)
        if samples.shape!=(4501,3) or not np.isfinite(samples).all() or c['sample_step_s']!=.02 or c['duration_s']!=90:
            raise ValueError('Invalid sweep samples')
        peaks=np.max(abs(samples),axis=0)
        if np.any(peaks>np.array([.30,.10,.20])):raise ValueError('Sweep exceeds configured bounds')
        if max(abs(samples[[0,-1],1]))>1e-8:raise ValueError('Sweep endpoint velocity is not zero')
    elif c['waveform']=='pilot_acceleration_multicosine':
        from .low_band import check_limits
        _,violations=check_limits(c)
        if violations:raise ValueError('Stimulus exceeds configured bounds: '+', '.join(violations))
    else:raise ValueError('Unsupported study waveform')
    return copy.deepcopy(c)

def available_conditions():
    """Only expose pilot conditions whose generated profiles are enabled."""
    result=dict(LABELS)
    result.update({'90秒 マルチコサイン30周波数（従来）':'multicosine',
        '30秒 注視点なし':'none','30秒 赤い固定点':'dot','30秒 ローディング表示':'spinner',
        '静止60秒：ノイズ確認':'pilot_rest'})
    for kind,title,T in [('single','単一0.20 Hz',40),('five','5周波数',100)]:
        for n in ((1,1.4,2,3) if kind=='five' else (1,2,3)):
            key=f'pilot_{kind}_x{str(n).replace(chr(46),chr(112))}'
            try:c=json.loads((ROOT.parent/'pilot_profiles'/(key+'.json')).read_text(encoding='utf-8'))
            except (OSError,ValueError):continue
            if c.get('enabled'):result[f'{title} {T}秒・成分振幅{n}倍']=key
    return result
