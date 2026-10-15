"""raslite — 一维恒定渐变流水面线推算（对标 HEC-RAS 1D Steady Flow）

模块:
    geometry    横断面几何与水力属性 (A/P/T/K/alpha)
    hydraulics  临界水深、正常水深、Froude 数
    profile     标准步进法 (Standard Step Method) 水面线推算

参考: HEC-RAS Hydraulic Reference Manual v6.0, Chapter 2
单位制: 公制 SI (m, m3/s)
"""

G = 9.80665  # 重力加速度 m/s2

from .geometry import CrossSection
from .hydraulics import critical_ws, normal_ws, froude_squared
from .profile import SteadyProfile, SectionResult, specific_force

__all__ = [
    "CrossSection", "critical_ws", "normal_ws", "froude_squared",
    "SteadyProfile", "SectionResult", "specific_force", "G",
]
