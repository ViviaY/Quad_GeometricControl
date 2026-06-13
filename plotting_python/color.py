import numpy as np
import matplotlib.pyplot as plt


class ColorsData:
    def __init__(self):
        self.colors_v1 = np.array([
            [0, 0.4470, 0.7410],  # Blue (Default MATLAB Blue)
            [0.8500, 0.3250, 0.0980],  # Orange (MATLAB Default)
            [0.9290, 0.6940, 0.1250],  # Yellow
            [0.4940, 0.1840, 0.5560],  # Purple
            [0.4660, 0.6740, 0.1880],  # Green
            [0.3010, 0.7450, 0.9330],  # Light Blue
            [0, 0, 0],  # Black
            [0.75, 0.75, 0.75],  # Gray
            [1, 0, 1],  # Magenta
            [0, 1, 1],  # Cyan
            [0, 0.5, 0],  # Dark Green (Olive)
            [0, 0.4470, 0.7410],  # Deep Blue (Repeated to ensure visibility)
            [0.5, 0, 0.5],  # Dark Purple
            [0.3, 0.3, 0.3],  # Dark Gray
            [0.6350, 0.0780, 0.1840],  # Dark Brown
            [0.2, 0.5, 0.3],  # Forest Green
            [0.2, 0.2, 0.6],  # Dark Navy Blue
            [0.3, 0.6, 0.8],  # Sky Blue
            [0.8, 0.4, 0.0],  # Burnt Orange
            [0.8, 0.6, 0.7]  # Soft Pink
        ])
        self.colors = {
            # bounds    
            "L": "#D62728",

            # GC
            "GC_1": "#1F77B4",
            "GC_2": "#547FB6",
            "GC_3": "#0B3C5D",

            # NMPC
            "NMPC10": "#76D7C4",
            "NMPC15": "#1ABC9C",
            "NMPC20": "#008B74",
            "NMPC25": "#005B4F",

            # CBF 
            "CBF_epi": "#C97DD8",
            "CBF_safety": "#9B4F96",
            "CBF_reach": "#5E2A84",
        }
