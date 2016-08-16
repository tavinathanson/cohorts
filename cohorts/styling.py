# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function

import matplotlib as mpl
import matplotlib.colors as colors
import seaborn as sb
import numpy as np

def set_styling():
    sb.set_style("white")
    red = colors.hex2color("#bb3f3f")
    blue = colors.hex2color("#5a86ad")
    deep_colors = sb.color_palette("deep")
    green = deep_colors[1]
    custom_palette = [red, blue, green]
    custom_palette.extend(deep_colors[3:])
    sb.set_palette(custom_palette)
    mpl.rcParams.update({"figure.figsize": np.array([6, 6]),
                         "legend.fontsize": 12,
                         "font.size": 16,
                         "axes.labelsize": 16,
                         "axes.labelweight": "bold",
                         "xtick.labelsize": 16,
                         "ytick.labelsize": 16})
