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

import pandas as pd
from nose.tools import eq_, ok_

from cohorts.plot import mann_whitney_plot

def test_mann_whitney():
    data = pd.DataFrame({"distribution": range(1, 7),
                         "condition": [True, False] * 3})
    U_default, p_default, _, _, _, _ = mann_whitney_plot(
        data, "condition", "distribution", skip_plot=True)
    U_two, p_two, _, _, _, _ = mann_whitney_plot(
        data, "condition", "distribution", alternative="two-sided", skip_plot=True)
    U_less, p_less, _, _, _, _ = mann_whitney_plot(
        data, "condition", "distribution", alternative="less", skip_plot=True)
    U_greater, p_greater, _, _, _, _ = mann_whitney_plot(
        data, "condition", "distribution", alternative="greater", skip_plot=True)
    eq_(p_default, p_two)
    ok_((p_default == p_less * 2) or (p_default == p_greater * 2))
