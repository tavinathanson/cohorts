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

from nose.tools import eq_

from cohorts.rounding import float_str

ROUNDING_CASES = {2: '2.00',
                  2.0: '2.00',
                  2.011: '2.01',
                  4.015: '4.02',
                  4.025: '4.03',
                  3.015: '3.02',
                  3.025: '3.03',
                  0.001: '0.001',
                  0.0010: '0.001',
                  0.00104: '0.0010',
                  0.00109: '0.0011',
                  0.00111: '0.0011',
                  4.5: '4.50',
                  1235345.132424: '1235345.13',
                  0.009: '0.009',
                  234.009: '234.01',
                  3.000001: '3.00',
                  13212324.0: '13212324.00',
                  -1.1: '-1.10'}

def test_rounding():
    for f, expected_str in ROUNDING_CASES.items():
        eq_(expected_str, float_str(f))
